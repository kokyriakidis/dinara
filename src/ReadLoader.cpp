// Dinara.
#include "ReadLoader.hpp"
#include "computeRunLengthRepresentation.hpp"
#include "filesystem.hpp"
#include "performanceLog.hpp"
#include "splitRange.hpp"
using namespace dinara;

// Standard library.
#include "chrono.hpp"
#include "iterator.hpp"
#include <filesystem>
#include "tuple.hpp"

// zlib for gzip decompression.
#include <zlib.h>

#include "MultithreadedObject.tpp"
namespace dinara {
template class MultithreadedObject<ReadLoader>;
}

// Load reads from a fastq or fasta file.
ReadLoader::ReadLoader(
    const string& fileName,
    uint64_t representation, // 0 = raw sequence, 1 = RLE sequence
    uint64_t minReadLength,
    bool noCache,
    uint64_t threadCount,
    const string& dataNamePrefix,
    size_t pageSize,
    Reads& reads):

    MultithreadedObject(*this),
    fileName(fileName),
    representation(representation),
    minReadLength(minReadLength),
    noCache(noCache),
    threadCount(threadCount),
    dataNamePrefix(dataNamePrefix),
    pageSize(pageSize),
    reads(reads)
{
    performanceLog << timestamp << "Loading reads from " << fileName << endl;

    adjustThreadCount();

    // Get the file extension.
    string extension;
    try {
        extension = filesystem::extension(fileName);
    } catch (...) {
        throw runtime_error("Input file " + fileName +
            " must have an extension consistent with its format.");
    }

    // Handle gzip-compressed files: strip .gz and record that we need decompression.
    if(extension == "gz" || extension == "GZ") {
        isGzipped = true;
        // Get the inner extension (e.g., "fastq" from "reads.fastq.gz").
        const string nameWithoutGz = fileName.substr(0, fileName.size() - 3);
        try {
            extension = filesystem::extension(nameWithoutGz);
        } catch (...) {
            throw runtime_error("Input file " + fileName +
                " must have a format extension before .gz (e.g., .fastq.gz, .fasta.gz).");
        }
    }

    // Fasta file. ReadLoader is more forgiving than OldFastaReadLoader.
    if(extension=="fasta" || extension=="fa" || extension=="FASTA" || extension=="FA") {
        processFastaFile();
        return;
    }

    // Fastq file.
    if(extension=="fastq" || extension=="fq" || extension=="FASTQ" || extension=="FQ") {
        processFastqFile();
        return;
    }

    // If getting here, the file extension is not supported.
    throw runtime_error("File extension " + extension + " is not supported. "
        "Supported file extensions are .fasta, .fa, .FASTA, .FA, .fastq, .fq, .FASTQ, .FQ "
        "(optionally with .gz suffix).");
}


// This is required because unique_ptr has a type completeness requirement
// during destruction.
ReadLoader::~ReadLoader() = default;


void ReadLoader::adjustThreadCount()
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
}



void ReadLoader::processFastaFile()
{

    // Read the entire fasta file.
    const auto t0 = std::chrono::steady_clock::now();
    allocateBufferAndReadFile();

    // Each thread stores reads in its own data structures.
    const auto t1 = std::chrono::steady_clock::now();
    allocatePerThreadDataStructures();
    runThreads(&ReadLoader::processFastaFileThreadFunction, threadCount);
    const auto t2 = std::chrono::steady_clock::now();
    buffer.remove();

    // Store the reads computed by each thread and free
    // the per-thread data structures.
    storeReads();
    const auto t3 = std::chrono::steady_clock::now();


    performanceLog << "Time to process this file:\n" <<
        "Allocate buffer + read: " << seconds(t1-t0) << " s.\n" <<
        "Parse: " << seconds(t2-t1) << " s.\n"
        "Store: " << seconds(t3-t2) << " s.\n"
        "Total: " << seconds(t3-t0) << " s." << endl;


}



// Each thread processes the reads whose initial ">" character
// is in the file block assigned to the read.
void ReadLoader::processFastaFileThreadFunction(size_t threadId)
{
    const char* bufferPointer = &buffer[0];
    const uint64_t bufferSize = buffer.size();

    // Allocate and access the data structures where this thread will store the
    // reads it finds.
    allocatePerThreadDataStructures(threadId);
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames = *threadReadNames[threadId];
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadMetaData = *threadReadMetaData[threadId];
    LongBaseSequences& thisThreadReads = *threadReads[threadId];
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
        *threadReadRepeatCounts[threadId];


    // Compute the file block assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, bufferSize, threadCount, threadId);
    if(begin == end) {
        return;
    }

    // Locate the first read that begins in this block.
    uint64_t offset = begin;
    if(offset == 0) {
        if(!fastaReadBeginsHere(offset)) {
            throw runtime_error("Fasta file " + fileName + " does not begin with a \">\".");
        }
    } else {
        while(!fastaReadBeginsHere(offset)) {
            ++offset;
            if(offset == end) {
                // We reached the end of the block assigned to this thread
                // without finding any reads.
                return;
            }
        }
    }

    // Ok, we are at the ">" for the first read assigned to this thread.

    // Main loop over the reads in the file block allocated to this thread.
    string readName;
    string readMetaData;
    vector<Base> read;
    vector<Base> runLengthRead;
    vector<uint8_t> readRepeatCount;
    while(offset < end) {
        DINARA_ASSERT(fastaReadBeginsHere(offset));

        // Consume the ">".
        ++offset;

        // Read the name.
        readName.clear();
        while(true) {
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing a read name.");
            }
            const char c = bufferPointer[offset];
            if(isspace(c)) {
                break;
            } else {
                readName.push_back(c);
            }
            offset++;
        }

        // Read the meta data.
        readMetaData.clear();
        while(true) {
            if(offset == bufferSize) {
                throw runtime_error("Reached end of file while processing header line for read "
                    + readName);
            }
            const char c = bufferPointer[offset];
            if(c == '\n') {
                break;
            }
            ++offset;
            if(isspace(c) and readMetaData.empty()) {
                // Do nothing. Only start storing after encountering the first non-space.
            } else {
                readMetaData.push_back(c);
            }
        }

        // Consume the line end at the end of the header line.
        ++offset;



        // Read the bases.
        // Note that here we can go past the file block assigned to this thread.
        read.clear();
        uint64_t invalidBaseCount = 0;
        while(offset != bufferSize) {
            const char c = bufferPointer[offset];

            // Skip white space.
            if(c==' ' || c=='\t' || c=='\n' || c=='\r') {
                ++offset;
                continue;
            }

            // If we reached the beginning of another read, stop here.
            if(c=='>' && fastaReadBeginsHere(offset))  {
                break;
            }

            // Get a base from this character.
            const Base base = Base::fromCharacterNoException(c);
            if(base.isValid()) {
                read.push_back(base);
            } else {
                ++invalidBaseCount;
            }
            ++offset;
        }


        // If we found invalid bases, skip this read.
        if(invalidBaseCount) {
            __sync_fetch_and_add(&discardedInvalidBaseReadCount, 1);
            __sync_fetch_and_add(&discardedInvalidBaseReadCount, read.size());
            continue;
        }

        // If the read is too short, skip it.
        if(read.size() < minReadLength) {
            __sync_fetch_and_add(&discardedShortReadReadCount, 1);
            __sync_fetch_and_add(&discardedShortReadBaseCount, read.size());
            continue;
        }

        // Store the read bases.
        if(representation == 1) {
            if(computeRunLengthRepresentation(read, runLengthRead, readRepeatCount)) {
                thisThreadReadNames.appendVector(readName.begin(), readName.end());
                thisThreadReadMetaData.appendVector(readMetaData.begin(), readMetaData.end());
                thisThreadReads.append(runLengthRead);
                thisThreadReadRepeatCounts.appendVector(readRepeatCount);
            } else {
                __sync_fetch_and_add(&discardedBadRepeatCountReadCount, 1);
                __sync_fetch_and_add(&discardedBadRepeatCountBaseCount, read.size());
            }
        } else {
            thisThreadReadNames.appendVector(readName.begin(), readName.end());
            thisThreadReadMetaData.appendVector(readMetaData.begin(), readMetaData.end());
            thisThreadReads.append(read);
        }
    }

}



// Function that returns true if a read begins
// at this position in Fasta format.
bool ReadLoader::fastaReadBeginsHere(uint64_t offset) const
{
    if(buffer[offset] == '>') {
        if(offset == 0) {
            return true;
        } else {
            return buffer[offset-1] == '\n';
        }
    } else {
        return false;
    }
}



// Process a fastq file under the following assumptions:
// - 4 lines per read.
// - For each read:
//   * The first line must begin with '@'.
//   * The third line must consist of just a '+' with nothing else.
//   * The second and fourth line must have the same number of bytes,
//     also equal to the number of bases in the read. The fourth
//     line is otherwise ignored.
//   * The second line can only contain ACGT. Invalid bases are not allowed.
// - No Windows line ends.
void ReadLoader::processFastqFile()
{

    // Read the entire fastq file.
    const auto t0 = std::chrono::steady_clock::now();
    allocateBufferAndReadFile();

    // Find all line ends in the file.
    const auto t1 = std::chrono::steady_clock::now();
    findLineEnds();
    // cout << "Found " << lineEnds.size() << " lines in this file." << endl;

    // Check that the number of lines is a multiple of 4
    // (there must be exactly 4 mlines per read per the above assumptions).
    if((lineEnds.size() %4) != 0) {
        throw runtime_error("File has " + to_string(lineEnds.size()) +
            " lines. Expected a multiple of 4. "
            "Only fastq files with each read on exactly 4 lines are supported.");
    }

    // Each thread stores reads in its own data structures.
    const auto t2 = std::chrono::steady_clock::now();
    allocatePerThreadDataStructures();
    runThreads(&ReadLoader::processFastqFileThreadFunction, threadCount);
    buffer.remove();


    // Store the reads computed by each thread and free
    // the per-thread data structures.
    const auto t3 = std::chrono::steady_clock::now();
    storeReads();
    const auto t4 = std::chrono::steady_clock::now();


    performanceLog << "Time to process this file:\n" <<
        "Allocate buffer + read: " << seconds(t1-t0) << " s.\n" <<
        "Locate: " << seconds(t2-t1) << " s.\n"
        "Parse: " << seconds(t3-t2) << " s.\n"
        "Store: " << seconds(t4-t3) << " s.\n"
        "Total: " << seconds(t4-t0) << " s." << endl;


}



void ReadLoader::processFastqFileThreadFunction(size_t threadId)
{

    // Allocate and access the data structures where this thread will store the
    // reads it finds.
    allocatePerThreadDataStructures(threadId);
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames = *threadReadNames[threadId];
    MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadMetaData = *threadReadMetaData[threadId];
    LongBaseSequences& thisThreadReads = *threadReads[threadId];
    MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
        *threadReadRepeatCounts[threadId];

    // Find the total number of reads in the file.
    DINARA_ASSERT((lineEnds.size() % 4) == 0); // We already checked for that.
    const ReadId readCountInFile = ReadId(lineEnds.size() / 4);

    // Compute the range of reads in the file assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, readCountInFile, threadCount, threadId);
    if(begin == end) {
        return;
    }



    // Loop over this range of reads.
    string readName;
    string readMetaData;
    vector<Base> read;
    vector<Base> runLengthRead;
    vector<uint8_t> readRepeatCount;
    const auto fileBegin = buffer.begin();
    for(uint64_t i=begin; i!=end; i++) {

        // Locate the 4 line ends corresponding to this read.
        const auto thisReadLineEnds = lineEnds.begin() + i * 4;

        // Locate the header line for this read.
        const auto headerBegin = fileBegin + ((i == 0) ? 0 : (1 + *(thisReadLineEnds - 1)));
        const auto headerEnd = fileBegin + thisReadLineEnds[0];
        DINARA_ASSERT(headerEnd > headerBegin);

        // Locate the sequence line for this read.
        const auto sequenceBegin = headerEnd + 1;
        const auto sequenceEnd = fileBegin + thisReadLineEnds[1];
        DINARA_ASSERT(sequenceEnd > sequenceBegin);

        // Locate the line containing the '+' for this read.
        const auto plusBegin = sequenceEnd + 1;
        const auto plusEnd = fileBegin + thisReadLineEnds[2];
        DINARA_ASSERT(plusEnd > plusBegin);

        // Locate the header line for this read.
        const auto scoresBegin = plusEnd + 1;
        const auto scoresEnd = fileBegin + thisReadLineEnds[3];
        DINARA_ASSERT(scoresEnd > scoresBegin);

        // Check the header line.
        if (headerEnd == headerBegin) {
            throw runtime_error("Empty header line for read at offset " +
                to_string(headerBegin - fileBegin) + ".");
        }
        if (*headerBegin != '@') {
            throw runtime_error("Read at offset " +
                to_string(headerBegin - fileBegin) +
                " does not begin with \"@\".");
        }

        // Extract the read name.
        // It starts immediately following the '@' and ends at
        // first white space.
        readName.clear();
        const auto nameBegin = headerBegin + 1;
        for (auto it = nameBegin; it != headerEnd; ++it) {
            const char c = *it;
            if (std::isspace(c)) {
                break;
            }
            readName.push_back(c);
        }
        if(readName.empty()) {
            throw runtime_error("Empty name for read at offset " +
                to_string(headerBegin - fileBegin) + ".");
        }

        // Extract the read meta data. It starts at the first non-space character
        // following the read name.
        readMetaData.clear();
        for (auto it = nameBegin + readName.size(); it != headerEnd; ++it) {
            const char c = *it;
            if (isspace(c) and readMetaData.empty()) {
                // Do nothing. Only start storing at the first non-space character.
            } else {
                readMetaData.push_back(c);
            }
        }


        // Check the line containing the plus.
        // It can contain arbitrary characters after the plus.
        // We already checked above that it is at least 1 character long.
        if (*plusBegin != '+') {
            throw runtime_error("Third line does not contain \"+\" for read " +
                                readName + " at offset " + to_string(headerBegin - fileBegin) + ".");
        }

        // Get the number of bases.
        const auto baseCount = sequenceEnd - sequenceBegin;
        if (scoresEnd - scoresBegin != baseCount) {
            throw runtime_error(
                    "Inconsistent numbers of bases and quality scores for read " +
                    readName + " at offset " +
                    to_string(headerBegin - fileBegin) + ": " +
                    to_string(baseCount) + " bases, " +
                    to_string(scoresEnd - scoresBegin) + " quality scores."
            );
        }

        // Get the bases.
        read.clear();
        for (auto it = sequenceBegin; it != sequenceEnd; ++it) {
            const char c = *it;
            const Base base = Base::fromCharacterNoException(c);
            if (!base.isValid()) {
                throw runtime_error("Invalid base " + string(1, c) + " for read " +
                                    readName + " at offset " + to_string(it - fileBegin) + ".");
            }
            read.push_back(base);
        }

        // If the read is too short, skip it.
        if (read.size() < minReadLength) {
            __sync_fetch_and_add(&discardedShortReadReadCount, 1);
            __sync_fetch_and_add(&discardedShortReadBaseCount, read.size());
            continue;
        }

        // Store the read.
        if(representation == 1) {
            if(computeRunLengthRepresentation(read, runLengthRead, readRepeatCount)) {
                thisThreadReadNames.appendVector(readName.begin(), readName.end());
                thisThreadReadMetaData.appendVector(readMetaData.begin(), readMetaData.end());
                thisThreadReads.append(runLengthRead);
                thisThreadReadRepeatCounts.appendVector(readRepeatCount);
            } else {
                __sync_fetch_and_add(&discardedBadRepeatCountReadCount, 1);
                __sync_fetch_and_add(&discardedBadRepeatCountBaseCount, read.size());
            }
        } else {
            thisThreadReadNames.appendVector(readName.begin(), readName.end());
            thisThreadReadMetaData.appendVector(readMetaData.begin(), readMetaData.end());
            thisThreadReads.append(read);
        }
    }
}



void ReadLoader::allocateBufferAndReadFile()
{
    if(isGzipped) {
        readGzipFile();
        return;
    }

    allocateBuffer();

    // Try reading using the requested setting of noCache/O_DIRECT.
    // If successful, all done.
    if(readFile(noCache)) {
        return;
    }

    // If there was failure and we are using noCache, try turning it off.
    // If successful, all done.
    if(noCache) {
        cout << "Turning off --Reads.noCache for " << fileName << endl;
        if(readFile(false)) {
            return;
        }
    }

    // If getting here, nothing worked.
    throw runtime_error("Error reading " + fileName);
}



// Read a gzip-compressed file, decompressing into the buffer.
// Uses the raw inflate API (like htslib's BGZF) instead of gzread
// to avoid double-buffering overhead.
void ReadLoader::readGzipFile()
{
    const auto t0 = std::chrono::steady_clock::now();

    const int64_t compressedSize = std::filesystem::file_size(fileName);

    // Read the entire compressed file into a temporary buffer.
    // This gives inflate contiguous input and avoids per-chunk read syscalls.
    vector<uint8_t> compressedData(compressedSize);
    {
        const int fd = ::open(fileName.c_str(), O_RDONLY);
        if(fd == -1) {
            throw runtime_error("Could not open " + fileName);
        }
        uint8_t* p = compressedData.data();
        int64_t remaining = compressedSize;
        while(remaining > 0) {
            const ssize_t n = ::read(fd, p, remaining);
            if(n <= 0) {
                ::close(fd);
                throw runtime_error("Error reading " + fileName);
            }
            p += n;
            remaining -= n;
        }
        ::close(fd);
    }

    // Estimate decompressed size and allocate output buffer.
    int64_t estimatedSize = compressedSize * 4;
    if(estimatedSize < 1024 * 1024) {
        estimatedSize = 1024 * 1024;
    }
    buffer.createNew(dataName("tmp-FastaBuffer"), pageSize);
    buffer.reserve(estimatedSize);
    buffer.resize(estimatedSize);

    // Initialize zlib for gzip decompression.
    // 15 + 16 = windowBits for gzip (auto-detect gzip/zlib header).
    z_stream zs = {};
    zs.next_in = compressedData.data();
    zs.avail_in = compressedSize;
    int ret = inflateInit2(&zs, 15 + 16);
    if(ret != Z_OK) {
        throw runtime_error("inflateInit2 failed for " + fileName);
    }

    int64_t totalOut = 0;
    // Decompress in 64 KB output chunks (matches htslib's BGZF_MAX_BLOCK_SIZE).
    const uint64_t outChunk = 64 * 1024;

    while(true) {
        // Grow output buffer if needed.
        if(totalOut + int64_t(outChunk) > int64_t(buffer.size())) {
            const int64_t newSize = max(int64_t(buffer.size()) * 2, totalOut + int64_t(outChunk));
            buffer.resize(newSize);
        }

        zs.next_out = reinterpret_cast<Bytef*>(&buffer[totalOut]);
        zs.avail_out = outChunk;

        ret = inflate(&zs, Z_SYNC_FLUSH);

        const int64_t produced = outChunk - zs.avail_out;
        totalOut += produced;

        if(ret == Z_STREAM_END) {
            break;
        }
        if(ret != Z_OK && ret != Z_BUF_ERROR) {
            inflateEnd(&zs);
            throw runtime_error("inflate failed for " + fileName +
                ": " + (zs.msg ? zs.msg : "unknown error"));
        }
        if(produced == 0 && zs.avail_in == 0) {
            break; // No progress, no input — done.
        }
    }

    inflateEnd(&zs);

    // Trim buffer to actual decompressed size.
    buffer.resize(totalOut);
    fileSize = totalOut;

    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = seconds(t1 - t0);

    performanceLog << "Decompressed gzip file: " << compressedSize << " -> "
        << totalOut << " bytes." << endl;
    performanceLog << "Decompress time: " << t01 << " s." << endl;
    performanceLog << "Decompress rate: " << double(totalOut) / t01 << " bytes/s." << endl;
}



void ReadLoader::allocateBuffer()
{
    const auto t0 = std::chrono::steady_clock::now();

    // Create a buffer to contain the entire file.
    fileSize = std::filesystem::file_size(fileName);
    buffer.createNew(dataName("tmp-FastaBuffer"), pageSize);

    // Do reserve before resize, to force using exactly the
    // amount of memory necessary and nothing more.
    buffer.reserve(fileSize);
    buffer.resize(fileSize);

    const auto t1 = std::chrono::steady_clock::now();

    performanceLog <<  "File size: " << fileSize << " bytes." << endl;
    performanceLog << "Buffer allocate time: " << seconds(t1 - t0) << " s." << endl;
}



// Read an entire file into the buffer.
bool ReadLoader::readFile(bool useODirect)
{
    const auto t0 = std::chrono::steady_clock::now();

    // Set up flags to open the file.
    int flags = O_RDONLY;
    if(useODirect) {
        flags |= O_DIRECT;
    }

    // Open the input file.
    const int fileDescriptor = ::open(fileName.c_str(), flags);
    if(fileDescriptor == -1) {
        return false;
    }

    // Read it in.
    char* bufferPointer = &buffer[0];
    uint64_t bufferCapacity = buffer.capacity();
    uint64_t bytesToRead = fileSize;
    while(bytesToRead) {
        const int64_t bytesRead = ::read(fileDescriptor, bufferPointer, bufferCapacity);
        if(bytesRead == -1) {
            ::close(fileDescriptor);
            return false;
        }
        bufferPointer += bytesRead;
        bytesToRead -= bytesRead;
        bufferCapacity -= bytesRead;
    }
    ::close(fileDescriptor);
    const auto t1 = std::chrono::steady_clock::now();
    const double t01 = seconds(t1 - t0);

    performanceLog << "Read time: " << t01 << " s." << endl;
    performanceLog << "Read rate: " << double(fileSize) / t01 << " bytes/s." << endl;
    return true;
}



// Find all of the line ends ('\n') in the buffer.
void ReadLoader::findLineEnds()
{
    // Each thread finds line ends in a block of the file
    // and stores them in its own vector.
    threadLineEnds.resize(threadCount);
    runThreads(&ReadLoader::findLineEndsThreadFunction, threadCount);

    // Combine the line ends founds by all threads.
    for(size_t threadId=0; threadId<threadCount; threadId++) {
        const vector<uint64_t>& thisThreadLineEnds = threadLineEnds[threadId];
        lineEnds.insert(lineEnds.end(),
            thisThreadLineEnds.begin(), thisThreadLineEnds.end());
    }
    threadLineEnds.clear();
}
void ReadLoader::findLineEndsThreadFunction(size_t threadId)
{
    // Access the vector where this thread will store the line ends it finds.
    vector<uint64_t>& thisThreadLineEnds = threadLineEnds[threadId];

    // Compute the file block assigned to this thread.
    uint64_t begin, end;
    tie(begin, end) = splitRange(0, buffer.size(), threadCount, threadId);
    if(begin == end) {
        return;
    }

    // Look for line ends in this block.
    for(uint64_t offset=begin; offset!=end; offset++) {
        if(buffer[offset] == '\n') {
            thisThreadLineEnds.push_back(offset);
        }
    }
}



// Create the name to be used for a MemoryMapped object.
string ReadLoader::dataName(
    const string& dataName) const
{
    if(dataNamePrefix.empty()) {
        return "";
    } else {
        return dataNamePrefix + dataName;
    }

}
string ReadLoader::threadDataName(
    size_t threadId,
    const string& dataName) const
{
    if(dataNamePrefix.empty()) {
        return "";
    } else {
        return dataNamePrefix + "tmp-ReadLoader-" + dataName + "-" + to_string(threadId);
    }

}



// Allocate space for the data structures where
// each thread stores the reads it found and their names.
void ReadLoader::allocatePerThreadDataStructures()
{
    threadReadNames.resize(threadCount);
    threadReadMetaData.resize(threadCount);
    threadReads.resize(threadCount);
    threadReadRepeatCounts.resize(threadCount); // Not actually used is representation==1
}


void ReadLoader::allocatePerThreadDataStructures(size_t threadId)
{
    threadReadNames[threadId] = make_unique< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    threadReadNames[threadId]->createNew(
        threadDataName(threadId, "ReadNames"), pageSize);
    threadReadMetaData[threadId] = make_unique< MemoryMapped::VectorOfVectors<char, uint64_t> >();
    threadReadMetaData[threadId]->createNew(
        threadDataName(threadId, "ReadMetaData"), pageSize);
    threadReads[threadId] = make_unique<LongBaseSequences>();
    threadReads[threadId]->createNew(
        threadDataName(threadId, "Reads"), pageSize);
    threadReadRepeatCounts[threadId] = make_unique< MemoryMapped::VectorOfVectors<uint8_t, uint64_t> >();
    if(representation == 1) {
        threadReadRepeatCounts[threadId]->createNew(
            threadDataName(threadId, "ReadRepeatCounts"), pageSize);
    }
}



// Store the reads computed by each thread and free
// the per-thread data structures.
void ReadLoader::storeReads()
{
    // Loop over all threads.
    for(size_t threadId=0; threadId<threadCount; threadId++) {

        // Access the names.
        MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadNames =
            *(threadReadNames[threadId]);

        // Access the meta data.
        MemoryMapped::VectorOfVectors<char, uint64_t>& thisThreadReadMetaData =
            *(threadReadMetaData[threadId]);

        // Access the reads.
        LongBaseSequences& thisThreadReads = *(threadReads[threadId]);
        const size_t n = thisThreadReadNames.size();
        DINARA_ASSERT(thisThreadReads.size() == n);
        DINARA_ASSERT(thisThreadReadMetaData.size() == n);

        // Access the repeat counts.
        MemoryMapped::VectorOfVectors<uint8_t, uint64_t>& thisThreadReadRepeatCounts =
            *threadReadRepeatCounts[threadId];
        if(representation == 1) {
            DINARA_ASSERT(thisThreadReadRepeatCounts.size() == n);
        } else {
            DINARA_ASSERT(not thisThreadReadRepeatCounts.isOpen());
        }

        // Store the reads.
        for(size_t i=0; i<n; i++) {
            reads.readNames.appendVector(thisThreadReadNames.begin(i), thisThreadReadNames.end(i));
            reads.readMetaData.appendVector(thisThreadReadMetaData.begin(i), thisThreadReadMetaData.end(i));
            reads.reads.append(thisThreadReads[i]);
            if(representation == 1) {
                const size_t j = reads.readRepeatCounts.size();
                reads.readRepeatCounts.appendVector(thisThreadReadRepeatCounts.size(i));
                copy(
                    thisThreadReadRepeatCounts.begin(i),
                    thisThreadReadRepeatCounts.end(i),
                    reads.readRepeatCounts.begin(j));
            }
        }

        // Remove the data structures used by this thread.
        thisThreadReadNames.remove();
        thisThreadReadMetaData.remove();
        thisThreadReads.remove();
        if(representation == 1) {
            thisThreadReadRepeatCounts.remove();
        }
    }

    // Clear the per-thread data structures.
    threadReadNames.clear();
    threadReadMetaData.clear();
    threadReads.clear();
    threadReadRepeatCounts.clear();

    // Free up unused allocated memory.
    reads.readNames.unreserve();
    reads.readMetaData.unreserve();
    if(representation == 1) {
        reads.readRepeatCounts.unreserve();
    }
    reads.reads.unreserve();

    // Allocate enough space for readFlags which are populated later.
    reads.readFlags.resize(reads.readCount());
}

