cmake_minimum_required(VERSION 3.16)
project(dinara_tests)

# Repository root (works both when configuring `tests/` directly and when added
# via `add_subdirectory(tests)` from the repo root).
set(DINARA_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")

# C++ dialect. Scoped to CXX only (add_definitions would also apply it to the
# project's one C source, ksw2_extd2_sse.c, which clang rejects -std=c++23 for).
add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-std=c++23>)

# Debug symbols for tests
add_definitions(-ggdb3 -O0)

# Enable test-only hooks in dinara sources.
add_definitions(-DDINARA_TESTING)

add_definitions(-DDINARA_ENABLE_VARIANT_CLUSTERING=0)

# Include directories
include_directories(${DINARA_ROOT}/src)
include_directories(${DINARA_ROOT}/srcMain)
include_directories(${DINARA_ROOT}/external/catch2)
include_directories(${DINARA_ROOT}/src/ksw2)
include_directories(${DINARA_ROOT}/external/hifiasmCandidates)

# Figure out if we are on Intel (as opposed to ARM).
set(X86_64 OFF)
if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    set(X86_64 ON)
endif()

# 16-byte compare and swap.
if(X86_64)
    add_definitions(-mcx16)
endif()

# Dinara build directory for libraries
set(DINARA_BUILD_DIR $ENV{HOME}/.dinaraBuild)
include_directories(SYSTEM ${DINARA_BUILD_DIR}/include)
link_directories(${DINARA_BUILD_DIR}/lib)

# Catch2 main (provides main())
add_library(catch_main OBJECT catch_main.cpp)

# ============================================================================
# INTEGRATION TESTS (link against dinara sources)
# ============================================================================

# Get all dinara source files EXCEPT main.cpp
file(GLOB DINARA_SRC_FILES CONFIGURE_DEPENDS
    ${DINARA_ROOT}/src/*.cpp
    ${DINARA_ROOT}/srcMain/*.cpp
    ${DINARA_ROOT}/src/ksw2/ksw2_extd2_sse.c
)
# Remove main.cpp (it has its own main())
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*main\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerVariantClustering\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/mode3-AnchorFromHetSites\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/ClusterGraph\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerClusterGraph\\.cpp$")
# Theseus disabled: these files use fork-only theseus APIs (align_from,
# multi-segment constructor) that upstream theseus-lib does not provide.
# Excluded from the build until theseus usage is reworked. Callers are stubbed.
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerPhasingMSA\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerTheseusReadWindowMSA\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerDirectOverlapMSA\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerTheseusMarkerGraphMSA\\.cpp$")
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerMultiSegmentMSA\\.cpp$")
# AssemblerHifiasmEC.cpp excluded: it was the only consumer of the removed
# AlignedEvidenceStore. Overlap phasing now runs off OverlapCigarStore; the
# hifiasm EC parity path is unused and kept in-tree for reference only.
list(FILTER DINARA_SRC_FILES EXCLUDE REGEX ".*/AssemblerHifiasmEC\\.cpp$")

# Find required libraries
find_library(ASTARPA_LIB astarpa_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(POASTA_LIB poasta_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(SIMD_MINIMIZERS_LIB simd_minimizers_c PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(ABPOA_LIB abpoa PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
find_library(THESEUS_LIB theseus PATHS ${DINARA_BUILD_DIR}/lib NO_DEFAULT_PATH)
if(NOT THESEUS_LIB)
    set(THESEUS_LIB "")
endif()
find_package(ZLIB REQUIRED)

# Integration test executable
add_executable(test_integration
    test_integration.cpp
    test_overlap_cigar_store.cpp
    test_phasing_coordinates.cpp
    test_msa_het_kernel.cpp
    ${DINARA_SRC_FILES}
    $<TARGET_OBJECTS:catch_main>
)
set_source_files_properties(
    ${DINARA_ROOT}/src/MurmurHash2.cpp
    PROPERTIES COMPILE_OPTIONS "-Wno-implicit-fallthrough")

# Some toolchains in this environment can produce non-executable output files
# (missing +x), which makes CTest fail with "permission denied".
add_custom_command(TARGET test_integration POST_BUILD
    COMMAND chmod +x $<TARGET_FILE:test_integration>
    VERBATIM
)

# Link all required libraries
if(APPLE)
    target_link_libraries(test_integration
        hifiasm_overlaps
        atomic boost_system boost_program_options boost_chrono boost_serialization
        spoa png lapack blas gfortran m pthread
        ${ABPOA_LIB} ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB}
        ZLIB::ZLIB
    )
else()
    # Most of our dependencies are static .a libraries, so we need either strict
    # ordering or a link group to resolve transitive symbols.
    target_link_libraries(test_integration
        hifiasm_overlaps
        atomic boost_system boost_program_options boost_chrono boost_serialization
        spoa png lapack blas gfortran m pthread
        -Wl,--start-group
        ${ABPOA_LIB} ${ASTARPA_LIB} ${POASTA_LIB} ${SIMD_MINIMIZERS_LIB} ${THESEUS_LIB}
        ${HTS_LINK_FLAGS} ${CURL_LINK_FLAGS}
        z
        -Wl,--end-group
    )
endif()

# Match the main executable's static link flags.
if(NOT APPLE)
    set_target_properties(test_integration PROPERTIES
        LINK_FLAGS "-static -Wl,--allow-multiple-definition")
endif()

# Ensure the test binary uses the same curated Dinara dependency builds at runtime.
set_target_properties(test_integration PROPERTIES
    BUILD_RPATH "${DINARA_BUILD_DIR}/lib"
)

# ============================================================================
# LONG-MARKERS UNIT TEST (compiled with -DDINARA_LONG_MARKERS)
# ============================================================================
# Validates the capacity-128 Kmer / 256-bit KmerId primitives (Uint256,
# ShortBaseSequence<__uint128_t>, id() round-trip, reverseComplement,
# extractKmer) in isolation. This is a focused unit target rather than a full
# rebuild of dinara under the flag: only the minimal sources the primitives
# depend on are compiled here.
add_executable(test_long_markers
    test_long_markers.cpp
    ${DINARA_ROOT}/src/extractKmer.cpp
    ${DINARA_ROOT}/src/DINARA_ASSERT.cpp
    $<TARGET_OBJECTS:catch_main>
)
target_compile_definitions(test_long_markers PRIVATE DINARA_LONG_MARKERS)
target_link_libraries(test_long_markers pthread)
add_custom_command(TARGET test_long_markers POST_BUILD
    COMMAND chmod +x $<TARGET_FILE:test_long_markers>
    VERBATIM
)

# ============================================================================
# PAF IMPORT UNIT TEST
# ============================================================================
# Validates the pure PAF import helpers in src/PafImport.hpp (line parsing,
# line-aligned chunk splitting, canonical keying and longest-overlap dedup) in
# isolation. Focused unit target: PafImport.hpp only depends on ReadId.hpp and
# span.hpp, so only DINARA_ASSERT.cpp is needed for the assert symbol.
add_executable(test_paf_import
    test_paf_import.cpp
    ${DINARA_ROOT}/src/DINARA_ASSERT.cpp
    $<TARGET_OBJECTS:catch_main>
)
target_link_libraries(test_paf_import pthread)
add_custom_command(TARGET test_paf_import POST_BUILD
    COMMAND chmod +x $<TARGET_FILE:test_paf_import>
    VERBATIM
)

# ============================================================================
# HIFIASM CIGAR IMPORT UNIT TEST
# ============================================================================
# Validates normalizeHifiasmCigar in src/HifiasmCigarImport.hpp: reframing the
# hifiasm alignment-frame CIGAR into dinara's canonical read0/read1 frame for
# all four query/target-id-order x strand cases. Pure logic over the
# header-only OverlapCigarStore, so only DINARA_ASSERT.cpp is needed.
add_executable(test_hifiasm_cigar_import
    test_hifiasm_cigar_import.cpp
    ${DINARA_ROOT}/src/DINARA_ASSERT.cpp
    $<TARGET_OBJECTS:catch_main>
)
target_link_libraries(test_hifiasm_cigar_import pthread)
add_custom_command(TARGET test_hifiasm_cigar_import POST_BUILD
    COMMAND chmod +x $<TARGET_FILE:test_hifiasm_cigar_import>
    VERBATIM
)

# Enable CTest
enable_testing()
add_test(NAME IntegrationTests COMMAND test_integration)
add_test(NAME LongMarkersTests COMMAND test_long_markers)
add_test(NAME PafImportTests COMMAND test_paf_import)
add_test(NAME HifiasmCigarImportTests COMMAND test_hifiasm_cigar_import)

# CI/dev environments can have LD_LIBRARY_PATH pointing at /usr/local/lib first.
# Force our curated dependency builds (in ${DINARA_BUILD_DIR}/lib) to take precedence.
set_tests_properties(IntegrationTests PROPERTIES
    ENVIRONMENT "LD_LIBRARY_PATH=${DINARA_BUILD_DIR}/lib"
)
