#!/usr/bin/python3

import os
import re
import shutil
import subprocess
import tempfile
import urllib.request

# Define the local build directory
HOME = os.path.expanduser("~")
DINARA_BUILD_DIR = os.path.join(HOME, ".dinaraBuild")
INCLUDE_DIR = os.path.join(DINARA_BUILD_DIR, "include")
LIB_DIR = os.path.join(DINARA_BUILD_DIR, "lib")

# Clean and recreate the build directory structure
def initializeBuildDirectory():
    print(f"Initializing build directory: {DINARA_BUILD_DIR}")
    if os.path.exists(DINARA_BUILD_DIR):
        print(f"Removing existing directory: {DINARA_BUILD_DIR}")
        shutil.rmtree(DINARA_BUILD_DIR)
    
    os.makedirs(INCLUDE_DIR, exist_ok=True)
    os.makedirs(LIB_DIR, exist_ok=True)
    print("Created include and lib directories.")

def isArm():
    return subprocess.getoutput("uname -p") == "aarch64"

def getRustToolchainHashFromLibrary(libPath):
    """Extract the Rust toolchain commit hash embedded in a static library.
    
    Rust static libraries contain paths like /rustc/<commit-hash>/library/...
    We extract this hash to determine which Rust version was used to build it.
    """
    if not os.path.exists(libPath):
        return None
    try:
        # Use strings to extract readable strings from the library, then grep for rustc path
        output = subprocess.getoutput("strings " + libPath + " | grep -oP '/rustc/[a-f0-9]{40}' | head -1")
        if output:
            # Extract the 40-character commit hash
            parts = output.strip().split('/')
            if len(parts) >= 3:
                return parts[2]
    except:
        pass
    return None

def checkRustLibrariesConsistency():
    """Check that all existing Rust libraries were built with the same toolchain.
    
    This detects the case where different Rust libraries were built at different
    times with different Rust versions, which causes linker errors.
    """
    rustLibraries = [
        (os.path.join(LIB_DIR, "libastarpa_c.a"), "libastarpa_c.a"),
        (os.path.join(LIB_DIR, "libpoasta_c.a"), "libpoasta_c.a"),
        (os.path.join(LIB_DIR, "libsimd_minimizers_c.a"), "libsimd_minimizers_c.a"),
    ]
    
    # Collect hashes for all existing libraries
    existingLibs = []
    for libPath, libName in rustLibraries:
        if os.path.exists(libPath):
            libHash = getRustToolchainHashFromLibrary(libPath)
            if libHash:
                existingLibs.append((libPath, libName, libHash))
    
    # If we have fewer than 2 libraries, no consistency check needed
    if len(existingLibs) < 2:
        return True
    
    # Check if all hashes match
    firstHash = existingLibs[0][2]
    mismatched = []
    for libPath, libName, libHash in existingLibs:
        if libHash != firstHash:
            mismatched.append((libName, libHash))
    
    if mismatched:
        # Build error message showing all libraries and their versions
        libDetails = "\n".join(
            "  - %s: %s" % (libName, libHash[:12]) 
            for _, libName, libHash in existingLibs
        )
        raise Exception(
            "\n" + "="*70 + "\n" +
            "ERROR: Rust libraries built with different toolchains!\n" +
            "="*70 + "\n" +
            "The following Rust libraries were built with different Rust versions:\n" +
            libDetails + "\n" +
            "\n" +
            "All Rust static libraries must be built with the SAME toolchain\n" +
            "to avoid linker errors (duplicate 'rust_eh_personality' symbols).\n" +
            "\n" +
            "To fix this, remove all existing Rust libraries and reinstall:\n" +
            "  rm -rf " + INCLUDE_DIR + "/astarpa " + INCLUDE_DIR + "/poasta " + INCLUDE_DIR + "/simd-minimizers\n" +
            "  rm -f " + LIB_DIR + "/libastarpa_c.a " + LIB_DIR + "/libpoasta_c.* " + LIB_DIR + "/libsimd_minimizers_c.*\n" +
            "  python3 scripts/InstallPrerequisites-Ubuntu.py\n" +
            "="*70
        )
    
    print("All existing Rust libraries were built with the same toolchain (%s)." % firstHash[:12])
    return True

def runCommand(command):
    if(os.system(command)):
        raise Exception("Error running command: " + command)


def patchRustPortableSimdLaneCountRemoval(cratePath):
    """Patch Rust sources that still use the old portable_simd LaneCount API.

    Newer nightlies expose `std::simd::Simd<T, LANES>` directly and removed
    `std::simd::{LaneCount, SupportedLaneCount}`. Some dependencies in the
    astar-pairwise-aligner workspace still import these names, which breaks
    compilation on recent toolchains.
    """
    if not os.path.isdir(cratePath):
        return 0

    patchedFiles = 0
    for root, _, files in os.walk(cratePath):
        for name in files:
            if not name.endswith(".rs"):
                continue
            path = os.path.join(root, name)
            try:
                with open(path, "r", encoding="utf-8") as f:
                    text = f.read()
            except Exception:
                continue

            original = text

            # Remove direct imports like: use std::simd::{LaneCount, SupportedLaneCount};
            text = re.sub(
                r"(?m)^\s*use\s+std::simd::\{\s*LaneCount\s*,\s*SupportedLaneCount\s*\};\s*\n",
                "",
                text,
            )
            text = re.sub(
                r"(?m)^\s*use\s+std::simd::\{\s*SupportedLaneCount\s*,\s*LaneCount\s*\};\s*\n",
                "",
                text,
            )

            # Remove nested imports like: simd::{LaneCount, SupportedLaneCount},
            text = re.sub(
                r"\bsimd::\{\s*LaneCount\s*,\s*SupportedLaneCount\s*\}\s*,\s*",
                "",
                text,
            )
            text = re.sub(
                r"\bsimd::\{\s*SupportedLaneCount\s*,\s*LaneCount\s*\}\s*,\s*",
                "",
                text,
            )

            # Remove bounds like: LaneCount<LANES>: SupportedLaneCount,
            text = re.sub(
                r"(?m)^\s*LaneCount\s*<\s*[^>]+\s*>\s*:\s*SupportedLaneCount\s*,?\s*\n",
                "",
                text,
            )
            text = re.sub(
                r"\bLaneCount\s*<\s*[^>]+\s*>\s*:\s*SupportedLaneCount\s*,?\s*",
                "",
                text,
            )

            # Clean up empty `where` clauses left behind (common formatting patterns).
            text = re.sub(r"(?m)^\s*where\s*\n\s*\{", "{", text)
            text = re.sub(r"\bwhere\s*\{", "{", text)

            # Remove now-empty simd import blocks if they occur.
            text = re.sub(r"(?m)^\s*use\s+std::simd::\{\s*\};\s*\n", "", text)

            if text != original:
                with open(path, "w", encoding="utf-8") as f:
                    f.write(text)
                patchedFiles += 1

    return patchedFiles


def patchRustObsoleteFeatureGates(workspaceRoot):
    """Patch obsolete Rust feature gates that fail on recent nightly toolchains.

    Upstream astar-pairwise-aligner currently enables `bigint_helper_methods`
    in pa-bitpacking, but that gate is no longer recognized by current nightlies.
    Removing just that gate keeps the crate compatible with modern nightly.
    """
    if not os.path.isdir(workspaceRoot):
        return 0

    patchedFiles = 0
    for root, _, files in os.walk(workspaceRoot):
        for name in files:
            if not name.endswith(".rs"):
                continue
            path = os.path.join(root, name)
            try:
                with open(path, "r", encoding="utf-8") as f:
                    text = f.read()
            except Exception:
                continue

            original = text

            # Common multiline feature-list style:
            # #![feature(
            #   bigint_helper_methods,
            #   ...
            # )]
            text = re.sub(
                r"(?m)^\s*bigint_helper_methods\s*,\s*\n",
                "",
                text,
            )

            # Single-line style fallbacks.
            text = re.sub(r",\s*bigint_helper_methods\b", "", text)
            text = re.sub(r"\bbigint_helper_methods\s*,\s*", "", text)

            if text != original:
                with open(path, "w", encoding="utf-8") as f:
                    f.write(text)
                patchedFiles += 1

    return patchedFiles
        
def installPackage(package):
    runCommand("sudo apt-get install --assume-yes " + package)

def installAptPackages():
    packages = [
    "git",
    "g++",
    "make",
    "cmake",
    "libboost-system-dev", 
    "libboost-program-options-dev",
    "libboost-graph-dev",
    "libboost-chrono-dev",
    "libpng-dev", 
    "libblas-dev", 
    "liblapack-dev",
    "gfortran",
    "ncbi-blast+",
    "graphviz",
    "gnuplot",
    "python3-dev", 
    "libsimde-dev",
    "zlib1g-dev",
    "libcli11-dev",
    "libhtscodecs-dev",
    "libbz2-dev",
    "liblzma-dev",
    "libdeflate-dev",
    "libssl-dev",
    "libnghttp2-dev",
    "libpsl-dev",
    "libssh-dev",
    "libbrotli-dev",
    "libzstd-dev",
    "librtmp-dev",
    "libldap-dev",
    "libidn2-dev",
    "libunistring-dev",
    ]
    runCommand("sudo apt-get install --assume-yes " + " ".join(packages))

    # libcurl4-openssl-dev and libcurl4-gnutls-dev conflict with each other.
    # Install whichever one is compatible with the current system.
    if os.system("sudo apt-get install --assume-yes libcurl4-openssl-dev") != 0:
        print("libcurl4-openssl-dev conflicts; trying libcurl4-gnutls-dev...")
        if os.system("sudo apt-get install --assume-yes libcurl4-gnutls-dev") != 0:
            print("WARNING: Could not install any libcurl dev package. "
                  "If one is already installed, this is fine.")



# We don't use Ubuntu package libseqan2-dev because
# it does not include the fix for this issue:
# https://github.com/seqan/seqan/issues/2524
# Instead, we clone the GitHub seqan/seqan repository, 
# then copy seqan/include/seqan to /usr/include/seqan.
def installSeqan():

    # The path where the include files will go.
    installPath = os.path.join(INCLUDE_DIR, "seqan")
    
    # First check that this path does not exist.
    if os.path.exists(installPath):
        print("The seqan install path %s already exists. Skipping installation." % installPath)
        return
        
    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building seqan library using temporary directory", temporaryDirectory)
        
        # Change to the temporary directory.
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone the Github repository.
        runCommand("git clone https://github.com/seqan/seqan.git")
        
        # Copy the include files.
        shutil.copytree("seqan/include/seqan", installPath)

        # Change back to the original directory.
        os.chdir(oldDirectory)



def installPybind11():
    try:
        # This works for Ubuntu 22.04 and older
        print("Attempting pybind11 installation using pip3")
        installPackage("python3-pip")
        runCommand("sudo pip3 install pybind11")
    except:
        # This works for Ubuntu 24.04 (and newer, presumably).
        print("Pybind11 installation using pip3 did not work, installing using apt-get.")
        installPackage("python3-pybind11")


        
def installSpoa():
    # The spoa library is available in the stable Ubuntu repository, but
    # without the static version.
    # So we have to build it from source.
    
    if os.path.exists(os.path.join(INCLUDE_DIR, "spoa/spoa.hpp")):
        print("spoa header found. Skipping installation.")
        return

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building spoa library using temporary directory", temporaryDirectory)
        
        # Change to the temporary directory.
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
                
        # Get the code.
        runCommand("sudo apt-get install curl")
        runCommand("curl -L https://github.com/rvaser/spoa/archive/refs/tags/4.0.8.tar.gz -o 4.0.8.tar.gz")
        runCommand("tar -xvf 4.0.8.tar.gz")
    
        # Set spoa build flags.
        if isArm():
            spoaBuildFlags="-DCMAKE_BUILD_TYPE=Release -Dspoa_build_tests=OFF"
        else:
            # The spoa dispatcher feature selects code at run time based on available hardware features,
            # which can improve performance.
            # However, in spoa v4.0.8 it introduces two additional dependencies:
            # - USCiLab/cereal
            # - google/cpu_features
            # To avoid these additional dependencies, we turn off the dispatcher feature for now.
            # We could turn it back on if we see significant performance degradation in this area.
            spoaBuildFlags = "-DCMAKE_BUILD_TYPE=Release -Dspoa_optimize_for_portability=ON -Dspoa_build_tests=OFF"

        # Add install prefix
        spoaBuildFlags += " -DCMAKE_INSTALL_PREFIX=" + DINARA_BUILD_DIR

        # Build the shared library.
        print("\n******** Building spoa shared library")
        os.mkdir("build")
        os.chdir("build")
        runCommand("cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=ON " + spoaBuildFlags)
        runCommand("make -j all")
        runCommand("make install")
        
        # Build the static library.
        os.chdir("..")
        print("\n******** Building spoa static library")
        os.mkdir("build-static")
        os.chdir("build-static")
        runCommand("cmake ../spoa-4.0.8 -DBUILD_SHARED_LIBS=OFF " + spoaBuildFlags)
        runCommand("make -j all")
        runCommand("make install")
        
        # Change back to the original directory.
        os.chdir(oldDirectory)


def setupRustToolchain():
    """Install and configure Rust with the latest nightly version.
    
    Astarpa requires nightly Rust for unstable features (#![feature(...)]).
    All Rust libraries are built with the same nightly toolchain to ensure
    they share the same std symbols, avoiding linker conflicts from
    duplicate rust_eh_personality definitions.
    """
    print("Setting up Rust toolchain...")
    installPackage("curl")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    rustupPath = os.path.expanduser("~/.cargo/bin/rustup")
    cbindgenPath = os.path.expanduser("~/.cargo/bin/cbindgen")
    
    if not os.path.exists(cargoPath):
        print("Rust not found. Installing Rust...")
        runCommand("curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y")
    
    # Prepend local cargo bin to PATH so that subprocesses (like cbindgen) find the correct cargo
    os.environ["PATH"] = os.path.dirname(cargoPath) + os.pathsep + os.environ["PATH"]
    
    # Update rustup and use latest nightly for consistency across all Rust libraries
    # Nightly is required because astarpa uses unstable Rust features
    print("Updating Rust to latest nightly...")
    runCommand(rustupPath + " update nightly")
    runCommand(rustupPath + " default nightly")
    
    if not os.path.exists(cbindgenPath):
        print("Installing cbindgen...")
        runCommand(cargoPath + " install cbindgen")
    else:
        print("cbindgen already installed.")
    
    return cargoPath, rustupPath, cbindgenPath


def installAstarpa():
    print("Installing astarpa-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    cbindgenPath = os.path.expanduser("~/.cargo/bin/cbindgen")

    installPath = os.path.join(INCLUDE_DIR, "astarpa")
    libPath = os.path.join(LIB_DIR, "libastarpa_c.a")
    
    if os.path.exists(installPath + "/astarpa.h") and os.path.exists(libPath):
        print("astarpa header found in %s/astarpa.h. Skipping installation." % installPath)
        return

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building astarpa library using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner.git")
        os.chdir("astar-pairwise-aligner/astarpa-c")
        
        # Patch for recent nightly Rust (requires unsafe(...) wrapper for no_mangle)
        print("Patching astarpa-c for recent nightly Rust compatibility...")
        runCommand("sed -i 's/#\\[no_mangle\\]/#\\[unsafe(no_mangle)\\]/g' src/lib.rs")
        runCommand("sed -i '1i #![allow(unsafe_op_in_unsafe_fn)]' src/lib.rs")

        # Patch workspace crates that still use the old portable_simd LaneCount API.
        # Recent nightlies removed `std::simd::{LaneCount, SupportedLaneCount}`.
        workspaceRoot = os.path.abspath(os.path.join(os.getcwd(), ".."))
        patched = patchRustPortableSimdLaneCountRemoval(os.path.join(workspaceRoot, "pa-bitpacking"))
        if patched:
            print(f"Patched portable_simd LaneCount API in {patched} Rust files.")

        patched = patchRustObsoleteFeatureGates(workspaceRoot)
        if patched:
            print(f"Patched obsolete Rust feature gates in {patched} Rust files.")

        # Build
        runCommand(cargoPath + " build --release")
        
        # Generate header
        runCommand(cbindgenPath + " --lang c --cpp-compat --crate astarpa-c -o astarpa.h")
        
        # Install
        # Install
        print("Installing astarpa to " + DINARA_BUILD_DIR + "...")
        if not os.path.exists(installPath):
            os.makedirs(installPath, exist_ok=True)
        runCommand("cp astarpa.h " + installPath)
        runCommand("cp ../target/release/libastarpa_c.a " + LIB_DIR)
        
        os.chdir(oldDirectory)


def installPoasta():
    print("Installing poasta-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    
    installPath = os.path.join(INCLUDE_DIR, "poasta")
    libPath = os.path.join(LIB_DIR, "libpoasta_c.a")
    
    if os.path.exists(installPath + "/poasta.h") and os.path.exists(libPath):
        print("poasta-c header found in %s/poasta.h. Skipping installation." % installPath)
        return
        
    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building poasta-c library using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/kokyriakidis/poasta-c.git")
        os.chdir("poasta-c")
        
        # Build
        runCommand("RUSTFLAGS='-C target-cpu=native' " + cargoPath + " build --release")
        
        # Install
        # Install
        print("Installing poasta-c to " + DINARA_BUILD_DIR + "...")
        if not os.path.exists(installPath):
            os.makedirs(installPath, exist_ok=True)
        runCommand("cp poasta.h " + installPath)
        runCommand("cp target/release/libpoasta_c.so " + LIB_DIR)
        runCommand("cp target/release/libpoasta_c.a " + LIB_DIR)
        
        os.chdir(oldDirectory)

def installSimdMinimizers():
    print("Installing simd-minimizers-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    
    installPath = os.path.join(INCLUDE_DIR, "simd-minimizers")
    libPath = os.path.join(LIB_DIR, "libsimd_minimizers_c.a")
    
    # Force reinstall to get latest version with Syncmers support
    if os.path.exists(installPath):
        shutil.rmtree(installPath)
    if os.path.exists(libPath):
        os.remove(libPath)
    if os.path.exists(os.path.join(LIB_DIR, "libsimd_minimizers_c.so")):
        os.remove(os.path.join(LIB_DIR, "libsimd_minimizers_c.so"))
        
    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building simd-minimizers-c library using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/kokyriakidis/simd-minimizers-c.git")
        os.chdir("simd-minimizers-c")
        
        # Build with target-cpu=native for best performance
        runCommand("RUSTFLAGS='-C target-cpu=native' " + cargoPath + " build --release")
        
        # Install
        print("Installing simd-minimizers-c to " + DINARA_BUILD_DIR + "...")
        if not os.path.exists(installPath):
            os.makedirs(installPath, exist_ok=True)
            
        # Robust header finding
        if os.path.exists("simd_minimizers.h"):
             runCommand("cp simd_minimizers.h " + installPath)
        elif os.path.exists("include/simd_minimizers.h"):
             runCommand("cp include/simd_minimizers.h " + installPath)
        else:
             print("Warning: simd_minimizers.h not found in root/include, searching...")
             runCommand(f"find . -name 'simd_minimizers.h' -exec cp {{}} {installPath} \\;")
             
        # Verify header was installed
        if not os.path.exists(os.path.join(installPath, "simd_minimizers.h")):
             raise Exception("Failed to install simd_minimizers.h")

        runCommand("cp target/release/libsimd_minimizers_c.so " + LIB_DIR)
        runCommand("cp target/release/libsimd_minimizers_c.a " + LIB_DIR)
        
        os.chdir(oldDirectory)


def installBasePrerequisites():
    """Install apt packages, seqan, pybind11, spoa, and set up Rust toolchain."""
    initializeBuildDirectory()
    checkRustLibrariesConsistency()
    installAptPackages()
    installSeqan()
    installPybind11()
    installSpoa()
    setupRustToolchain()

def installHtslib():
    """Build and install htslib from source without curl support.

    Dinara does not use htslib's network features, and statically linking
    libcurl pulls in a large chain of transitive dependencies (ssh, gnutls,
    gssapi, brotli, etc.) that are difficult to satisfy with static archives.
    Building without --enable-libcurl avoids this entirely.
    """
    htslibVersion = "1.21"
    htslibTarball = f"htslib-{htslibVersion}.tar.bz2"
    htslibUrl = f"https://github.com/samtools/htslib/releases/download/{htslibVersion}/{htslibTarball}"
    htslibInstallMarker = "/usr/local/lib/libhts.a"

    if os.path.exists(htslibInstallMarker):
        # Check if existing htslib was built with curl (which we want to avoid).
        import subprocess
        nmOut = subprocess.run(["nm", htslibInstallMarker],
                               capture_output=True, text=True).stdout
        if "curl_" not in nmOut:
            print("htslib already installed without curl at " + htslibInstallMarker + ". Skipping.")
            return
        print("Existing htslib has curl support. Rebuilding without curl...")

    else:
        print(f"Installing htslib {htslibVersion} (without libcurl)...")

    with tempfile.TemporaryDirectory() as tmpDir:
        tarballPath = os.path.join(tmpDir, htslibTarball)
        urllib.request.urlretrieve(htslibUrl, tarballPath)
        runCommand(f"tar xjf {tarballPath} -C {tmpDir}")

        sourceDir = os.path.join(tmpDir, f"htslib-{htslibVersion}")
        runCommand(f"cd {sourceDir} && ./configure --prefix=/usr/local "
                   f"--disable-libcurl --disable-gcs --disable-s3 "
                   f"--enable-bz2 --enable-lzma --with-libdeflate")
        runCommand(f"make -C {sourceDir} -j")
        runCommand(f"sudo make -C {sourceDir} install")

    print("htslib installed to /usr/local (without libcurl).")


def installAbpoa():
    # abPOA is built directly here (previously it was a byproduct of the
    # now-removed shasta2 install). The distro package and abPOA releases only
    # ship the executable, and neither the upstream CMakeLists (unconditional
    # -march=native, non-portable) nor the upstream Makefile (no shared lib)
    # fit our needs, so we compile the sources directly. Recipe mirrors
    # shasta2's BuildAbpoa.py: portable flags, both static (linked by the
    # executable, needs -fPIC) and shared libraries.
    print("Installing abPOA (static + shared) into " + DINARA_BUILD_DIR + "...")

    abpoaIncludeDir = os.path.join(INCLUDE_DIR, "abpoa")

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building abPOA using temporary directory", temporaryDirectory)
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)

        runCommand("git clone https://github.com/yangao07/abPOA.git")
        # Pin to an exact release for reproducible builds.
        runCommand("git -C abPOA checkout v1.5.6")
        # abpoa.h -> simd_instruction.h -> "simde/simde/x86/*.h", where simde is
        # a git submodule. Without this init the simde tree is empty and every
        # TU that includes <abpoa.h> fails to compile.
        runCommand("git -C abPOA submodule update --init --recursive")
        os.chdir("abPOA")

        # Common compile flags (portable; no -march=native).
        commonFlags = ("-I ../include -O3 -Wall -Wno-unused-function "
                       "-Wno-misleading-indentation -Wno-stringop-overflow "
                       "-fno-tree-vectorize")

        os.mkdir("lib")
        os.chdir("src")

        # Static library (with -fPIC so it can be linked into shared objects).
        print("Building the abPOA static library.")
        runCommand("cc -c -fPIC " + commonFlags + " *.c")
        runCommand("ar -csr ../lib/libabpoa.a *.o")
        runCommand("rm *.o")

        # Shared library.
        print("Building the abPOA shared library.")
        runCommand("cc -c -fPIC " + commonFlags + " *.c")
        runCommand("cc -shared -o ../lib/libabpoa.so *.o")
        runCommand("rm *.o")

        os.chdir("..")

        # Install libraries.
        runCommand("cp lib/libabpoa.a " + os.path.join(LIB_DIR, "libabpoa.a"))
        runCommand("cp lib/libabpoa.so " + os.path.join(LIB_DIR, "libabpoa.so"))

        # Install headers. Non-excluded TUs include both <abpoa.h> and
        # "abpoa/abpoa.h", so populate both include/abpoa/ and include/.
        # Copy the FULL header set (not just abpoa.h) to BOTH locations:
        # abpoa.h pulls in "simd_instruction.h" via a quoted include, which in
        # turn includes the "simde/" tree. Copying only abpoa.h leaves those
        # siblings unresolved and breaks any TU that does #include <abpoa.h>.
        if not os.path.exists(abpoaIncludeDir):
            os.makedirs(abpoaIncludeDir, exist_ok=True)
        runCommand("cp include/*.h " + abpoaIncludeDir)
        runCommand("cp -r include/simde " + abpoaIncludeDir)
        runCommand("cp include/*.h " + INCLUDE_DIR)
        runCommand("cp -r include/simde " + INCLUDE_DIR)

        os.chdir(oldDirectory)

    print("abPOA installed.")



def installTheseusLib():
    # kokyriakidis fork of theseus-lib, pericles branch. The fork is a SUPERSET
    # of upstream: it keeps the API shasta2's theseusWrapper.cpp expects
    # (default Heuristics() and the 6-arg align) AND adds the fork-only API that
    # dinara's MSA code uses (align_from + the multi-segment TheseusMSA ctor).
    # Using the fork therefore satisfies both shasta2 and dinara from a single
    # install, so dinara's multi-segment MSA sources can stay enabled.
    #
    # The fork emits GFA segment names as numeric node ids directly from
    # TheseusAlignerImpl::print_as_gfa_internal, so the legacy graph.h GFA
    # name patch upstream needed is not required here.
    print("Installing theseus-lib (kokyriakidis fork, pericles branch)...")

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building theseus-lib using temporary directory", temporaryDirectory)

        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)

        # Clone the fork's pericles branch, then pin to an exact commit
        # for reproducible builds. The repository is public, so no
        # authentication is required.
        runCommand(
            "git clone -b pericles "
            "https://github.com/kokyriakidis/theseus-lib-multi-segment.git "
            "theseus-lib"
        )
        os.chdir("theseus-lib")
        runCommand("git checkout 13a4321937b9bf36314c114c066efb13fa1d4dbc")

        # Build and install static library into DINARA_BUILD_DIR
        os.mkdir("build")
        os.chdir("build")
        runCommand(
            "cmake .. "
            "-DCMAKE_BUILD_TYPE=Release "
            "-DBUILD_SHARED_LIBS=OFF "
            f"-DCMAKE_INSTALL_PREFIX={DINARA_BUILD_DIR}"
        )
        runCommand("cmake --build . -j")
        runCommand("cmake --install .")

        os.chdir(oldDirectory)

    print("theseus-lib installed.")


# All installable targets in dependency order.
# Each entry: (name, function, description)
INSTALL_TARGETS = [
    ("base",             installBasePrerequisites, "apt packages, seqan, pybind11, spoa, Rust toolchain"),
    ("htslib",           installHtslib,            "htslib (without libcurl, for static linking)"),
    ("astarpa",          installAstarpa,           "astarpa alignment library (Rust)"),
    ("poasta",           installPoasta,            "poasta alignment library (Rust)"),
    ("simd-minimizers",  installSimdMinimizers,    "SIMD minimizers library (Rust)"),
    ("theseus",          installTheseusLib,        "theseus-lib POA aligner (C++)"),
    ("abpoa",            installAbpoa,             "abPOA POA aligner (static + shared)"),
]

def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Install Dinara build prerequisites.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="Examples:\n"
               "  Install everything:          ./InstallPrerequisites-Ubuntu.py\n"
               "  Install only abpoa:          ./InstallPrerequisites-Ubuntu.py --only abpoa\n"
               "  Install multiple targets:    ./InstallPrerequisites-Ubuntu.py --only abpoa theseus\n"
               "  List available targets:      ./InstallPrerequisites-Ubuntu.py --list\n")
    parser.add_argument("--only", nargs="+", metavar="TARGET",
                        help="Install only the specified target(s). Skips base prerequisites unless 'base' is listed.")
    parser.add_argument("--list", action="store_true",
                        help="List available install targets and exit.")
    args = parser.parse_args()

    if args.list:
        print("Available install targets:")
        for name, _, desc in INSTALL_TARGETS:
            print(f"  {name:20s} {desc}")
        return

    target_map = {name: func for name, func, _ in INSTALL_TARGETS}

    if args.only:
        for name in args.only:
            if name not in target_map:
                print(f"Error: unknown target '{name}'. Use --list to see available targets.")
                return
            print(f"\n{'='*60}")
            print(f"Installing: {name}")
            print(f"{'='*60}")
            target_map[name]()
    else:
        # Full install: run everything in order.
        for name, func, _ in INSTALL_TARGETS:
            func()

    print("\nInstallation completed successfully.")

if __name__ == "__main__":
    main()
