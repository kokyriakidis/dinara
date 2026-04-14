#!/usr/bin/python3

import os
import re
import shutil
import subprocess
import tempfile
import urllib.request
import json

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
    ]
    runCommand("sudo apt-get install --assume-yes " + " ".join(packages))



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


def installBubbleFinder():
    print("Installing BubbleFinder...")

    installBinDir = os.path.join(HOME, ".local", "bin")
    installBinary = os.path.join(installBinDir, "BubbleFinder")
    sourceDir = os.path.join(HOME, "Downloads", "BubbleFinder")

    os.makedirs(installBinDir, exist_ok=True)

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building BubbleFinder using temporary directory", temporaryDirectory)

        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)

        if os.path.isdir(sourceDir):
            print("Updating existing BubbleFinder source at", sourceDir)
            runCommand("git -C " + sourceDir + " pull --ff-only")
        else:
            if os.path.exists(sourceDir):
                # A non-directory file exists at the source path (e.g. an old binary).
                # Remove it so git clone can create the directory.
                os.remove(sourceDir)
            print("Cloning BubbleFinder source to", sourceDir)
            runCommand("git clone https://github.com/algbio/BubbleFinder.git " + sourceDir)

        # BubbleFinder requires several git submodules (gbz, sdsl-lite, etc.).
        runCommand("git -C " + sourceDir + " submodule update --init --recursive")

        buildDir = os.path.join(sourceDir, "build")
        if os.path.exists(buildDir):
            shutil.rmtree(buildDir)

        runCommand("cmake -S " + sourceDir + " -B " + buildDir + " -DCMAKE_BUILD_TYPE=Release")
        runCommand("cmake --build " + buildDir + " -j")
        runCommand("cp " + os.path.join(buildDir, "BubbleFinder") + " " + installBinary)
        runCommand("chmod +x " + installBinary)

        os.chdir(oldDirectory)

    print("BubbleFinder installed at " + installBinary)



# Initialize build directory (clean start)
initializeBuildDirectory()

# Install non-Rust prerequisites first
installAptPackages() 
installSeqan()
installPybind11() 
installSpoa()

# Setup Rust toolchain with pinned nightly version (must be done before Rust libraries)
setupRustToolchain()

def installAbpoa():
    print("Installing abPOA (Shared Library)...")
    
    if os.path.exists("/usr/local/lib/libabpoa.so"):
        print("libabpoa.so found in /usr/local/lib. Skipping installation.")
        return

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building abPOA using temporary directory", temporaryDirectory)
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        runCommand("git clone https://github.com/yangao07/abPOA.git")
        os.chdir("abPOA")

        # Patch Makefile as requested by user
        # 1. Comment out SIMDE flags
        runCommand("sed -i 's/^EXTRA_FLAGS += -DUSE_SIMDE/# EXTRA_FLAGS += -DUSE_SIMDE/' Makefile")
        
        # 2. Add -fpic -shared to CFLAGS
        # We append it to the definition of CFLAGS
        runCommand("sed -i 's/^CFLAGS =/CFLAGS += -fpic -shared /' Makefile")

        # 3. Change library name to .so and build command
        # Replace .a with .so
        runCommand("sed -i 's/libabpoa.a/libabpoa.so/g' Makefile")
        # Replace ar command with compiler for shared lib
        runCommand("sed -i 's/ar -csru/$(CC) -shared -o/g' Makefile")

        runCommand("make")

        # Install
        print("Installing abPOA to /usr/local...")
        if os.path.exists("lib/libabpoa.so"):
            runCommand("sudo cp lib/libabpoa.so /usr/local/lib/")
            # Also copy headers
            if not os.path.exists("/usr/include/abpoa"):
                runCommand("sudo mkdir -p /usr/include/abpoa")
            # Usually abpoa.h is in include/. Copying content of include/*
            runCommand("sudo cp include/*.h /usr/include/abpoa/")
            # Also copy to /usr/include/abpoa.h for compatibility if needed?
            # shasta code might use <abpoa.h> or <abpoa/abpoa.h>
            runCommand("sudo cp include/abpoa.h /usr/include/")
        else:
            raise Exception("Build failed: lib/libabpoa.so not found")

        os.chdir(oldDirectory)


# Check that any existing Rust libraries were built with the same toolchain
checkRustLibrariesConsistency()

def installShasta2():
    print("Installing shasta2 and dependencies (abPOA)...")
    
    installPath = os.path.join(INCLUDE_DIR, "shasta2")
    libPath = os.path.join(LIB_DIR, "libshasta2.so")
    staticLibPath = os.path.join(LIB_DIR, "libshasta2.a")
    abpoaLibPath = os.path.join(LIB_DIR, "libabpoa.so")
    abpoaStaticLibPath = os.path.join(LIB_DIR, "libabpoa.a")
    
    # Force rebuild: Remove existing shasta2 files to ensure patch is applied
    print("Removing existing shasta2 to force rebuild with patch...")
    if os.path.exists(installPath):
        shutil.rmtree(installPath)
    if os.path.exists(libPath):
        os.remove(libPath)
    if os.path.exists(staticLibPath):
        os.remove(staticLibPath)

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building shasta2 using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/paoloshasta/shasta2.git")
        
        # Cleanup existing build directory to prevent BuildAbpoa.py failure
        home = os.path.expanduser("~")
        shastaBuildDir = home + "/.shasta2Build"
        if os.path.exists(shastaBuildDir):
             print("Cleaning up existing shasta2 build directory: " + shastaBuildDir)
             shutil.rmtree(shastaBuildDir)

        # Install dependencies (This builds abPOA in ~/.shasta2Build/abpoa)
        # It handles the -fPIC shared library build for us!
        runCommand("python3 shasta2/scripts/InstallBuildPrerequisites.py")

        # Copy abPOA headers and library from shasta2 build
        home = os.path.expanduser("~")
        abpoaBuildDir = home + "/.shasta2Build/abpoa"
        
        # Install abPOA Shared
        if os.path.exists(abpoaBuildDir + "/abPOA/lib/libabpoa.so"):
            print("Installing abPOA shared from shasta2 build...")
            runCommand("cp " + abpoaBuildDir + "/abPOA/lib/libabpoa.so " + abpoaLibPath)
        else:
            print("Warning: libabpoa.so not found in shasta2 build directory.")

        # Install abPOA Static
        if os.path.exists(abpoaBuildDir + "/abPOA/lib/libabpoa.a"):
            print("Installing abPOA static from shasta2 build...")
            runCommand("cp " + abpoaBuildDir + "/abPOA/lib/libabpoa.a " + abpoaStaticLibPath)
        else:
            print("Warning: libabpoa.a not found in shasta2 build directory.")

        # Install Headers
        abpoaIncludeDir = os.path.join(INCLUDE_DIR, "abpoa")
        if not os.path.exists(abpoaIncludeDir):
            os.makedirs(abpoaIncludeDir, exist_ok=True)
        runCommand("cp " + abpoaBuildDir + "/abPOA/include/*.h " + abpoaIncludeDir)
        # Also copy abpoa.h to /include/abpoa.h for compatibility
        runCommand("cp " + abpoaBuildDir + "/abPOA/include/abpoa.h " + INCLUDE_DIR)
        
        os.chdir("shasta2")
        

        
        # Build shasta2 library (Python Module + Static Lib)
        # The option -DBUILD_STATIC_LIBRARY=ON was added recently to shasta2.
        if not os.path.exists("build"):
            os.mkdir("build")
        os.chdir("build")

        # Use absolute path for install prefix to avoid ambiguity
        installTmpDir = os.path.abspath("../install_tmp")
        
        # Run cmake with install prefix
        # We need to explicitly include the paths to our locally built dependencies (poasta, astarpa, etc)
        cxx_flags = f"-I{INCLUDE_DIR}/poasta -I{INCLUDE_DIR}/astarpa -I{INCLUDE_DIR}/simd-minimizers"
        runCommand(f"cmake .. -DCMAKE_CXX_FLAGS=\"{cxx_flags}\" -DBUILD_EXECUTABLE=OFF -DBUILD_PYTHON_MODULE=ON -DBUILD_STATIC_LIBRARY=ON -DCMAKE_INSTALL_PREFIX={installTmpDir}")
        runCommand("make -j")
        # Run make install to populate install_tmp
        runCommand("make install")
        
        # Install shasta2 library (Shared)
        print("Installing shasta2 shared to " + LIB_DIR + "...")
        # Check standard install location first, then fallback
        # Shared might be in lib or bin depending on platform/cmake
        possibleSharedNames = [
            os.path.join(installTmpDir, "lib/shasta2.so"), 
            os.path.join(installTmpDir, "bin/shasta2.so"), 
            os.path.join(installTmpDir, "shasta2-install/bin/shasta2.so"),
            "shasta2-install/bin/shasta2.so", # Fallback if prefix ignored
            "PythonModule/shasta2.so", 
            "src/shasta2.so"
        ]
        
        foundShared = False
        for name in possibleSharedNames:
            if os.path.exists(name):
                 runCommand("cp " + name + " " + libPath)
                 foundShared = True
                 break
                 
        if not foundShared:
            print("Warning: shasta2.so not found even after make install.")

        # Install shasta2 library (Static)
        print("Installing shasta2 static to " + LIB_DIR + "...")
        
        # User confirmed DESTINATION is shasta2-install/bin
        possibleStaticNames = [
            os.path.join(installTmpDir, "shasta2-install/bin/shasta2.a"),
            os.path.join(installTmpDir, "bin/shasta2.a"),
            os.path.join(installTmpDir, "lib/shasta2.a"),
            "shasta2-install/bin/shasta2.a", # Fallback if prefix ignored (installed in build)
            "src/shasta2.a"
        ]
        
        foundStatic = False
        sourceStaticLib = ""
        
        for name in possibleStaticNames:
             if os.path.exists(name):
                 sourceStaticLib = name
                 foundStatic = True
                 break
                 
        if foundStatic:
            print(f"Found static library at: {sourceStaticLib}")
            # Rename to libshasta2.a for Dinara
            runCommand("cp " + sourceStaticLib + " " + os.path.join(LIB_DIR, "libshasta2.a"))
        else:
             print("Error: shasta2.a not found in likely locations. Listing files:")
             runCommand("find . -name '*.a'")
             if os.path.exists(installTmpDir):
                 runCommand(f"find {installTmpDir} -name '*.a'")

        
        # Install shasta2 headers
        if not os.path.exists(installPath):
            os.makedirs(installPath, exist_ok=True)
        # Check install_tmp first, then src
        tmpInclude = os.path.join(installTmpDir, "include")
        if os.path.exists(tmpInclude):
             runCommand(f"cp -r {tmpInclude}/* {installPath}")
        else:
             runCommand("cp -r ../src/* " + installPath)
        
        os.chdir(oldDirectory)


def installTheseusLib():
    print("Installing theseus-lib...")

    installPath = os.path.join(INCLUDE_DIR, "theseus")
    libPath = os.path.join(LIB_DIR, "libtheseus.a")

    if os.path.exists(installPath) and os.path.exists(libPath):
        print("theseus headers and library found. Skipping installation.")
        return

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building theseus-lib using temporary directory", temporaryDirectory)

        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)

        # Clone repo
        runCommand("git clone https://github.com/albertjimenezbl/theseus-lib.git")
        os.chdir("theseus-lib")

        # Patch graph.h: POA-built vertices never have their name field set,
        # so print_as_gfa() writes empty names -> invalid GFA (Bandage shows nothing).
        # Fix: fall back to 1-based index when name is empty.
        graphH = "theseus/graph.h"
        with open(graphH, "r") as f:
            src = f.read()

        old_gfa = (
            "            // Print all nodes as segments\n"
            "            for (const auto &vtx : _vertices)\n"
            "            {\n"
            "                gfa_output << \"S\\t\" << vtx.name << \"\\t\" << vtx.value << \"\\n\";\n"
            "            }\n"
            "\n"
            "            // Print all edges as links\n"
            "            for (const auto &vtx : _vertices)\n"
            "            {\n"
            "                // Go through all incoming vertices (with this you cover all possible edges,\n"
            "                // since the graph is directed)\n"
            "                for (const auto &edge : vtx.in_edges)\n"
            "                {\n"
            "                    gfa_output << \"L\\t\" << _vertices[edge.from_vertex].name << \"\\t+\\t\"\n"
            "                        << vtx.name << \"\\t+\\t\"\n"
            "                        << edge.overlap << \"M\\n\";\n"
            "                }\n"
            "            }"
        )
        new_gfa = (
            "            // Build a per-vertex name: use the stored name if present,\n"
            "            // otherwise fall back to the 1-based vertex index.\n"
            "            // (POA-built graphs never set vtx.name, so it is always empty.)\n"
            "            auto vtx_name = [&](size_t idx) -> std::string {\n"
            "                return _vertices[idx].name.empty()\n"
            "                    ? std::to_string(idx + 1)\n"
            "                    : _vertices[idx].name;\n"
            "            };\n"
            "\n"
            "            // Print all nodes as segments (skip empty sentinel nodes)\n"
            "            for (size_t i = 0; i < _vertices.size(); ++i)\n"
            "            {\n"
            "                if (_vertices[i].value.empty()) continue;\n"
            "                gfa_output << \"S\\t\" << vtx_name(i) << \"\\t\" << _vertices[i].value << \"\\n\";\n"
            "            }\n"
            "\n"
            "            // Print all edges as links (skip any link touching an empty sentinel)\n"
            "            for (size_t i = 0; i < _vertices.size(); ++i)\n"
            "            {\n"
            "                if (_vertices[i].value.empty()) continue;\n"
            "                for (const auto &edge : _vertices[i].in_edges)\n"
            "                {\n"
            "                    if (_vertices[edge.from_vertex].value.empty()) continue;\n"
            "                    gfa_output << \"L\\t\" << vtx_name(edge.from_vertex) << \"\\t+\\t\"\n"
            "                        << vtx_name(i) << \"\\t+\\t\"\n"
            "                        << edge.overlap << \"M\\n\";\n"
            "                }\n"
            "            }"
        )
        if old_gfa not in src:
            print("Warning: could not apply theseus GFA patch (source may have changed). Continuing anyway.")
        else:
            src = src.replace(old_gfa, new_gfa)
            with open(graphH, "w") as f:
                f.write(src)
            print("Applied theseus GFA name patch.")

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


def installVg():
    print("Installing vg (variation graph toolkit)...")

    installBinDir = os.path.join(HOME, ".local", "bin")
    installBinary = os.path.join(installBinDir, "vg")

    os.makedirs(installBinDir, exist_ok=True)

    if os.path.exists(installBinary) and os.access(installBinary, os.X_OK):
        print("vg binary found at " + installBinary + ". Skipping installation.")
        return

    # Fetch the latest release asset URL from the GitHub API.
    apiUrl = "https://api.github.com/repos/vgteam/vg/releases/latest"
    req = urllib.request.Request(apiUrl, headers={"Accept": "application/vnd.github+json",
                                                   "User-Agent": "dinara-install"})
    with urllib.request.urlopen(req) as response:
        release = json.loads(response.read().decode())

    tag = release["tag_name"]
    print("Latest vg release: " + tag)

    # The prebuilt static Linux x86_64 binary is always named "vg" in the release assets.
    assetUrl = None
    for asset in release["assets"]:
        if asset["name"] == "vg":
            assetUrl = asset["browser_download_url"]
            break

    if assetUrl is None:
        raise Exception("Could not find the 'vg' static binary asset in release " + tag)

    print("Downloading vg from " + assetUrl + " ...")
    urllib.request.urlretrieve(assetUrl, installBinary)
    os.chmod(installBinary, 0o755)

    print("vg " + tag + " installed at " + installBinary)


def installMinipoa():
    print("Installing minipoa (fast SIMD POA MSA tool)...")

    installBinDir = os.path.join(HOME, ".local", "bin")
    installBinary = os.path.join(installBinDir, "minipoa")

    os.makedirs(installBinDir, exist_ok=True)

    if os.path.exists(installBinary) and os.access(installBinary, os.X_OK):
        print("minipoa binary found at " + installBinary + ". Skipping installation.")
        return

    sourceDir = os.path.join(HOME, "Downloads", "minipoa")

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building minipoa using temporary directory", temporaryDirectory)

        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)

        if os.path.isdir(sourceDir):
            print("Updating existing minipoa source at", sourceDir)
            runCommand("git -C " + sourceDir + " pull --ff-only")
        else:
            if os.path.exists(sourceDir):
                os.remove(sourceDir)
            print("Cloning minipoa source to", sourceDir)
            runCommand("git clone https://github.com/NCl3-lhd/minipoa.git " + sourceDir)

        buildDir = os.path.join(sourceDir, "build")
        if os.path.exists(buildDir):
            shutil.rmtree(buildDir)

        runCommand("cmake -S " + sourceDir + " -B " + buildDir +
                   " -DCMAKE_BUILD_TYPE=Release")
        runCommand("cmake --build " + buildDir + " -j")

        # The binary is placed in the build directory.
        builtBinary = os.path.join(buildDir, "minipoa")
        if not os.path.exists(builtBinary):
            # Some CMake setups put it directly in the source tree build dir.
            raise Exception("minipoa binary not found after build. Check build output.")

        runCommand("cp " + builtBinary + " " + installBinary)
        os.chmod(installBinary, 0o755)

        os.chdir(oldDirectory)

    print("minipoa installed at " + installBinary)


# Install all Rust libraries
installAstarpa()
installPoasta()
installSimdMinimizers()
installBubbleFinder()

# Install shasta2 (and abpoa via shasta2 scripts)
installShasta2()

# Install theseus-lib (C++23 POA / sequence-to-graph aligner)
installTheseusLib()

installVg()
installMinipoa()

# Make sure the newly created libraries are immediately visible to the loader.
# For local install, we don't need ldconfig, but we might need to set LD_LIBRARY_PATH environment variable
# runCommand("sudo ldconfig")

print("Installation of Dinara prerequisites completed successfully.")
