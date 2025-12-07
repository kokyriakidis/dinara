#!/usr/bin/python3

import os
import shutil
import subprocess
import tempfile

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
        ("/usr/local/lib/libastarpa_c.a", "libastarpa_c.a"),
        ("/usr/local/lib/libpoasta_c.a", "libpoasta_c.a"),
        ("/usr/local/lib/libsimd_minimizers_c.a", "libsimd_minimizers_c.a"),
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
            "  sudo rm -rf /usr/include/astarpa /usr/include/poasta /usr/include/simd-minimizers\n" +
            "  sudo rm -f /usr/local/lib/libastarpa_c.a /usr/local/lib/libpoasta_c.* /usr/local/lib/libsimd_minimizers_c.*\n" +
            "  python3 scripts/InstallPrerequisites-Ubuntu.py\n" +
            "="*70
        )
    
    print("All existing Rust libraries were built with the same toolchain (%s)." % firstHash[:12])
    return True

def runCommand(command):
    if(os.system(command)):
        raise Exception("Error running command: " + command)
        
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
    ]
    runCommand("sudo apt-get install --assume-yes " + " ".join(packages))



# We don't use Ubuntu package libseqan2-dev because
# it does not include the fix for this issue:
# https://github.com/seqan/seqan/issues/2524
# Instead, we clone the GitHub seqan/seqan repository, 
# then copy seqan/include/seqan to /usr/include/seqan.
def installSeqan():

    # The path where the include files will go.
    installPath = "/usr/include/seqan"
    
    # First check that this path does not exist.
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
    
    if os.path.exists("/usr/local/include/spoa/spoa.hpp"):
        print("spoa header found in /usr/local/include/spoa/spoa.hpp. Skipping installation.")
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

    installPath = "/usr/include/astarpa"
    libPath = "/usr/local/lib/libastarpa_c.a"
    
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

        # Build
        runCommand(cargoPath + " build --release")
        
        # Generate header
        runCommand(cbindgenPath + " --lang c --cpp-compat --crate astarpa-c -o astarpa.h")
        
        # Install
        print("Installing astarpa to /usr/local...")
        if not os.path.exists(installPath):
            runCommand("sudo mkdir -p " + installPath)
        runCommand("sudo cp astarpa.h " + installPath)
        runCommand("sudo cp ../target/release/libastarpa_c.a /usr/local/lib/")
        
        os.chdir(oldDirectory)


def installPoasta():
    print("Installing poasta-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    
    installPath = "/usr/include/poasta"
    libPath = "/usr/local/lib/libpoasta_c.a"
    
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
        print("Installing poasta-c to /usr/local...")
        if not os.path.exists(installPath):
            runCommand("sudo mkdir -p " + installPath)
        runCommand("sudo cp poasta.h " + installPath)
        runCommand("sudo cp target/release/libpoasta_c.so /usr/local/lib/")
        runCommand("sudo cp target/release/libpoasta_c.a /usr/local/lib/")
        
        os.chdir(oldDirectory)

def installSimdMinimizers():
    print("Installing simd-minimizers-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    
    installPath = "/usr/include/simd-minimizers"
    libPath = "/usr/local/lib/libsimd_minimizers_c.a"
    
    if os.path.exists(installPath + "/simd_minimizers.h") and os.path.exists(libPath):
        print("simd-minimizers-c header found in %s/simd_minimizers.h. Skipping installation." % installPath)
        return
        
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
        print("Installing simd-minimizers-c to /usr/local...")
        if not os.path.exists(installPath):
            runCommand("sudo mkdir -p " + installPath)
        runCommand("sudo cp simd_minimizers.h " + installPath)
        runCommand("sudo cp target/release/libsimd_minimizers_c.so /usr/local/lib/")
        runCommand("sudo cp target/release/libsimd_minimizers_c.a /usr/local/lib/")
        
        os.chdir(oldDirectory)


# Install non-Rust prerequisites first
installAptPackages() 
installSeqan()
installPybind11() 
installSpoa()

# Setup Rust toolchain with pinned nightly version (must be done before Rust libraries)
setupRustToolchain()

# Check that any existing Rust libraries were built with the same toolchain
checkRustLibrariesConsistency()

# Install all Rust libraries (all use the same toolchain to avoid symbol conflicts)
installAstarpa()
installPoasta()
installSimdMinimizers()
  
# Make sure the newly created libraries are immediately visible to the loader.
runCommand("sudo ldconfig")

print("Installation of Dinara prerequisites completed successfully.")
