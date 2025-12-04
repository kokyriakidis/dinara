#!/usr/bin/python3

import os
import shutil
import subprocess
import tempfile

def isArm():
    return subprocess.getoutput("uname -p") == "aarch64"

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


def installAstarpa():
    print("Checking for Rust installation...")
    installPackage("curl")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    rustupPath = os.path.expanduser("~/.cargo/bin/rustup")
    cbindgenPath = os.path.expanduser("~/.cargo/bin/cbindgen")
    
    if not os.path.exists(cargoPath):
        print("Rust not found. Installing Rust...")
        runCommand("curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y")
    
    # Switch to nightly as recommended by astarpa docs
    print("Switching to nightly toolchain...")
    runCommand(rustupPath + " install nightly")
    runCommand(rustupPath + " default nightly")

    if not os.path.exists(cbindgenPath):
        print("Installing cbindgen...")
        runCommand(cargoPath + " install cbindgen")
    else:
        print("cbindgen already installed.")

    installPath = "/usr/include/astarpa"
    if os.path.exists(installPath + "/astarpa.h"):
        print("astarpa header found in %s/astarpa.h. Skipping installation." % installPath)
        return

    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building astarpa library using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner.git")
        os.chdir("astar-pairwise-aligner/astarpa-c")
        
        # Patch for recent nightly Rust (unsafe attributes)
        print("Patching astarpa-c for nightly Rust compatibility...")
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


installAptPackages() 
installSeqan()
installPybind11() 
installSpoa()
installAstarpa()

def installSimdMinimizers():
    installPath = "/usr/include/simd-minimizers"
    if os.path.exists(installPath + "/simd_minimizers.h"):
        print("simd-minimizers-c header found in %s/simd_minimizers.h. Skipping installation." % installPath)
        return

    print("Installing simd-minimizers-c...")
    
    cargoPath = os.path.expanduser("~/.cargo/bin/cargo")
    
    # Ensure Rust is installed (installAstarpa might have done it, but check to be safe)
    if not os.path.exists(cargoPath):
        print("Rust not found. Installing Rust...")
        runCommand("curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y")
        
    with tempfile.TemporaryDirectory() as temporaryDirectory:
        print("Building simd-minimizers-c library using temporary directory", temporaryDirectory)
        
        oldDirectory = os.getcwd()
        os.chdir(temporaryDirectory)
        
        # Clone repo
        runCommand("git clone https://github.com/kokyriakidis/simd-minimizers-c.git")
        os.chdir("simd-minimizers-c")
        
        # Build
        # We use target-cpu=native for best performance as requested
        runCommand("RUSTFLAGS='-C target-cpu=native' " + cargoPath + " build --release")
        
        # Install
        print("Installing simd-minimizers-c to /usr/local...")
        if not os.path.exists(installPath):
            runCommand("sudo mkdir -p " + installPath)
        runCommand("sudo cp simd_minimizers.h " + installPath)
        runCommand("sudo cp target/release/libsimd_minimizers_c.so /usr/local/lib/")
        
        os.chdir(oldDirectory)

installSimdMinimizers()
  
# Make sure the newly created libraries are immediately visible to the loader.
runCommand("sudo ldconfig")

print("Installation of Dinara prerequisies completed successfully.")
