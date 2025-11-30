#!/usr/bin/python3

import dinara
import argparse
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Parse the command line arguments.
parser = argparse.ArgumentParser(description=
    'Run Mode 3 assembly starting from the marker graph.')
            
parser.add_argument(
    "--debug",
    dest="debug",
    action="store_true",
)    
        
arguments = parser.parse_args() 

# Open a performance log.
dinara.openPerformanceLog('Mode3Assembly.log')

# Create the Assembler object and access what we need.
options = dinara.AssemblerOptions('dinara.conf')
a = dinara.Assembler()
a.accessMarkers()

# Run Mode 3 assembly using existing Anchors.
a.mode3Reassembly(0, options.assemblyOptions.mode3Options, arguments.debug)
 
