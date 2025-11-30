#!/usr/bin/python3

import dinara
import GetConfig
import os
import sys

# Check that we have what we need.
if not os.path.lexists('Data'):
    raise Exception('Missing: Data. Use SetupRunDirectory.py to set up the run directory.')
if not os.path.lexists('dinara.conf'):
    raise Exception('Missing: configuration file dinara.conf. Sample available in dinara-install/conf.')


# Read the config file.
config = GetConfig.getConfig()

# Create the assembler.
a = dinara.Assembler(createNew=True)


