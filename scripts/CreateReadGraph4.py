#!/usr/bin/python3

import os
import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = dinara.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
a.accessCompressedAlignments()

a.createReadGraph4(
    maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']))


