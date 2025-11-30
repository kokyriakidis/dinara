#!/usr/bin/python3

import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Create the Assembler.
a = dinara.Assembler()

# Set up the consensus caller.
a.setupConsensusCaller(config['Assembly']['consensusCaller'])

# Access what we need.
a.accessKmers()
a.accessMarkers()
a.accessMarkerGraphVertices()

# Do it.
a.assembleMarkerGraphVertices()



