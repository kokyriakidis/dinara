#!/usr/bin/python3

import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()


# Initialize the assembler and access what we need.
a = dinara.Assembler()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite=True)
a.pruneMarkerGraphStrongSubgraph(
    iterationCount = int(config['MarkerGraph']['pruneIterationCount']))



