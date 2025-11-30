#!/usr/bin/python3

import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()
        
a = dinara.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(True)
a.accessDisjointSetsHistogram()
a.flagPrimaryMarkerGraphEdges(
    int(config['Assembly']['mode3.minPrimaryCoverage']), 
    int(config['Assembly']['mode3.maxPrimaryCoverage']), 
    0)
 
