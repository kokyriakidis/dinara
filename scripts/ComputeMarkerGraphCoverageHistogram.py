#!/usr/bin/python3

import dinara
import GetConfig

config = GetConfig.getConfig()

a = dinara.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges(accessEdgesReadWrite = True)
a.computeMarkerGraphCoverageHistogram()

