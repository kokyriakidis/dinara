#!/usr/bin/python3

import dinara
import GetConfig

config = GetConfig.getConfig()

a = dinara.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphReverseComplementVertex()
a.accessMarkerGraphEdges(True, True)

a.createMarkerGraphSecondaryEdges(
    secondaryEdgeMaxSkip = int(config['MarkerGraph']['secondaryEdges.maxSkip']))


