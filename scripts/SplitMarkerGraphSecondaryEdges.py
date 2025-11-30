#!/usr/bin/python3

import dinara
import GetConfig

config = GetConfig.getConfig()

a = dinara.Assembler()

a.accessMarkers()
a.accessMarkerGraphEdges(True, True)
a.accessMarkerGraphReverseComplementEdge()

a.splitMarkerGraphSecondaryEdges(
    errorRateThreshold = float(config['MarkerGraph']['secondaryEdges.split.errorRateThreshold']),
    minCoverage = int(float(config['MarkerGraph']['secondaryEdges.split.minCoverage'])
    )


