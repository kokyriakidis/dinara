#!/usr/bin/python3

import dinara
import GetConfig

config = GetConfig.getConfig()


a = dinara.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.createMarkerGraphEdgesStrict(
    minEdgeCoverage = int(config['MarkerGraph']['minEdgeCoverage']),
    minEdgeCoveragePerStrand = int(config['MarkerGraph']['minEdgeCoveragePerStrand']))

