#!/usr/bin/python3

import dinara

a = dinara.Assembler()

a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessAssemblyGraphEdgeLists()

a.gatherOrientedReadsByAssemblyGraphEdge()
a.writeOrientedReadsByAssemblyGraphEdge()




