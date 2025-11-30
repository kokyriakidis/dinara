#!/usr/bin/python3

import dinara

import GetConfig
config = GetConfig.getConfig()

a = dinara.Assembler()
a.accessMarkers()
a.accessMarkerGraphEdges()
a.accessMarkerGraphVertices()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdgeLists()
a.accessAssemblyGraphEdges()
a.analyzeAssemblyGraphBubbles()


