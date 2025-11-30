#!/usr/bin/python3

import dinara

import GetConfig
config = GetConfig.getConfig()

a = dinara.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
a.accessMarkerGraphEdges()
a.accessMarkerGraphVertices()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdgeLists()
a.accessAssemblyGraphEdges()

a.analyzeAssemblyGraphBubbles(debug=True)

a.createReadGraphMode1(
    maxAlignmentCount = int(config['ReadGraph']['maxAlignmentCount']))


