#!/usr/bin/python3

import dinara

# Initialize the assembler and access what we need.
a = dinara.Assembler()
a.accessMarkers()
a.accessMarkerGraphVertices()
a.accessMarkerGraphEdges()
a.accessMarkerGraphReverseComplementVertex()
a.findMarkerGraphReverseComplementEdges()

