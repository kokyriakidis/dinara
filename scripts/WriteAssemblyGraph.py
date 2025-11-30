#!/usr/bin/python3

import dinara
import sys

if not len(sys.argv) == 2:
     raise Exception('Call with one argument, the file name.')
fileName = sys.argv[1];

a = dinara.Assembler()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphEdgeLists()
a.writeAssemblyGraph(fileName)



