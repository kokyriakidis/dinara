#!/usr/bin/python3

import dinara

a = dinara.Assembler()
a.accessAssemblyGraphVertices()
a.accessAssemblyGraphEdges()
a.accessAssemblyGraphSequences()
a.writeGfa1('Assembly.gfa')



