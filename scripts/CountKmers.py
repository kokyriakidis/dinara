#!/usr/bin/python3

import dinara

dinara.openPerformanceLog('CountKmers.log')

a = dinara.Assembler()
a.accessMarkers()
a.countKmers()

