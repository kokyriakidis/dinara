#!/usr/bin/python3

import dinara
import GetConfig

config = GetConfig.getConfig()

a = dinara.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.computeReadGraphConnectedComponents(
    int(config['ReadGraph']['minComponentSize']))

