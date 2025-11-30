#!/usr/bin/python3

import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

a = dinara.Assembler()
a.accessAlignmentData()
a.accessReadGraph()
a.removeReadGraphBridges(
	maxDistance = int(config['Assembly']['iterative.bridgeRemovalMaxDistance']))

