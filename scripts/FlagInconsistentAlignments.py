#!/usr/bin/python3

import dinara
import GetConfig

# Read the config file.
config = GetConfig.getConfig()

# Initialize the assembler and access what we need.
a = dinara.Assembler()
a.accessMarkers()
a.accessAlignmentDataReadWrite()
# a.accessCompressedAlignments()
a.accessReadGraphReadWrite()
a.flagInconsistentAlignments(
    triangleErrorThreshold = 200,
    leastSquareErrorThreshold = 200,
    leastSquareMaxDistance = 1)


