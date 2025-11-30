#!/usr/bin/python3

"""
This checks that all built-in configurations work with
the current dinara executable.
It must run from the dinara-install/bin directory.
"""

import dinara
import os

badCount = 0;
for p in dinara.configurationTable:
    configurationName = p[0]
    command = './dinara --command listConfiguration --config ' + configurationName
    returnCode = os.system(command + ' > /dev/null')
    if not (returnCode == 0):
    	print(configurationName, 'does not work.')
    	badCount = badCount + 1

if badCount== 0:
    print('All built-in configurations work.')
else:
    print(badCount, 'built-in configurations are not working.')



