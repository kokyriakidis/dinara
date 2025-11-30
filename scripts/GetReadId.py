#!/usr/bin/python3

import dinara
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('readName', type=str)
arguments = parser.parse_args()

a = dinara.Assembler()
print(a.getReads().getReadId(arguments.readName))

