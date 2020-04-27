#!/usr/bin/python

import mrcfile as mrc
import sys 

if len(sys.argv) < 2:
	print("usage: "+sys.argv[0]+" mrcfile ")
	exit(1)

s = mrc.open(sys.argv[1]).data.shape

print("mrc shape: "+str(s[0])+" "+str(s[1])+" "+str(s[2]))
