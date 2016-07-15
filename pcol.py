#!/usr/bin/env python
import os,sys

fileh=[i for i in sys.argv[1:] if os.path.isfile(i)]
if fileh: file=open(fileh[0])
else:file=sys.stdin

if sys.argv.count('evens'):
  for jj in file: print " ".join( [jj.split()[i] for i in range(0,len(jj.split()),2) if jj.split()]) 
elif sys.argv.count('odds'):
  for jj in file: print " ".join( [jj.split()[i] for i in [0]+range(1,len(jj.split()),2) if jj.split()]) 
else:
  for jj in file: print " ".join( [jj.split()[i] for i in [int(k)-1 for k in sys.argv[1:] if k.isdigit()] if jj.split()]) 

