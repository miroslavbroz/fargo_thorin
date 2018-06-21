#!/usr/bin/env python

### An auxiliary python script to convert binary
### output files from the first example calculation into ascii
### input files for the second example calculation.
###
### Author: Ondrej Chrenko
### chrenko@sirrah.troja.mff.cuni.cz
### 2017

import numpy as np
import sys

### Set the grid parameters
Nrad=1024
#Nsec=1
Nsec=4

### Select the files' id
#nout=0
nbegin=0
nout=120

def bin2ascii(filetype, nout):
    filename = "../out_relax/" + filetype + str(nout) + ".dat"
    try:
        fargogrid1D = np.fromfile(filename, dtype=np.double)
    except IOError as ex:
        print ex
        return
    outfile = open("%s%04d.cfg" % (filetype, nout), 'w')
    for i in range (0,Nrad):
        for j in range (0,Nsec): 
            l =  j + i*Nsec
            outfile.write ("%#.18g\n" % fargogrid1D[l])
    outfile.close()

#filetypes = ["gasdens", "gastemper", "gasvrad", "gasvtheta", "gaspebbledens", "gaspebblevrad", "gaspebblevtheta"]
#filetypes = ["gasdens", "gastemper", "gasvrad", "gasvtheta"]
filetypes = ["gasdens", "gastemper"]

for filetype in filetypes:
    for i in xrange(nbegin,nout+1):
        print i
        bin2ascii(filetype, i)


