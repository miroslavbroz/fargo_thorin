#!/usr/bin/gnuplot

set term png small
set out "gaspebblevrad.png"

set xl "r [au]"
set yl "gaspebblevrad"

load "config.plt"

f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

tmp=0.0005
set yr [-tmp:tmp]
set zeroaxis

p "gaspebblevrad.cfg" u (f($0)):1 w l lt 3 lw 3


