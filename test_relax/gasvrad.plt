#!/usr/bin/gnuplot

set term png small
set out "gasvrad.png"

set xl "r [au]"
set yl "gasvrad"

load "config.plt"

f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

tmp=0.0005
tmp=0.00005
set yr [-tmp:tmp]
set zeroaxis

p "gasvrad.cfg" u (f($0)):1 w l lt 3 lw 3


