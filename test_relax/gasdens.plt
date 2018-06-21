#!/usr/bin/gnuplot

set term png small
set out "gasdens.png"

set xl "r [au]"
set yl "gasdens"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

set yr [0:0.00060]
set yr [0:0.00006]

p "gasdens.cfg" u (f($0)):1 w l lw 3


