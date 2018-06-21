#!/usr/bin/gnuplot

set term png small
set out "gaspebbledens.png"

set xl "r [au]"
set yl "gaspebbledens"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

set yr [0:1.5e-6]

p "gaspebbledens.cfg" u (f($0)):1 w l lw 3

q
  "gasdens.cfg" u (f($0)):($1/10.) w l lw 1

