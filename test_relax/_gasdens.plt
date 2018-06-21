#!/usr/bin/gnuplot

set term x11

set xl "r [au]"
set yl "gasdens"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

#set yr [0:0.0012]

p "gasdens0000.cfg" u (f($0)):1 w l lw 1,\
  "gasdens0100.cfg" u (f($0)):1 w l lw 1,\
  "gasdens0200.cfg" u (f($0)):1 w l lw 1
pa -1


