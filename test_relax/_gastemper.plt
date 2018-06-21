#!/usr/bin/gnuplot

set term x11

set xl "r [au]"
set yl "gastemper"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

#set yr [0:0.0012]

p "gastemper000.cfg" u (f($0)):1 w l lw 1,\
  "gastemper000.cfg_h0.04" u (f($0)):1 w l lw 1,\
  "gastemper150.cfg" u (f($0)):1 w l lw 1
pa -1


