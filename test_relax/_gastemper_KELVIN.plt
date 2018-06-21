#!/usr/bin/gnuplot

GM_SI = 1.32712440018e20
AU_SI = 149597870700.0
R_STANDARD = 8.3144598
MMW = 2.4
R_SI = R_STANDARD/(MMW*0.001)
T2SI = GM_SI/R_SI/AU_SI

T_MRI = 1000

set term x11

set xl "r [au]"
set yl "gastemper [K]"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr
g(T) = T*T2SI

set yr [0:1500.]

p "gastemper0000.cfg" u (f($0)):(g($1)) w l lw 1,\
  "gastemper0100.cfg" u (f($0)):(g($1)) w l lw 1,\
  "gastemper0200.cfg" u (f($0)):(g($1)) w l lw 1,\
  T_MRI w l lt 0
pa -1

set term png small
set out "_gastemper_KELVIN.png"
rep


