#!/usr/bin/gnuplot

GM_SI = 1.32712440018e20
AU_SI = 149597870700.0
R_STANDARD = 8.3144598
MMW = 2.4
R_SI = R_STANDARD/(MMW*0.001)
T2SI = GM_SI/R_SI/AU_SI

T_MRI = 1000

set term png small
set out "gastemper.png"

set xl "r [au]"
set yl "gastemper"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

#set yr [0:0.006]
#set yr [0:0.0006]

p "gastemper.cfg" u (f($0)):1 w l lt 2 lw 3,\
  T_MRI/T2SI w l lt 0

