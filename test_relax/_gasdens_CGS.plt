#!/usr/bin/gnuplot

AU_SI = 149597870700.0
MSOL_SI = 1.98855e30
SIGMA2CGS = MSOL_SI*1000.0/(AU_SI*100.0)**2

set term x11

set xl "r [au]"
set yl "gasdens [g cm^-2]"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr
g(Sigma) = Sigma*SIGMA2CGS

#set yr [0:1500.]

p "gasdens0000.cfg" u (f($0)):(g($1)) w l lw 1,\
  "gasdens0100.cfg" u (f($0)):(g($1)) w l lw 1,\
  "gasdens0200.cfg" u (f($0)):(g($1)) w l lw 1,\
  750. w l lt 0
pa -1

set term png small
set out "_gasdens_CGS.png"
rep


