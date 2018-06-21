#!/usr/bin/gnuplot

set term png small
set out "gasvtheta.png"

set xl "r [au]"
set yl "gasvtheta - v_kepl"

load "config.plt"

r(i) = Rmin+(Rmax-Rmin)*int(i/Nsec)/Nr
v_kepl(r) = 1./sqrt(r)

#tmp=0.005
#set yr [-tmp:tmp]
set zeroaxis

set arrow from 1,graph 0 to 1,graph 1 nohead lt 0

p "gasvtheta.cfg" u (r($0)):($1-v_kepl(r($0))) w l lt 3 lw 3,\


