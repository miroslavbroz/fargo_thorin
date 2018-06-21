#!/usr/bin/gnuplot

GM_SI = 1.32712440018e20
AU_SI = 149597870700.0
R_STANDARD = 8.3144598
MMW = 2.4
R_SI = R_STANDARD/(MMW*0.001)
T2SI = GM_SI/R_SI/AU_SI

TMRI = 1000.0
TWIDTH = 25.0
ALPHAMRI = 1.9e-2
ALPHADEAD = 1.e-3

alpha(T) = (ALPHAMRI - ALPHADEAD)*(1-tanh((TMRI-T)/TWIDTH))/2. + ALPHADEAD

########################################################################

set term png small
set out "alpha.png"

set xl "r [au]"
set yl "alpha"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

set yr [1.e-4:0.1]
set logscale y

p "gastemper.cfg" u (f($0)):(alpha($1*T2SI)) w l lw 1


