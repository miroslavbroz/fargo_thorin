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

ADIABIND = 1.4
MOLWEIGTH = 1.0
GASCONST = 1.0
SQRT_ADIABIND_INV = 1.0/sqrt(ADIABIND)
SQRT2PI_INV = 1.0/sqrt(2.0*pi)

Omega(r) = r**(-1.5)
cs(T) = sqrt(ADIABIND*T/MOLWEIGTH*GASCONST)
H(r,T) = cs(T)/Omega(r)*SQRT_ADIABIND_INV

nu(r,Sigma,T) = alpha(T*T2SI)*cs(T)*H(r,T)
Mdot(r,Sigma,T) = 3.*pi*nu(r,Sigma,T)*Sigma

########################################################################

set term png small
set out "mdot.png"

set xl "r [au]"
set yl "Mdot [code units]"

load "config.plt"
r(i) = Rmin+(Rmax-Rmin)*int(i/Nsec)/Nr

set yr [0:2.5e-8]

p "<paste gasdens.cfg gastemper.cfg" u (r($0)):(Mdot(r($0),$1,$2)) w l lt 8 lw 3


