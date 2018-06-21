#!/usr/bin/gnuplot

load "opacity_LP.plt"
load "opacity_BL.plt"
load "opacity_ZHU.plt"

########################################################################

GM_SI = 1.32712440018e20
AU_SI = 149597870700.0
R_STANDARD = 8.3144598
MMW = 2.4
R_SI = R_STANDARD/(MMW*0.001)
T2SI = GM_SI/R_SI/AU_SI

ADIABIND = 1.4
MOLWEIGTH = 1.0
GASCONST = 1.0
SQRT_ADIABIND_INV = 1.0/sqrt(ADIABIND)
SQRT2PI_INV = 1.0/sqrt(2.0*pi)

MSOL_SI = 1.98855e30
RHO2CGS = MSOL_SI*1000.0/(AU_SI*100.0)**3

Omega(r) = r**(-1.5)
cs(T) = sqrt(ADIABIND*T/MOLWEIGTH*GASCONST)
H(r,T) = cs(T)/Omega(r)*SQRT_ADIABIND_INV

rho(r,Sigma,T) = Sigma/H(r,T)*SQRT2PI_INV
opacity_LP(r,Sigma,T) = kappa_LP(rho(r,Sigma,T)*RHO2CGS,T*T2SI)
opacity_BL(r,Sigma,T) = kappa_BL(rho(r,Sigma,T)*RHO2CGS,T*T2SI)
opacity_ZHU(r,Sigma,T) = kappa_ZHU(rho(r,Sigma,T)*RHO2CGS,T*T2SI)

########################################################################

set term png small
set out "opacity.png"

set xl "r [au]"
set yl "kappa [cm^2 g^-1]"

load "config.plt"
r(i) = Rmin+(Rmax-Rmin)*int(i/Nsec)/Nr

#set yr [0:8]

p \
  "<paste gasdens.cfg gastemper.cfg" u (r($0)):(opacity_ZHU(r($0),$1,$2)) w l lw 3,\
  "<paste gasdens.cfg gastemper.cfg" u (r($0)):(opacity_BL(r($0),$1,$2))  w l lw 2,\
  "<paste gasdens.cfg gastemper.cfg" u (r($0)):(opacity_LP(r($0),$1,$2))  w l lt 0

q


