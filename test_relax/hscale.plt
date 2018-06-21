#!/usr/bin/gnuplot

# temperature[l] = energy[l]*MOLWEIGHT*(ADIABIND-1.0)/(rho[l]*GASCONST);
# energy[l] = temperature[l] / (MOLWEIGHT*(ADIABIND-1.0)) * rho[l]*GASCONST;
# cs[l] = sqrt( ADIABIND*(ADIABIND-1.0)*energy[l]/rho[l] );
# cs[l] = sqrt( ADIABIND*temperature[l]/MOLWEIGHT*GASCONST );
# H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;

ADIABIND = 1.4
MOLWEIGTH = 1.0
GASCONST = 1.0
SQRT_ADIABIND_INV = 1.0/sqrt(ADIABIND)

Omega(r) = r**(-1.5)
cs(T) = sqrt(ADIABIND*T/MOLWEIGTH*GASCONST)
H(r,T) = cs(T)/Omega(r)*SQRT_ADIABIND_INV

########################################################################

set term png small
set out "hscale.png"

set xl "r [au]"
set yl "H"

load "config.plt"
f(i) = Rmin+(Rmax-Rmin)*i/Nsec/Nr

#set yr [0:0.0004]
#set size ratio -1

p "gastemper.cfg" u (f($0)):(H(f($0),$1)) w l lw 1


