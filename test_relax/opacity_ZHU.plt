
# Zhu etal. (2012), http://www.physics.unlv.edu/~zhzhu/Opacity.html

xlp(rho,T) = log10(rho*T*8.314472*1.e7/2.4)
xlt(T) = log10(T)

kappa_ZHU1(rho,T) = 1.5*(xlt(T)-1.16331)-0.736364					# water ice grains
kappa_ZHU2(rho,T) = -3.53154212*xlt(T)+8.767726-(7.24786-8.767726)*(xlp(rho,T)-5.)/16.	# water ice evaporation
kappa_ZHU3(rho,T) = 1.5*(xlt(T)-2.30713)+0.62 						# non-water grains
kappa_ZHU4(rho,T) = -5.832*xlt(T)+17.7							# graphite corrosion
kappa_ZHU5(rho,T) = 2.129*xlt(T)-5.9398							# no-graphite grains
kappa_ZHU6(rho,T) = 129.88071-42.98075*xlt(T)+(142.996475-129.88071)*0.1*(xlp(rho,T)+4)	# grain evaporation
kappa_ZHU7(rho,T) = -15.0125+4.0625*xlt(T)						# water vapour
kappa_ZHU8(rho,T) = 58.9294-18.4808*xlt(T)+(61.6346-58.9294)*xlp(rho,T)/4.		# water dissociation
kappa_ZHU9(rho,T) = -12.002+2.90477*xlt(T)+(xlp(rho,T)-4)/4.*(13.9953-12.002)		# molecules
kappa_ZHU10(rho,T) = -39.4077+10.1935*xlt(T)+(xlp(rho,T)-4)/2.*(40.1719-39.4077)	# H scattering
kappa_ZHU11(rho,T) = 17.5935-3.3647*xlt(T)+(xlp(rho,T)-6)/2.*(17.5935-15.7376)		# bound-free, free-free
kappa_ZHU12(rho,T) = -0.48								# e- scattering
kappa_ZHU13(rho,T) = 3.586*xlt(T)-16.85							# molecules and H scattering

xlop(rho,T) = \
  xlt(T) < 2.23567+0.01899*(xlp(rho,T)-5.) ? kappa_ZHU1(rho,T) : \
  xlt(T) < 2.30713+0.01899*(xlp(rho,T)-5.) ? kappa_ZHU2(rho,T) : \
  xlt(T) < (17.7-0.62+1.5*2.30713)/(1.5+5.832) ? kappa_ZHU3(rho,T) : \
  xlt(T) < (5.9398+17.7)/(5.832+2.129) ? kappa_ZHU4(rho,T) : \
  xlt(T) < (129.88071+5.9398 + (142.996475-129.88071)*0.1*(xlp(rho,T)+4))/(2.129+42.98075) ? kappa_ZHU5(rho,T) : \
  xlt(T) < (129.88071+15.0125 + (142.996475-129.88071)*0.1*(xlp(rho,T)+4))/(4.0625+42.98075) ? kappa_ZHU6(rho,T) : \
  xlt(T) < 3.28+xlp(rho,T)/4.*0.12 ? kappa_ZHU7(rho,T) : \
  xlt(T) < 3.41+0.03328*xlp(rho,T)/4. ? kappa_ZHU8(rho,T) : \
  xlt(T) < 3.76+(xlp(rho,T)-4)/2.*0.03 ? kappa_ZHU9(rho,T) : \
  xlt(T) < 4.07+(xlp(rho,T)-4)/2.*0.08 ? kappa_ZHU10(rho,T) : \
  xlt(T) < 5.3715+(xlp(rho,T)-6)/2.*0.5594 ? kappa_ZHU11(rho,T) : \
  kappa_ZHU12(rho,T)

xlop_(rho,T) = xlop(rho,T) < 3.586*xlt(T)-16.85 ? ( xlt(T)<4. ? kappa_ZHU13(rho,T) : xlop(rho,T) ) : xlop(rho,T)

kappa_ZHU(rho,T) = 10.**(xlop_(rho,T))

