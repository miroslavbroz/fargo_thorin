
# Keith and Wardle (2014), Table 1, left column (denoted as Bell & Lin 1994)

kappa1_BL(rho,T) = 2.e-4  * rho**(0)     * T**(2)		# ice grains
kappa2_BL(rho,T) = 2.e16  * rho**(0)     * T**(-7)		# ice evaporation
kappa3_BL(rho,T) = 0.1    * rho**(0)     * T**(1./2.)		# metal (!) grains
kappa4_BL(rho,T) = 2.e81  * rho**(1)     * T**(-24)		# metal evaporation
kappa5_BL(rho,T) = 1.e-8  * rho**(2./3.) * T**(3)		# molecules (!)
kappa6_BL(rho,T) = 1.e-36 * rho**(1./3.) * T**(10)		# H scattering
kappa7_BL(rho,T) = 1.5e20 * rho**(1)     * T**(-5./2.)		# bound-free, free-free
kappa8_BL(rho,T) = 0.348  * rho**(0)     * T**(0)		# electron scattering

kappa_BL(rho,T) = T < (2.e16/2.e-4)**(1./9.) ? kappa1_BL(rho,T) : \
  T < (2.e16/0.1)**(1./7.5) ? kappa2_BL(rho,T) : \
  T < 2286.8*rho**(1./24.5) ? kappa3_BL(rho,T) : \
  T < 2029.8*rho**(1./81.) ? kappa4_BL(rho,T) : \
  T < 1.e4*rho**(1./21.) ? kappa5_BL(rho,T) : \
  T < 3.11952e4*rho**(4./75.) ? kappa6_BL(rho,T) : \
  T < 1.7939e8*rho**(2./5.) ? kappa7_BL(rho,T) : \
  kappa8_BL(rho,T)

