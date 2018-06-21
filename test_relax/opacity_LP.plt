
# Muller & Kley (2012), Table 1 (denoted as Lin & Papaloizou 1985)

kappa1_LP(rho,T) = 2.e-4  * rho**(0)     * T**(2)		# ice grains
kappa2_LP(rho,T) = 2.e16  * rho**(0)     * T**(-7)		# ice evaporation
kappa3_LP(rho,T) = 5.e-3  * rho**(0)     * T**(1)		# silicate grains
kappa4_LP(rho,T) = 2.e34  * rho**(2./3.) * T**(-9)		# silicate evaporation
kappa5_LP(rho,T) = 2.e-8  * rho**(2./3.) * T**(3)		# molecules
kappa6_LP(rho,T) = 1.e-36 * rho**(1./3.) * T**(10)		# H scattering
kappa7_LP(rho,T) = 1.5e20 * rho**(1)     * T**(-5./2.)		# bound-free, free-free

kappa_LP(rho,T) = T < (2.e16/2.e-4)**(1./9.) ? kappa1_LP(rho,T) : \
  T < (2.e16/5.e-3)**(1./8.) ? kappa2_LP(rho,T) : \
  T < 4.6e3*rho**(1./15.) ? kappa3_LP(rho,T) : \
  T < 3000. ? kappa4_LP(rho,T) : \
  T < 1.1e4*rho**(1./21.) ? kappa5_LP(rho,T) : \
  T < 3.e4*rho**(4./75.) ? kappa6_LP(rho,T) : \
  kappa7_LP(rho,T)

