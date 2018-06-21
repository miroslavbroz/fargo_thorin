/** \file Planet.c

Accretion of disk material onto the planets, and solver of planetary
orbital elements.  The prescription used for the accretion is the one
designed by W. Kley.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

Additional modifications due to 3D planetary orbits
by Miroslav Broz (miroslav.broz@email.cz), Aug 28th 2017


                                             ^ z
              RRoche   .....       z2        |.
                      .     .                | .
                     .   .   .     ZPlanet   |   .
Sun ------------------.     .------          |    . rho =~ exp(-z^2)
        midplane       .....       z1        |   . 
                                             | .
                                             |.
                                             |
*/

#include "fargo.h"

#define QTRAPEPS 1.e-4
#define QTRAPMAX 16

#define MIN(a,b) ((a) < (b) ? (a) : (b));
#define MAX(a,b) ((a) > (b) ? (a) : (b));

void AccreteOntoPlanets (Rho, Vrad, Vtheta, dt, sys)
real dt;
PolarGrid *Rho, *Vrad, *Vtheta;
PlanetarySystem *sys;
{
  real RRoche, Rplanet, distance, dx, dy, deltaM, angle, temp;
  int i_min,i_max, j_min, j_max, i, j, l, jf, ns, nr, lip, ljp, k;
  real Xplanet, Yplanet, Zplanet, Mplanet, VXplanet, VYplanet;
  real facc, facc1, facc2, frac1, frac2; /* We adopt the same notations as W. Kley */
  real *dens, *abs, *ord, *vrad, *vtheta, *cs;
  real PxPlanet, PyPlanet, vrcell, vtcell, vxcell, vycell, xc, yc;
  real dPxPlanet, dPyPlanet, dMplanet;
  real Rperp, dz, z1, z2, H, factor;
  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  vrad   = Vrad->Field;
  vtheta = Vtheta->Field;
  cs     = SoundSpeed->Field;

  for (k=0; k < sys->nb; k++) {
    if (sys->acc[k] > 1e-10) {
      dMplanet = dPxPlanet = dPyPlanet = 0.0;
				/* Hereafter : initialization of W. Kley's parameters */
      facc = fabs(dt)*(sys->acc[k]);
      facc1 = 1.0/3.0*facc;
      facc2 = 2.0/3.0*facc;
      frac1 = 0.75;
      frac2 = 0.45;
				/* W. Kley's parameters initialization finished */
      if ((facc1 > 1.0) || (facc2 > 1.0)) {
        printf("Warning: Gas accretion too fast! facc1 = %.8e, facc2 = %.8e\n", facc1, facc2);
      }
      facc1 = MIN(facc1,1.0);
      facc2 = MIN(facc2,1.0);
      Xplanet = sys->x[k];
      Yplanet = sys->y[k];
      Zplanet = sys->z[k];
      VXplanet = sys->vx[k];
      VYplanet = sys->vy[k];
      Mplanet = sys->mass[k];
      Rperp = sqrt(Xplanet*Xplanet+Yplanet*Yplanet);
      Rplanet = sqrt(Xplanet*Xplanet+Yplanet*Yplanet+Zplanet*Zplanet);
      RRoche = pow((1.0/3.0*Mplanet*MassTaper),(1.0/3.0))*Rplanet; /* Central mass is 1.0; use MassTaper only here! */
      i_min=0;
      i_max=nr-1;
      while ((Rsup[i_min] < Rperp-RRoche) && (i_min < nr)) i_min++;
      while ((Rinf[i_max] > Rperp+RRoche) && (i_max > 0)) i_max--;
      angle = atan2 (Yplanet, Xplanet);
      j_min =(int)((real)ns/2.0/PI*(angle - 2.0*RRoche/Rperp));
      j_max =(int)((real)ns/2.0/PI*(angle + 2.0*RRoche/Rperp));
      PxPlanet = Mplanet*VXplanet;
      PyPlanet = Mplanet*VYplanet;

/* #THORIN 'shared' & 'atomic' construct is now replaced with 'reduction' construct */
#pragma omp parallel for private(j,jf,vrcell,vtcell,vxcell,vycell,l,lip,ljp,xc,yc,dx,dy,distance,deltaM) reduction(+ : dPxPlanet, dPyPlanet, dMplanet)
      for (i = i_min; i <= i_max; i++) {
	for (j = j_min; j <= j_max; j++) {
	  jf = j;
	  while (jf <  0)  jf += ns;
	  while (jf >= ns) jf -= ns;
	  l   = jf+i*ns;
	  lip = l+ns;
	  ljp = l+1;
	  if (jf == ns-1) ljp = i*ns;
	  xc = abs[l];
	  yc = ord[l];
	  dx = Xplanet-xc;
	  dy = Yplanet-yc;
	  distance = sqrt(dx*dx+dy*dy);
	  vtcell=0.5*(vtheta[l]+vtheta[ljp])+Rmed[i]*OmegaFrame;
	  vrcell=0.5*(vrad[l]+vrad[lip]);
	  vxcell=(vrcell*xc-vtcell*yc)/Rmed[i];
	  vycell=(vrcell*yc+vtcell*xc)/Rmed[i];

          H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;

	  if (distance < frac1*RRoche) {
            dz = sqrt(pow(frac1*RRoche,2)-pow(distance,2));
            z1 = Zplanet - dz;
            z2 = Zplanet + dz;
            factor = qtrap(rho_func, z1/(sqrt(2.0)*H), z2/(sqrt(2.0)*H), QTRAPEPS);
            factor *= 1.0/(sqrt(2.0*PI)*H);

	    deltaM = facc1*factor*dens[l]*Surf[i];
	    if (i < Zero_or_active) deltaM = 0.0;
	    if (i >= Max_or_active) deltaM = 0.0;
	    dens[l] *= (1.0 - facc1*factor);
	    dPxPlanet    += deltaM*vxcell;	/* #THORIN 'atomic' directives used to be here */
	    dPyPlanet    += deltaM*vycell;
	    dMplanet     += deltaM;
	  }
	  if (distance < frac2*RRoche) {
            dz = sqrt(pow(frac2*RRoche,2)-pow(distance,2));
            z1 = Zplanet - dz;
            z2 = Zplanet + dz;
            factor = qtrap(rho_func, z1/(sqrt(2.0)*H), z2/(sqrt(2.0)*H), QTRAPEPS);
            factor *= 1.0/(sqrt(2.0*PI)*H);

	    deltaM = facc2*factor*dens[l]*Surf[i];
	    if (i < Zero_or_active) deltaM = 0.0;
	    if (i >= Max_or_active) deltaM = 0.0;
	    dens[l] *= (1.0 - facc2);
	    dPxPlanet    += deltaM*vxcell;
	    dPyPlanet    += deltaM*vycell;
	    dMplanet     += deltaM;
	  }
	}
      }
      MPI_Allreduce (&dMplanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dMplanet = temp;
      MPI_Allreduce (&dPxPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dPxPlanet = temp;
      MPI_Allreduce (&dPyPlanet, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      dPyPlanet = temp;

      if (GasAccretHeating) {
        Heating(Mplanet, dMplanet, Xplanet, Yplanet, dt, nr, ns, k);
      }

      PxPlanet += dPxPlanet;
      PyPlanet += dPyPlanet;
      Mplanet  += dMplanet;
      if (sys->FeelDisk[k] == YES) {
	sys->vx[k] = PxPlanet/Mplanet;
	sys->vy[k] = PyPlanet/Mplanet;
      }
      sys->mass[k] = Mplanet;
    }
  }
}

/* Equlibrium vertical density profile rho(z) */

real rho_func(real z) {
  return exp(-pow(z,2));
}

/* A driver for Simpson-rule integration (as in NR) */

#define FUNC(x) ((*func)(x))

real trapzd(real (*func)(real), real a, real b, int n) {
  real x, tnm, sum, del;
  static real s;
  int it, j;

  if (n==1) {
    return (s = 0.5*(b-a)*(FUNC(a)+FUNC(b)));
  } else {
    for (it=1, j=1; j<n-1; j++) it <<= 1;
    tnm = it;
    del = (b-a)/tnm;
    x = a+0.5*del;
    for (sum=0.0, j=1; j<=it; j++, x+=del) sum += FUNC(x);
    s = 0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

/* Integration with required relative precision (eps) */

real qtrap(real (*func)(real), real a, real b, real eps)
{
  int j;
  real s, olds=0.0;

  for (j=1; j<=QTRAPMAX; j++) {
    s = trapzd(func, a, b, j);
    if (j > 5)
      if (fabs(s-olds) < eps*fabs(olds) || (s == 0.0 && olds == 0.0)) return s;
    olds = s;
  }
  printf("qtrap(): Warning: too many integration steps! QTRAPMAX = %d, eps = %.8e\n", QTRAPMAX, eps);
  return s;
}

/*

Accretion heating (can be applied for gas, pebbles, or both).
Heat only 1 cell containing the planet.

Note: heatsrc[] and heatsrc_index[] are globals, to be used in EnergySource.c

*/

void Heating1(Mplanet, dMplanet, Xplanet, Yplanet, dt, nr, ns, k)
real Mplanet, dMplanet, Xplanet, Yplanet, dt;
int nr, ns, k;
{
  int i, j, l;
  real Rperp, tdelay, taper, plrho, plradius, L, ang, dE;

  Rperp = sqrt(Xplanet*Xplanet + Yplanet*Yplanet);
  if (Rperp < Rinf[0] || Rperp > Rsup[nr-1]) return;
  tdelay = 4.0*DT*(real)HEATINGDELAY;
  taper = PhysicalTime/tdelay;
  if (taper > 1.0) taper = 1.0;
  plrho = PLANETARYDENSITY/RHO2CGS;
  plradius = pow (3.0*Mplanet/(4.0*PI*plrho), 1.0/3.0);
  L = Mplanet*dMplanet/dt/plradius;  // see Benitez et al. (2015); G=1; and we work in terms of LUMINOSITY
  ang = atan2(Yplanet,Xplanet);
  if (ang < 0.0) ang += 2.0*PI;
  i = 0;
  while (Rsup[i] < Rperp) i++;
  j = floor((real)ns/2.0/PI*ang + 0.5);
  if (j==ns) j=0;
  l = j + i*ns;
  heatsrc_index[k] = l;
  dE = L*taper/Surf[i];  // final dimension: W/m^2 -> average power density in cells closest to the planet
  heatsrc[k] += dE;
}

/* Which heating to use, btw.? */

void Heating(Mplanet, dMplanet, Xplanet, Yplanet, dt, nr, ns, k)
real Mplanet, dMplanet, Xplanet, Yplanet, dt;
int nr, ns, k;
{
  Heating1(Mplanet, dMplanet, Xplanet, Yplanet, dt, nr, ns, k);
}


void FindOrbitalElements (x,y,vx,vy,m,n)
real x,y,vx,vy,m;
int n;
{
  real Ax, Ay, e, d, h, a, E, M, V;
  real PerihelionPA;
  FILE *output;
  char name[256];
  if (CPU_Rank != CPU_Number-1) return;
  sprintf (name, "%sorbit%d.dat", OUTPUTDIR, n);
  output = fopenp (name, "a");
  h = x*vy-y*vx;
  d = sqrt(x*x+y*y);
  Ax = x*vy*vy-y*vx*vy-G*m*x/d;
  Ay = y*vx*vx-x*vx*vy-G*m*y/d;
  e = sqrt(Ax*Ax+Ay*Ay)/G/m;
  a = h*h/G/m/(1-e*e);
  if (e != 0.0) {
    E = acos((1.0-d/a)/e);
  } else {
    E = 0.0;
  }
  if ((x*y*(vy*vy-vx*vx)+vx*vy*(x*x-y*y)) < 0) E= -E;
  M = E-e*sin(E);
  if (e != 0.0) {
    V = acos ((a*(1.0-e*e)/d-1.0)/e);
  } else {
    V = 0.0;
  }
  if (E < 0.0) V = -V;
  if (e != 0.0) {
    PerihelionPA=atan2(Ay,Ax);
  } else {
    PerihelionPA=atan2(y,x);
  }
  fprintf (output, "%.12g\t%.12g\t%.12g\t%.12g\t%.12g\t%.12g\n", PhysicalTime, e, a, M, V,\
	   PerihelionPA);
  fclose (output);
}
 
