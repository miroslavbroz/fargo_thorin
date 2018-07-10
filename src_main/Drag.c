/**
 * @file Drag.c
 *
 * @brief Aerodynamic drag on planets (Stokes, not Epstein):
 *
 *          1   S            2
 * \vec a = - C - rho_g v_rel  \hat v_rel
 *          2   m
 *
 * @author Miroslav Broz <miroslav.broz@email.cz>
 *
 */

#include "fargo.h"

#define COEFF 0.48	/* for sphere, low Re */

void StokesDrag(Rho, Vrad, Vtheta, m, x, y, z, vx, vy, vz, axpl, aypl, azpl)
PolarGrid *Rho, *Vrad, *Vtheta;
real m, x, y, z, vx, vy, vz;
real *axpl, *aypl, *azpl;
{
  int nr, ns, i, j, l, lip, ljp;
  real *dens, *vrad, *vtheta, *abs, *ord, *cs;
  real Rperp, ang;
  real xc, yc, vtcell, vrcell, vxcell, vycell, vzcell;
  real H, factor;
  real plrho, plradius, dvx, dvy, dvz, dvsq, tmp;
  real ax, ay, az;

  nr     = Rho->Nrad;
  ns     = Rho->Nsec;
  dens   = Rho->Field;
  vrad   = Vrad->Field;
  vtheta = Vtheta->Field;
  abs    = CellAbscissa->Field;
  ord    = CellOrdinate->Field;
  cs     = SoundSpeed->Field;

  Rperp = sqrt(x*x + y*y);
  if (Rperp >= Rinf[0] && Rperp <= Rsup[nr-1]) {

    i = 0;
    while (Rsup[i] < Rperp) i++;
    ang = atan2(y,x);
    if (ang < 0.0) ang += 2.0*PI;
    j = floor((real)ns/2.0/PI*ang + 0.5);
    if (j==ns) j=0;
    l = j + i*ns;
    lip = l+ns;
    ljp = l+1;
 
    xc = abs[l];
    yc = ord[l];
    vtcell = 0.5*(vtheta[l]+vtheta[ljp])+Rmed[i]*OmegaFrame;
    vrcell = 0.5*(vrad[l]+vrad[lip]);
    vxcell = (vrcell*xc-vtcell*yc)/Rmed[i];
    vycell = (vrcell*yc+vtcell*xc)/Rmed[i];
    vzcell = 0.0;
 
    H = cs[l]*OmegaInv[i]*SQRT_ADIABIND_INV;
    factor = rho_func(z/(sqrt(2.0)*H));
    factor *= 1.0/(sqrt(2.0*PI)*H);
 
    plrho = PLANETARYDENSITY/RHO2CGS;
    plradius = pow(3.0*m/(4.0*PI*plrho), 1.0/3.0);
    dvx = vx-vxcell;
    dvy = vy-vycell;
    dvz = vz-vzcell;
    dvsq = dvx*dvx + dvy*dvy + dvz*dvz;
    tmp = -0.5*COEFF*PI*plradius*plradius/m * factor*dens[l] * dvsq * 1.0/sqrt(dvsq);
    ax = tmp*dvx;
    ay = tmp*dvy;
    az = tmp*dvz;

  } else {
    ax = 0.0;
    ay = 0.0;
    az = 0.0;
  }

  /* other CPUs need this (local) acceleration too... */
  MPI_Allreduce (&ax, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *axpl = tmp;
  MPI_Allreduce (&ay, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *aypl = tmp;
  MPI_Allreduce (&az, &tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *azpl = tmp;
}


void AdvanceSystemFromDrag (Rho, Vrad, Vtheta, sys, dt)
PolarGrid *Rho, *Vrad, *Vtheta;
PlanetarySystem *sys;
real dt;		       
{
  int npl, k;
  real m, x, y, z, vx, vy, vz;
  real ax, ay, az;

  npl = sys->nb;
  for (k = 0; k < npl; k++) {
    if (sys->FeelDisk[k] == YES) {
      m = sys->mass[k];
      x = sys->x[k];
      y = sys->y[k];
      z = sys->z[k];
      vx = sys->vx[k];
      vy = sys->vy[k];
      vz = sys->vz[k];

      StokesDrag (Rho, Vrad, Vtheta, m, x, y, z, vx, vy, vz, &ax, &ay, &az);

      sys->vx[k] += dt * ax;
      sys->vy[k] += dt * ay;
      sys->vz[k] += dt * az;
    }
  }
}


