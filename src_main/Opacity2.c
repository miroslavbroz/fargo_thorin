/**
 * @file Opacity2.c
 *
 * @brief More complex opacity tables.
 * Originally, opacity.f Fortran 77 code by D. Semenov,
 * and Energy3D.c code from Fargoca by E. Lega.
 * Substantially modified, cleaned, and corrected.
 *
 */

#include "fargo.h"

#define c_max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define c_min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

int rows = 126, cols = 94;
real *gas_tgrid, *gas_pgrid, *gas_roselland, *gas_planck;

/* Simple linear interpolation */
// Error: There was index overflow, and unnecceseray for-cycle (rewritten)

real interp(real* x, real* y, int n, real x0) {

  if (n <= 1)
    return y[0];
  if (x0 < x[0])
    return y[0];
  if (x0 > x[n-1])
    return y[n-1];

  int i = 1;
  while ((x[i] < x0) && (i < n-1))
    i++;
  return y[i-1] + (y[i]-y[i-1]) * (x0-x[i-1])/(x[i]-x[i-1]);
}

/* Malygin etal. (2014) gas opacity */

real opacity_malygin(int rosseland, real rho, real T_gas) {
  
  real tP, tT, denom, fac, kappa;
  static real mul1, mul2, mul3, mul4;
  static int Plow, Phigh;
  int Tlow, Thigh;
  int i;
  real *array;
  
  if (rosseland) 
    array = &gas_roselland[0];
  else
    array = &gas_planck[0];

  tP = R_STANDARD*1.0e7/MMW * rho * T_gas;
  tT = T_gas;

  if (tT < gas_tgrid[0])
    return 0.0;
  else if (tT > gas_tgrid[cols-1])
    tT = gas_tgrid[cols-1];
  if (tP < gas_pgrid[0])
    tP = gas_pgrid[0];
  else if (tP > gas_pgrid[rows-1])
    tP = gas_pgrid[rows-1];

  for (i=0; i<cols-2; i++) {  // Error: index overflow (-2 was missing)
    if (tT >= gas_tgrid[i]) {
      Tlow = i;
      Thigh = i+1;
    }
  }
  
  fac = ((real) rows) / log10( gas_pgrid[rows-1]/gas_pgrid[0] );
  Plow = (int) (fac * log10(tP/gas_pgrid[0]));
  Plow = c_min(c_max(Plow, 0), rows-2); // Error: index underflow or overflow (check was missing)
  Phigh = Plow + 1;

  // Interpolate on pressure-temperature grid
  denom = 1.0/((gas_tgrid[Thigh] - gas_tgrid[Tlow]) * (gas_pgrid[Phigh] - gas_pgrid[Plow]));
  mul1 = (gas_tgrid[Thigh] - tT) * (gas_pgrid[Phigh] - tP);
  mul2 = (tT - gas_tgrid[Tlow])  * (gas_pgrid[Phigh] - tP);
  mul3 = (gas_tgrid[Thigh] - tT) * (tP - gas_pgrid[Plow]);
  mul4 = (tT - gas_tgrid[Tlow])  * (tP - gas_pgrid[Plow]);

  kappa = denom * (mul1 * array[Tlow + Plow*cols] + mul2 * array[Thigh + Plow*cols] 
    + mul3 * array[Tlow + Phigh*cols] + mul4 * array[Thigh + Phigh*cols]);

  return kappa;
}

/* Semenov etal. (2003) dust opacity */

real opacity_semenov_malygin(int rosseland, real rho, real T_dust) 
{
  const static real dust_rosseland[6][6] = {
    { 0.244153190819473e-03,-0.191588103585833e-03, 0.229132294521843e-03, 0.454088516995579e-06,-0.389280273285906e-08, 0.319006401785413e-11},
    { 0.292517586013911e+01,-0.801305197566973e-01, 0.923001367150898e-03,-0.349287851393541e-05, 0.587151269508960e-08,-0.371333580736933e-11},
    {-0.125893536782728e+02, 0.160197849829283e+00,-0.582601311634824e-03, 0.110792133181219e-05,-0.105506373600323e-08, 0.404997080931974e-12},
    {-0.192550086994197e+01, 0.354301116460647e-01,-0.127355043519226e-03, 0.233773506650911e-06,-0.207587683535083e-09, 0.728005983960845e-13},
    { 0.165244978116957e+01,-0.260479348963077e-02, 0.590717655258634e-05,-0.300202906476749e-08, 0.803836688553263e-12,-0.916709776405312e-16},
    {0.00,0.00,0.00,0.00,0.00,0.00}
  };

  const static real dust_planck[6][6] = {
    { 0.423062838413742e-03,-0.191556936842122e-02, 0.130732588474620e-02,-0.108371821805323e-04, 0.427373694560880e-07,-0.664969066656727e-10},
    {-0.173587657890234e+01, 0.302186734201724e-01, 0.311604311624702e-03,-0.210704968520099e-05, 0.472321713029246e-08,-0.371134628305626e-11},
    {-0.638046050114383e+01, 0.120954274022502e+00,-0.436299138970822e-03, 0.784339065020565e-06,-0.681138400996940e-09, 0.233228177925061e-12},
    {-0.345042341508906e+01, 0.664915248499724e-01,-0.240971082612604e-03, 0.417313950448598e-06,-0.349467552126090e-09, 0.115933380913977e-12},
    { 0.667823366719512e+01,-0.166789974601131e-01, 0.238329845360675e-04,-0.140583908595144e-07, 0.417753047583439e-11,-0.503398974713655e-15}, 
    {0.00,0.00,0.00,0.00,0.00,0.00}
  };

  // evaporation temperature T_ev(rho)
  //       water ice        iron             orthopyroxene    olivine
  const static real logtt[8][4] = {
    {log10(1.090E+02),log10(8.350E+02),log10(9.020E+02),log10(9.290E+02)},
    {log10(1.180E+02),log10(9.080E+02),log10(9.800E+02),log10(9.970E+02)},
    {log10(1.290E+02),log10(9.940E+02),log10(1.049E+03),log10(1.076E+03)},
    {log10(1.430E+02),log10(1.100E+03),log10(1.129E+03),log10(1.168E+03)},
    {log10(1.590E+02),log10(1.230E+03),log10(1.222E+03),log10(1.277E+03)},
    {log10(1.800E+02),log10(1.395E+03),log10(1.331E+03),log10(1.408E+03)},
    {log10(2.070E+02),log10(1.612E+03),log10(1.462E+03),log10(1.570E+03)},
    {log10(2.440E+02),log10(1.908E+03),log10(1.621E+03),log10(1.774E+03)}
  };

  // corresponding densities rho
  const static real logrr[8] = {-18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0};
      
  real *eD;
  real aKext;
  real T_ev[4], temp[8], dT[5], tmax1, tmax2, tmin, T1, T2, TD;
  real T[5], aKrL, aKrR, AA, BB, FF, kappa_gas, kappa_dust;
  real logrho;
  real T_gas;
  int i, j;
  int KK;
  int smooth;

  // T must be more than few [K]
  if (T_dust <= 0.9) 
    return 0.9;

  aKext = 0.0;

  if (rosseland) 
    eD = &dust_rosseland[0];
  else
    eD = &dust_planck[0];
      
  // T_ev interpolation for given rho
  logrho = log10(rho);
  for (i = 0; i < 4 ; i++) {
    for (j = 0; j < 8 ; j++) {
      temp[j] = logtt[j][i];
    }
    // Error: this is not bilinear at all!
    // Error: should be logarithmic (otherwise, there are uneccessary steps)!
    T_ev[i] = pow(10.0, interp(&logrr[0], &temp[0], 8, logrho));
  }
      
  // additional transition temperatures T
  T[0] = T_ev[0];  // water ice
  T[1] = 275.0;    // volatile organics
  T[2] = 425.0;    // refractory organics
  T[3] = 680.0;    // troilite
  tmax1 = c_max(T_ev[1], T_ev[2]);
  tmax2 = c_max(T_ev[2], T_ev[3]);
  tmin  = c_min(tmax1, tmax2);
  T[4] = tmin;     // iron, orthopyroxene, and olivine 

  // corresponding smoothing intervals
  dT[0] = 5.0;    // water ice
  dT[1] = 5.0;    // volatile organics
  dT[2] = 15.0;   // refractory organics
  dT[3] = 5.0;    // troilite
  dT[4] = 100.0;  // iron, orthopyroxen, and olivine

  // determine temperature regime
  KK = 5;
  if (T_dust <= (T[0]+dT[0]))
    KK = 0;    
  for (i = 1; i < 5; i++)
    if ((T_dust > (T[i-1]+dT[i-1])) && (T_dust <= (T[i]+dT[i]))) 
      KK = i;

  if (KK == 5) { // Error: this condition was missing! (nan; cf. dT[KK] below)
    aKext = 0.0;
  } else {
 
    smooth = 0;
    for (i=0; i<5; i++) 
      if (abs(T_dust - T[i]) < dT[i]) // Error: smoothing did NOT work properly when '<=' was used
        smooth = 1; 
 
    if (smooth == 1) {
      T1 = T[KK] - dT[KK];
      T2 = T[KK] + dT[KK];
      TD = T_dust - T[KK];

      // y = a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
      j = (KK+1)*6;
      aKrL = ((((eD[j-1]*T1 + eD[j-2])*T1 + eD[j-3])*T1 + eD[j-4])*T1 + eD[j-5])*T1 + eD[j-6];

      if (KK == 4) {
        // Error: aKg_ext was totally undefined (in the original Semenov F77 code, gas opacities were called)
        aKrR = 0.0;
      } else {
        j = (KK+2)*6;
        aKrR = ((((eD[j-1]*T2 + eD[j-2])*T2 + eD[j-3])*T2 + eD[j-4])*T2 + eD[j-5])*T2 + eD[j-6];
      }

      AA = 0.5e0*(aKrL-aKrR);
      BB = 0.5e0*(aKrL+aKrR);
      FF = PI/2.0/dT[KK];
      aKext = BB-AA*sin(FF*TD);

    } else {  // smooth != 1
      j = (KK+1)*6;
      aKext = ((((eD[j-1]*T_dust + eD[j-2])*T_dust + eD[j-3])*T_dust + eD[j-4])*T_dust + eD[j-5])*T_dust + eD[j-6];
    }
  }  // KK != 5

  kappa_dust = aKext;
  T_gas = T_dust;
  kappa_gas = opacity_malygin(rosseland, rho, T_gas);

  if (rosseland == 0)
    return kappa_dust + kappa_gas;  // Planck opacities are additive
  else
    return c_max(kappa_dust, kappa_gas);  // Rosseland opacities are NOT additive
 
}

/* Semenov etal. (2003) for dust & Malygin etal. (2014) for gas */
/* In our model, we only need Rosseland opacities... */

real opacity_SM03(real rho, real temperature) {
  return opacity_semenov_malygin(1, rho, temperature);
}

/* Read tabulated gas opacities */
/* See also ReadfromFile() in Init.c */

void InitOpacities()
{
  int i, foo=0;
  char name[256];
  FILE *input;

  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &fargostat);
  
  sprintf(name, "%s%s", INPUTDIR, "data_gasopacity_malygin.txt");
  input = fopen (name, "r");

  if (input == NULL) {
    fprintf (stderr, "WARNING ! Opacity file '%s' not found. Stopping.", name); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    prs_exit (1);
  }
  
  fscanf(input, "%i", &rows);
  fscanf(input, "%i", &cols);
  
  gas_pgrid = (real *) malloc(sizeof(real) * rows);
  gas_tgrid = (real *) malloc(sizeof(real) * cols);
  gas_roselland = (real *) malloc(sizeof(real) * rows*cols);
  gas_planck = (real *) malloc(sizeof(real) * rows*cols);
  
  if ((gas_pgrid==NULL) || (gas_tgrid==NULL) || (gas_roselland==NULL) || (gas_planck==NULL)) {
    fprintf (stderr, "Error: Opacities allocation failed.\n"); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    prs_exit (1);
  }
  
  for (i=0; i<rows; i++)
    fscanf(input, "%lf", &gas_pgrid[i]);
  for (i=0; i<cols; i++)
    fscanf(input, "%lf", &gas_tgrid[i]);
  for (i=0; i<rows*cols; i++)
    fscanf(input, "%lf", &gas_planck[i]);
  for (i=0; i<rows*cols; i++)
    fscanf(input, "%lf", &gas_roselland[i]);
  
  fclose (input);
  
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);
  
  mastererr("Opacities were initialized.");
  fflush (stdout);
}

void FreeOpacities() {
  mastererr("Freeing opacities...");
  free(gas_pgrid);
  free(gas_tgrid);
  free(gas_roselland);
  free(gas_planck);
}


