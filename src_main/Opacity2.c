
#include "fargo.h"

#define c_max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define c_min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

int opacity_gas_rows = 126, opacity_gas_cols = 94;
real *opa_gas_tscale, *opa_gas_pscale, *opa_gas_ross, *opa_gas_planck;

/* Simple linear interpolation */

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

// This function uses the data for *opa_gas_tscale, *opa_gas_pscale, *opa_gas_ross, *opa_gas_planck
// and interpolates them bilinearly on the given grid

real opacity_malygin(int rosseland, real rho, real T_gas) {
  
  real tP, tT, denom, fac, tempKappa;
  static real mul1, mul2, mul3, mul4;
  static int Plow, Phigh;
  int Tlow, Thigh;
  int i;
  real *array;
  
  if (rosseland) 
    array = &opa_gas_ross[0];
  else
    array = &opa_gas_planck[0];

  // Equation of state for pressure, no conversion needed for temperature
  tP = R_STANDARD*1.0e7/MMW * rho * T_gas;
  tT = T_gas;

  // Search for value in T grid
  if (tT < opa_gas_tscale[0])
    return 0.0;
  else if (tT > opa_gas_tscale[opacity_gas_cols-1])
    tT = opa_gas_tscale[opacity_gas_cols-1];
  if (tP < opa_gas_pscale[0])
    tP = opa_gas_pscale[0];
  else if (tP > opa_gas_pscale[opacity_gas_rows-1])
    tP = opa_gas_pscale[opacity_gas_rows-1];

  for (i=0; i<opacity_gas_cols-2; i++) {  // Error: index overflow (-2 was missing)
    if (tT >= opa_gas_tscale[i]) {
      Tlow = i;
      Thigh = i+1;
    }
  }
  
  // Alternative grid search in P:
  fac = ((real) opacity_gas_rows) / log10( opa_gas_pscale[opacity_gas_rows-1]/opa_gas_pscale[0] );
  Plow = (int) (fac * log10(tP/opa_gas_pscale[0]));
  Plow = c_min(c_max(Plow, 0), opacity_gas_rows-2); // Error: index underflow or overflow (check was missing)
  Phigh = Plow + 1;

  // This method is much faster, but should we use a different P-T grid the 7.0 would be wrong
  // Interpolate on pressure-temperature grid
  denom = 1.0/((opa_gas_tscale[Thigh] - opa_gas_tscale[Tlow]) * (opa_gas_pscale[Phigh] - opa_gas_pscale[Plow]));
  mul1 = (opa_gas_tscale[Thigh] - tT) * (opa_gas_pscale[Phigh] - tP);
  mul2 = (tT - opa_gas_tscale[Tlow])  * (opa_gas_pscale[Phigh] - tP);
  mul3 = (opa_gas_tscale[Thigh] - tT) * (tP - opa_gas_pscale[Plow]);
  mul4 = (tT - opa_gas_tscale[Tlow])  * (tP - opa_gas_pscale[Plow]);

  tempKappa = denom * (mul1 * array[Tlow + Plow * opacity_gas_cols] + mul2 * array[Thigh + Plow * opacity_gas_cols] 
    + mul3 * array[Tlow + Phigh * opacity_gas_cols] + mul4 * array[Thigh + Phigh * opacity_gas_cols]);

  return tempKappa;
}

// The Malygin gas-opacities live as data on a grid, thus we create a grid
// The Semenov dust-opacities live as coefficients in a analytic function, thus we need only those
// We want this data to be initialized and a pointer to it already lying around, cuz we'll call the opacity function often

real opacity_semenov_malygin(int rosseland, real rho, real temperature) 
{
  // Rosseland mean opacities
  const static real opa_ross[6][6] = {
    {0.244153190819473e-03,-0.191588103585833e-03,0.229132294521843e-03,0.454088516995579e-06,-0.389280273285906e-08,0.319006401785413e-11},
    {0.292517586013911e+01,-0.801305197566973e-01,0.923001367150898e-03,-0.349287851393541e-05,0.587151269508960e-08,-0.371333580736933e-11},
    {-0.125893536782728e+02,0.160197849829283e+00,-0.582601311634824e-03,0.110792133181219e-05,-0.105506373600323e-08,0.404997080931974e-12},
    {-0.192550086994197e+01,0.354301116460647e-01,-0.127355043519226e-03,0.233773506650911e-06,-0.207587683535083e-09,0.728005983960845e-13},
    {0.165244978116957e+01,-0.260479348963077e-02,0.590717655258634e-05,-0.300202906476749e-08,0.803836688553263e-12,-0.916709776405312e-16},
    {0.00,0.00,0.00,0.00,0.00,0.00}
  };

  // Planck mean opacities
  const static real opa_planck[6][6] = {
    {0.423062838413742e-03,-0.191556936842122e-02,0.130732588474620e-02,-0.108371821805323e-04,0.427373694560880e-07,-0.664969066656727e-10},
    {-0.173587657890234e+01,0.302186734201724e-01,0.311604311624702e-03,-0.210704968520099e-05,0.472321713029246e-08,-0.371134628305626e-11},
    {-0.638046050114383e+01,0.120954274022502e+00,-0.436299138970822e-03,0.784339065020565e-06,-0.681138400996940e-09,0.233228177925061e-12},
    {-0.345042341508906e+01,0.664915248499724e-01,-0.240971082612604e-03,0.417313950448598e-06,-0.349467552126090e-09,0.115933380913977e-12},
    {0.667823366719512e+01,-0.166789974601131e-01,0.238329845360675e-04,-0.140583908595144e-07,0.417753047583439e-11,-0.503398974713655e-15}, 
    {0.00,0.00,0.00,0.00,0.00,0.00}
  };

  // Water ice evaporation temperature, then metallic iron, orthopyroxene and olivine; depending on density array 'ro'
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

  // The density values for which the evaporation temperatures are given
  const static real logro[8] = {-18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0};
      
  // Global variable(s):
  real *eD;
  real aKext;

  // Local variable(s):
  real T_ev[4], temp[8], dT[5], tmax1, tmax2, tmin, T1, T2, TD;
  real T[5], aKrL, aKrR, AA, BB, FF, kappa_gas, kappa_dust;
  real logrho;
  int KK;
  int i, j, iter;
  int smooth;

  // Initialization of the output parameters:
  // T must be more that few [K]
  if (temperature <= 0.9) 
    return 0.9;

  aKext = 0.0;

  if (rosseland) 
    eD = &opa_ross[0];
  else
    eD = &opa_planck[0];
      
  // Interpolation of a matrix of evaporation temperatures for a given density 'rho_in':
  logrho = log10(rho);
  for (i = 0; i < 4 ; i++) {
    for (j = 0; j < 8 ; j++) {
      temp[j] = logtt[j][i];
    }
  ; // Error: this is not bilinear at all!
    // Error: this should be logarithmic (otherwise, there are uneccessary steps)!
    T_ev[i] = pow(10.0, interp(&logro[0], &temp[0], 8, logrho));
  }
      
  // Set up the evaporation temperature array 'T(1:5)':
  T[0] = T_ev[0];  // water ice
  T[1] = 275.0;    // volatile organics
  T[2] = 425.0;    // refractory organics
  T[3] = 680.0;    // troilite
  tmax1 = c_max(T_ev[1], T_ev[2]);
  tmax2 = c_max(T_ev[2], T_ev[3]);
  tmin  = c_min(tmax1, tmax2);
  T[4] = tmin;     // average evaporation temperatures of iron, olivine, and orthopyroxene

  // Determination of a temperature regime where smoothing is necessary 
  dT[0] = 5.0e0;    // an interval of T where ice is evaporating
  dT[1] = 5.0e0;    // volatile organics
  dT[2] = 15.0e0;   // refractory organics
  dT[3] = 5.0e0;    // troilite
  dT[4] = 100.0e0;  // a wide interval of T where iron, pyroxe, and olivine are evaporating

  // Determination of a temperature regime for a given temperature:
  KK = 5;
  if (temperature <= (T[0]+dT[0]))
    KK = 0;    
  for (iter = 1; iter < 5; iter++)
    if ((temperature > ( T[iter-1]+dT[iter-1] )) && (temperature <= ( T[iter]+dT[iter] ))) 
      KK = iter;

  if (KK == 5) { // Error: this condition was missing! (nan; cf. dT[KK] below)
    aKext = 0.0;
  } else {
 
    // The dust-dominated opacity:
    // Switch to smoothing of the opacity if a temperature is near an evaporation temperature of a dust material:
    smooth = 0;
    for (i=0; i<5; i++) 
      if (abs(temperature - T[i]) < dT[i]) // Error: smoothing did NOT work properly when '<=' was used
        smooth = 1; 
 
    // If 'T_in' inside of (Tev-dT, Tev+dT) then start smoothing:
    if (smooth == 1) {
      T1 = T[KK] - dT[KK];
      T2 = T[KK] + dT[KK];
      TD = temperature - T[KK];

      // Calculation of a mean opacity (extinction):
      aKrL = ((((eD[(KK+1)*6-1]*T1 + eD[(KK+1)*6-2])*T1 + eD[(KK+1)*6-3])*T1 + eD[(KK+1)*6-4])*T1 + eD[(KK+1)*6-5])*T1 + eD[(KK+1)*6-6];

      if (KK == 4) {
        // Error: aKg_ext was totally undefined (in the original Semenov code, gas opacities were called)
        aKrR = 0.0;
      } else 
        aKrR = ((((eD[(KK+2)*6-1]*T2 + eD[(KK+2)*6-2])*T2 + eD[(KK+2)*6-3])*T2 + eD[(KK+2)*6-4])*T2 + eD[(KK+2)*6-5])*T2 + eD[(KK+2)*6-6];

      AA = 0.5e0*(aKrL-aKrR);
      BB = 0.5e0*(aKrL+aKrR);
      FF = PI/2.0/dT[KK];
      aKext = BB-AA*sin(FF*TD);

    } else {
      // Smoothing is not necessary, direct calculation by a fit polynom of fifth degree:
      // y = a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
      aKext = ((((eD[(KK+1)*6-1]*temperature + eD[(KK+1)*6-2])*temperature + eD[(KK+1)*6-3])*temperature + eD[(KK+1)*6-4])*temperature + eD[(KK+1)*6-5])*temperature + eD[(KK+1)*6-6];
    }
  }  // KK != 5

  kappa_dust = aKext;
  kappa_gas = opacity_malygin(rosseland, rho, temperature);

  if (rosseland == 0)
    return kappa_dust + kappa_gas;  // Planck opacities are additive
  else
    return c_max(kappa_dust,kappa_gas);  // Rosseland opacities are NOT additive
 
}

/* Semenov (2003) for dust & Malygin etal. (2014) for gas */
/* In our model, we only need Rosseland opacities... */

real opacity_SM03(real rho, real temperature) {
  return opacity_semenov_malygin(1, rho, temperature);
}

// This is essentially just a copy-paste of the ReadfromFile function, specialized for the needs of our opacity data
// * It allocates the small real 5x6 arrays for Semenov Dust
// * the real 126x100 for the Malygin Gas opacities
// * And returns a pointer to the allocated memory on the local RAM

void InitOpacities()
{
  int c, foo=0;
  char dataname[256];
  FILE *datafile;

  /* Simultaneous read access to the same file have been observed to give wrong results. */
  /* A sequential reading is imposed below. */
  /* If current CPU has a predecessor, wait for a message from him */
  if (CPU_Rank > 0) MPI_Recv (&foo, 1, MPI_INT, CPU_Rank-1, 10, MPI_COMM_WORLD, &fargostat);
  
  sprintf(dataname, "%s%s", INPUTDIR, "data_gasopacity_malygin.txt");
  datafile = fopen (dataname, "r");

  if (datafile == NULL) {
    fprintf (stderr, "WARNING ! Opacity file not found. Stopping. Searched for %s\n", dataname); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
    prs_exit (1);
  }
  
  fscanf(datafile, "%i", &opacity_gas_rows);
  fscanf(datafile, "%i", &opacity_gas_cols);
  
  fprintf(stdout,"Initializing opacity data with %i rows and %i cols.\n", opacity_gas_rows, opacity_gas_cols);
  fflush(stdout);

  // opacity_rows & cols: global variables. TO DO if necessary: read their length in dynamically
  // allocate memory for *opa_tscale, *opa_pscale, *opa_ross, *opa_planck;
  opa_gas_pscale = (real *)malloc(sizeof(real) * opacity_gas_rows); // Pressure grid
  opa_gas_tscale = (real *)malloc(sizeof(real) * opacity_gas_cols); // Temperature grid
  
  opa_gas_ross   = (real *)malloc(sizeof(real) * opacity_gas_rows * opacity_gas_cols); // Rosseland opacities
  opa_gas_planck = (real *)malloc(sizeof(real) * opacity_gas_rows * opacity_gas_cols); // Planck opacities
  
  if ((opa_gas_pscale == NULL) || (opa_gas_tscale == NULL) || (opa_gas_ross == NULL) || (opa_gas_planck == NULL)) {
    fprintf (stderr, "WARNING ! Opacity files failed to initialize. \n"); 
    if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
      prs_exit (1);
  }
  
  for (c = 0; c < opacity_gas_rows; c++) {
    fscanf(datafile, "%lf", &opa_gas_pscale[c]);
  }
  for (c = 0; c < opacity_gas_cols; c++) {
    fscanf(datafile, "%lf", &opa_gas_tscale[c]);
  }
  // WHICH ORDER IS CORRECT?! planck -> rosseland?
  for (c = 0; c < opacity_gas_rows * opacity_gas_cols; c++) {
    fscanf(datafile, "%lf", &opa_gas_planck[c]);
  }
  for (c = 0; c < opacity_gas_rows * opacity_gas_cols; c++) {
    fscanf(datafile, "%lf", &opa_gas_ross[c]);
  }
  
  fclose (datafile);
  
  /* Next CPU is waiting. Tell it to start now by sending the message that it expects */
  if (CPU_Rank < CPU_Number-1) MPI_Send (&foo, 1, MPI_INT, CPU_Rank+1, 10, MPI_COMM_WORLD);
  MPI_Barrier (MPI_COMM_WORLD);  /* previous CPUs do not touch anything meanwhile */
  
  mastererr("Finished initializing opacities.");
  fflush (stdout);
}

void FreeOpacities() {
  mastererr("Freeing kappa opacities...");
  free(opa_gas_pscale);
  free(opa_gas_tscale);
  free(opa_gas_ross);
  free(opa_gas_planck);
}


