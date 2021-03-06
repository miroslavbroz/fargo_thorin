/** \file var.c

Contains the function that connects the string of the parameter file to global variables.
The var() function is found in Interpret.c

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#define __LOCAL
#include "fargo.h"
#undef __LOCAL

void
InitVariables()
{
  var("DT", &DT, REAL, YES, "1.");
  var("SIGMA0", &SIGMA0, REAL, YES, "173.");
  var("NINTERM", &NINTERM, INT, YES, "10.");
  var("NTOT", &NTOT, INT, YES, "1501.");
  var("OUTPUTDIR", OUTPUTDIR, STRING, YES, "out");
  var("INPUTDIR", INPUTDIR, STRING, NO, "in");
  var("INNERBOUNDARY", INNERBOUNDARY, STRING, NO, "WALL");
  var("OUTERBOUNDARY", OUTERBOUNDARY, STRING, NO, "WALL");
  var("LABELADVECTION", ADVLABEL, STRING, NO, "NO");
  var("TRANSPORT", TRANSPORT, STRING, NO, "FAST");
  var("PLANETCONFIG", PLANETCONFIG, STRING, NO, "Systems/SolarSystem.cfg");
  var("MASSTAPER", &MASSTAPER, REAL, NO, "0.0000001");
  var("RADIALSPACING", GRIDSPACING, STRING, NO, "ARITHMETIC");
  var("NRAD", &NRAD, INT, YES, "64.0");
  var("NSEC", &NSEC, INT, YES, "64.0");
  var("RMIN", &RMIN, REAL, YES, "1.0");
  var("RMAX", &RMAX, REAL, YES, "1.0");
  var("THICKNESSSMOOTHING", &THICKNESSSMOOTHING, REAL, NO, "0.0");
  var("ROCHESMOOTHING", &ROCHESMOOTHING, REAL, NO, "0.0");
  var("ASPECTRATIO", &ASPECTRATIO, REAL, YES, "0.05");
  var("VISCOSITY", &VISCOSITY, REAL, NO, "0.0");  
  var("ALPHAVISCOSITY", &ALPHAVISCOSITY, REAL, NO, "0.0");  
  var("ALPHAFLOCK", ALPHAFLOCK, STRING, NO, "NO");
  var("TMRI", &TMRI, REAL, NO, "1000.0");
  var("TWIDTH", &TWIDTH, REAL, NO, "25.0");
  var("ALPHAMRI", &ALPHAMRI, REAL, NO, "0.0");
  var("ALPHADEAD", &ALPHADEAD, REAL, NO, "0.0");
  var("SIGMASLOPE", &SIGMASLOPE, REAL, YES, "0.0");  
  var("RELEASERADIUS", &RELEASERADIUS, REAL, NO, "0.0");  
  var("RELEASEDATE", &RELEASEDATE, REAL, NO, "0.0");  
  var("OMEGAFRAME", &OMEGAFRAME, REAL, NO, "0.0");
  var("DISK", DISK, STRING, NO, "YES");
  var("FRAME", FRAME, STRING, NO, "FIXED");
  var("OUTERSOURCEMASS", OUTERSOURCEMASS, STRING, NO, "NO");
  var("WRITEDENSITY", WRITEDENSITY, STRING, NO, "YES");
  var("WRITEVELOCITY", WRITEVELOCITY, STRING, NO, "YES");
  var("INDIRECTTERM", INDIRECTTERM, STRING, NO, "YES");
  var("EXCLUDEHILL", EXCLUDEHILL, STRING, NO, "NO");
  var("IMPOSEDDISKDRIFT", &IMPOSEDDISKDRIFT, REAL, NO, "0.0");
  var("FLARINGINDEX", &FLARINGINDEX, REAL, NO, "0.0");
  var("ECCENTRICITY", &ECCENTRICITY, REAL, NO, "0.0");
  var("CAVITYRADIUS", &CAVITYRADIUS, REAL, NO, "0.0");
  var("CAVITYRATIO", &CAVITYRATIO, REAL, NO, "1.0");
  var("CAVITYWIDTH", &CAVITYWIDTH, REAL, NO, "1.0");
  var("TRANSITIONRADIUS", &TRANSITIONRADIUS, REAL, NO, "0.0");
  var("TRANSITIONRATIO", &TRANSITIONRATIO, REAL, NO, "1.0");
  var("TRANSITIONWIDTH", &TRANSITIONWIDTH, REAL, NO, "1.0");
  var("LAMBDADOUBLING", &LAMBDADOUBLING, REAL, NO, "0.0");
/* #THORIN: new params assoc.w. energy eq. implementation */
  var("ENERGYEQUATION", ENERGYEQUATION, STRING, NO, "NO");
  var("WRITETEMPERATURE", WRITETEMPERATURE, STRING, NO, "NO");
  var("WRITEENERGY", WRITEENERGY, STRING, NO, "NO");
  var("WRITEDIVV", WRITEDIVV, STRING, NO, "NO");
  var("WRITEQPLUS", WRITEQPLUS, STRING, NO, "NO");
  var("WRITEQBALANCE", WRITEQBALANCE, STRING, NO, "NO");
  var("WRITEOPACITY", WRITEOPACITY, STRING, NO, "NO");
  var("ADIABIND", &ADIABIND, REAL, NO, "1.4");
  var("COOLINGTIME", &COOLINGTIME, REAL, NO, "-1.0");
  var("STELLARIRRADIATION", STELLARIRRADIATION, STRING, NO, "NO");
  var("OPACITYDROP", &OPACITYDROP, REAL, NO, "0.6");
  var("EFFECTIVETEMPERATURE", &EFFECTIVETEMPERATURE, REAL, NO, "5656.0");
  var("STELLARRADIUS", &STELLARRADIUS, REAL, NO, "3.0");
  var("DISCALBEDO", &DISCALBEDO, REAL, NO, "0.5");
  var("PARAMETRICOPACITY", &PARAMETRICOPACITY, REAL, NO, "0.0");
/* #THORIN: a new initialisation option */
  var("INITIALIZEFROMFILE", INITIALIZEFROMFILE, STRING, NO, "NO");
  var("DENSINFILE", DENSINFILE, STRING, NO, "in/gasdens.cfg");
  var("VRADINFILE", VRADINFILE, STRING, NO, "in/gasvrad.cfg");
  var("VTHETAINFILE", VTHETAINFILE, STRING, NO, "in/gasvtheta.cfg");
  var("TEMPERINFILE", TEMPERINFILE, STRING, NO, "in/gastemper.cfg");
/* #THORIN: setup of the damping boundary condition */
  var("DAMPTOWARDS", DAMPTOWARDS, STRING, NO, "ZEROVRAD");
  var("DAMPINGRMINFRAC", &DAMPINGRMINFRAC, REAL, NO, "1.25");
  var("DAMPINGRMAXFRAC", &DAMPINGRMAXFRAC, REAL, NO, "0.84");
  var("DAMPINGPERIODFRAC", &DAMPINGPERIODFRAC, REAL, NO, "1.0");
/* #THORIN: new params assoc.w. the REBOUND code implementation */
  var("NOUTELEMENTS", &NOUTELEMENTS, INT, NO, "1");
  var("PLANETARYDENSITY", &PLANETARYDENSITY, REAL, NO, "1.0");
  var("RESOLVECOLLISIONS", RESOLVECOLLISIONS, STRING, NO, "NO");
  var("TARGETNPL", &TARGETNPL, INT, NO, "-1.0");
  var("IAS15PRECISSION", &IAS15PRECISSION, REAL, NO, "1.e-9");
  var("IAS15MINDT", &IAS15MINDT, REAL, NO, "0.0");
  var("DISCARDPLANETS", DISCARDPLANETS, STRING, NO, "NO");
/* #THORIN: disc-planet interaction control */
  var("WRITETORQUEFILES", WRITETORQUEFILES, STRING, NO, "YES");
  var("HILLCUT", &HILLCUT, REAL, NO, "0.8");
  var("VERTICALDAMPING", &VERTICALDAMPING, REAL, NO, "0.1");
  var("PLANETSFEELDISK", PLANETSFEELDISK, STRING, NO, "NO");
  var("ACCRETIONRATE", &ACCRETIONRATE, REAL, NO, "0.0");
  var("GASACCRETIONHEATING", GASACCRETIONHEATING, STRING, NO, "NO");
/* #THORIN: new params assoc.w. pebble accretion */
  var("PEBBLEACCRETION", PEBBLEACCRETION, STRING, NO, "NO");
  var("BACKREACTION", BACKREACTION, STRING, NO, "NO");
  var("ACCRETIONALHEATING", ACCRETIONALHEATING, STRING, NO, "NO");
  var("WRITEETA", WRITEETA, STRING, NO, "NO");
  var("PEBBLEFLUX", &PEBBLEFLUX, REAL, NO, "2.0e-4");
  var("PEBBLEALPHA", &PEBBLEALPHA, REAL, NO, "1.0e-4");
  var("PEBBLECOAGULATION", &PEBBLECOAGULATION, REAL, NO, "0.5");
  var("PEBBLEBULKDENS", &PEBBLEBULKDENS, REAL, NO, "1.0");
  var("SCHMIDTNUMBER", &SCHMIDTNUMBER, REAL, NO, "1.0");
  var("PARTICLEDIFFUSION", PARTICLEDIFFUSION, STRING, NO, "NO");
  var("HEATINGDELAY", &HEATINGDELAY, INT, NO, "100");
/* #THORIN: tools */
  var("PARAMETRICACCRETION", &PARAMETRICACCRETION, REAL, NO, "0.0");
  var("TORQUEMAPINFILE", TORQUEMAPINFILE, STRING, NO, "NO");
  var("GETTORQUEFORPLANET", &GETTORQUEFORPLANET, INT, NO, "-1");
  var("PEBBLETOGASMAX", &PEBBLETOGASMAX, REAL, NO, "0.01");
  var("FRAGMENTFACTOR", &FRAGMENTFACTOR, REAL, NO, "1.0");
  var("FRAGMENTTHRESHOLD", &FRAGMENTTHRESHOLD, REAL, NO, "1.0");
  var("FRAGMENTALPHA", &FRAGMENTALPHA, REAL, NO, "1.0e-3");
  var("OPACITYTABLE", OPACITYTABLE, STRING, NO, "BL94");
  var("GASFLUX", &GASFLUX, REAL, NO, "0.0");
  var("DAMPING", DAMPING, STRING, NO, "NO");
  var("EVAPORATIONTEMPERATURE", &EVAPORATIONTEMPERATURE, REAL, NO, "1000.0");
  var("EVAPORATIONRATE", &EVAPORATIONRATE, REAL, NO, "1.0e-3");
  var("FIRSTSTEP", &FIRSTSTEP, REAL, NO, "-1.0");
  var("MINIMUMSTEP", &MINIMUMSTEP, REAL, NO, "1.0e-3");
  var("MAXIMUMSTEP", &MAXIMUMSTEP, REAL, NO, "1.0e6");
  var("AERODYNAMICDRAG", AERODYNAMICDRAG, STRING, NO, "NO");
  var("ELEMENTSWITHDISK", ELEMENTSWITHDISK, STRING, NO, "NO");
  var("PEBBLEGRAVITY", PEBBLEGRAVITY, STRING, NO, "NO");
}

