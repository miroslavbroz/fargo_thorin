/** \file types.h

Definition of the structures used in the FARGO code.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

#include <sys/times.h>

/** The boolean type will be used mainly for the variables
    corresponding to the command line switches  */
typedef int     boolean;	
/** Definition of the type 'real' used throughout the code. You can
    force FARGO to work in single precision by redefining real to
    float here. This is an untested feature of FARGO, however, and in
    practice it should be of little use.  */
typedef double	real;

/** Set of two reals. It is used whenever a set of two reals is needed,
and it usually represents a vector in the (x,y) plane (e.g., a %force), but not
only: it is for instance used to store the mass inside and outside the orbit
in MassInOut(). */
struct pair {
  real            x;
  real            y;
};

typedef struct pair Pair;

/** A structure used to store any scalar fied on the computational
domain. In addition to the size (Nrad*Nsec) and field pointer, it also
has a name, which is used by WritePolarGrid () to name appropriately
the output file, automatically. */
struct polargrid {
  int             Nrad; /**< Radial size of the grid, in number of zones */
  int             Nsec; /**< Azimuthal size of the grid, in number of zones */
  real           *Field; /** Pointer to the array of Nrad*Nsec reals (e.g., density, etc.) */
  char           *Name;  /**< Name of the PolarGrid (can be "dens", "vrad", "vtheta" or "label"). */
};

typedef struct polargrid PolarGrid;

#define		YES	1
#define		NO	0
#define		REAL	1
#define		INT	0
#define		STRING  2
#define 	SINE	0
#define		COSINE	1
#define		ABSCISSA	0
#define		ORDINATE	1
#define		HEIGHT		2
#define		INF		0
#define 	SUP		1
#define         GET             0
#define         MARK            1
#define         FREQUENCY       2
#define         COM_DENSITY     0
#define         COM_VRAD        1
#define         COM_VTHETA      2

#define		MAX1D	16384

#define		MAXPLANETS	120

/** The Param structure handles the parameters of the parameter
file. It allows to associate the values found in the file to the
strings defining the parameters, and ultimately to the global
variables associated to these strings. It is used by the function
var(). */
struct param {
  char name[80]; /**< Name of the parameter. This is the (case insensitive) string found in the parameter file */
  int  type;    /**< Type of the parameter (e.g. INT, REAL, or STRING), see var.c */
  char *variable;  /**< A pointer to the corresponding variable */
  int read;	   /**< This variable is set to YES if and only if the parameter has been found in the parameter file  */
  int necessary;   /**< Tell whether defining the parameter is optional or mandatory.  */
};

typedef struct param Param;

/** This structure is used for monitoring CPU time usage. It is used
    only if -t is specified on the command line. */
struct timeprocess {
  char name[80];
  clock_t clicks;
};

typedef struct timeprocess TimeProcess;

/** Contains all the information about a planetary system at a given
    instant in time. 
    #THORIN: 3rd dimension added, acceleration from the disk added. */
struct planetary_system {
  int nb;			/**< Number of planets */
  real *mass;			/**< Masses of the planets */
  real *x;			/**< x-coordinate of the planets */
  real *y;			/**< y-coordinate of the planets */
  real *z;			/**< z-coordinate of the planets */
  real *vx;			/**< x-coordinate of the planets'velocities */
  real *vy;		        /**< y-coordinate of the planets'velocities */
  real *vz;			/**< z-coordinate of the planets'velocities */
  real *ax;			/**< ax-coordinate of the planets' acceleration from the disk */
  real *ay;			/**< ay-coordinate of the planets' acceleration from the disk */
  real *az;			/**< az-coordinate of the planets' acceleration from the disk */
  real *acc;			/**< The planets' accretion times^-1 */
  char **name;			/**< The planets' names */
  boolean *FeelDisk;		/**< For each planet tells if it feels the disk (ie migrates) */
  boolean *FeelOthers;		/**< For each planet tells if it feels
				   the other planet's gravity */
};

typedef struct planetary_system PlanetarySystem;
