/** \file fargo.h 

Contains all the include directives requested by the code.
In addition, it contains a preprocessor test to know whether
we should build a sequential or MPI executable.

@author THORIN modifications by
Ondřej Chrenko <chrenko@sirrah.troja.mff.cuni.cz>, Copyright (C) 2017;
original code by Frédéric Masset

*/

/* #THORIN: includes rebound.h and omp.h */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "rebound.h"
#include "integrator_ias15.h"
#include "fondam.h"
#include "types.h"
#include "proto.h"
#ifdef _PARALLEL
#include <mpi.h>
#else
#include "mpi_dummy.h"
#endif
#ifndef __LOCAL
#include "param.h"
#include "global_ex.h"
#else
#include "param_noex.h"
#include "global.h"
#endif
#include <sys/times.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>
#include <string.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#ifdef _TRAP_FPE
#include <signal.h>
#include <fenv.h>
#endif
#ifdef OPENMP
#include <omp.h>
#endif
