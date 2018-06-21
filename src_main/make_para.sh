#!/bin/bash

# FARGO_THORIN compilation shell script
# for a parallel MPI build.
#
# Copyright (C) 2017 Ondřej Chrenko
# email: chrenko@sirrah.troja.mff.cuni.cz
#
# The script works with 'makefile' and 'makefile.reb'.
# The default setup in the makefiles is intended for
# Linux machines with gcc, MPI and OpenMP support.
# We also tested the code on CPU clusters with
# the PGI compilers but this option is not included
# in the public makefiles to keep them simple and
# ready for free use. Different architectures
# were not tested and users must modify the makefiles
# in case they aim to test various compilers.

# The script first issues 'makefile.reb' to build shared
# REBOUND library 'librebound.so'. Next it issues 'makefile'
# to compile the code using the THORIN build.

module add intelcdk-17.1
module add openmpi-2.0.1-intel

make -f makefile.reb
export FARGO_ARCH=INTEL
make BUILD=parallel

