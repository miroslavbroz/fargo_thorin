#!/bin/sh

# when exiting the script, clean the scratch
#trap 'clean_scratch' EXIT
# if the task is terminated, copy important files first and clean the scratch after that
#trap 'cp -r $SCRATCHDIR/out $SCRATCHDIR/thorin.out $SCRATCHDIR/thorin.err $DATADIR && clean_scratch' TERM

DATADIR="/storage/praha1/home/$LOGNAME/fargo_study/opacity_Zhu12_test1_OLD_BL94"

cp -r \
  $DATADIR/job8.sh \
  $DATADIR/in \
  $DATADIR/out \
  $DATADIR/src_main \
  $DATADIR/src_reb \
  $SCRATCHDIR

module add intelcdk-17.1
module add openmpi-2.0.1-intel
#module add openmpi-2.0.1-gcc
#module add gcc-4.8.4
#module add cuda
#module add mvapich-3.1.4-gcc
#module add openmpi
#nvidia-smi > fargo.out

cd $SCRATCHDIR/src_main
./makeclean.sh
./make_para.sh

export LD_LIBRARY_PATH=$SCRATCHDIR/src_main:$LD_LIBRARY_PATH

cd $SCRATCHDIR
mkdir out

module list > thorin.out
ls -l * >> thorin.out
mpirun ./thorin -vm in/in.par >> thorin.out 2> thorin.err

cp -r $SCRATCHDIR/out $SCRATCHDIR/thorin.out $SCRATCHDIR/thorin.err $DATADIR || export CLEAN_SCRATCH=false


