#!/bin/sh

qsub -q default -l select=2:ncpus=14:ngpus=0:mem=2gb:scratch_shared=10gb:cl_doom=True -l walltime=01:00:00 -m abe -M mira@sirrah.troja.mff.cuni.cz job8.sh


