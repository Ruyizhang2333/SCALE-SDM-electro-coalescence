#!/bin/bash
#PBS -q S
#PBS -l select=2:ncpus=40:mpiprocs=40
#PBS -l walltime=4:00:00
#PBS -N scale-sdm

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}
#module load intel/16.0.1 mpt hdf5/1.8.16 netcdf/4.4.0
module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1

# run
mpiexec_mpt dplace -s1 ./scale-les_init init.conf || exit
mpiexec_mpt dplace -s1 ./scale-les  run.conf || exit
