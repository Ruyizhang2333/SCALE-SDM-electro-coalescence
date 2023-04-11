#! /bin/sh -x
#PBS -q mid
#PBS -l nodes=4:ppn=20
#PBS -l walltime=48:00:00
#PBS -N scale-sdm

cd ${PBS_O_WORKDIR}
module load intel/2015.6.233 szip/2.1.1 hdf5/1.8.9 netcdf/4.1.3

#define variables
n_proc=$(cat $PBS_NODEFILE | wc -l)

# run
mpirun -np ${n_proc} ./scale-les_init  init.conf  || exit
mpirun -np ${n_proc} ./scale-les  run.conf  || exit

################################################################################
