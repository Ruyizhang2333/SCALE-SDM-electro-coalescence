Introduction of electro-coalescence treatment of droplets SCALE-SDM 0.2.5-2.3.0 for reproduce the results of Zhang et al. (2023, GMD)

Author: Ruyi Zhang (ruyizhang96@gmail.com)

# Model description

SCALE is a library of weather and climate models of the Earth and planets (Detailed instruction of SCALE (without SDM) are avaiable from the SCALE web page: http://r-ccs-climate.riken.jp/scale/doc/index.html). 

Shima et al. (2009, 2020) constructed a particle-based cloud model SCALE-SDM by implementing the SDM into SCALE. 

The analytic expression for electro-coalescence with the accurate electrostatic force for a pair of droplets with opposite sign charges is established by treating the droplets as conducting spheres (CSs) developed by Zhou et al. (2009) were implmented into SCALE-SDM.

# Simulation setup
## Required software and supported environment

Fortran and C compiler are required to compile SCALE-SDM. MPI, NetCDF4, and HDF5 libraries are are also required.

The numerical experiments were conducted by using Intel Fortran/C compiler 17.0.0, SGI MPT 2.15, HDF5 1.8.12, and NetCDF 4.4.1. For data analysis, R 3.2.3 and python 3.9.0 was used.

## Set environment variable
```ruby
export SCALE_SYS=Linux64-intel-impi
```

## Set compiler options
Modified the compiler options in Makefile according to your computer, here's sample:
```ruby
vi sysdep/Makedef.Linux64-intel-mpich2
...
-NETCDF_INCLUDE ?= -I$(NETCDF3)/include
-NETCDF_LIBS    ?= -L$(NETCDF3)/lib -lnetcdff -lnetcdf
+NETCDF_INCLUDE ?= $(NETCDF_INC) -lnetcdf $(NETCDF_FORTRAN_INC) -lnetcdff
+NETCDF_LIBS    ?= $(NETCDF_LIB) -lnetcdf $(NETCDF_FORTRAN_LIB) -lnetcdff $(HFD5_LIB) -lhdf5_hl -lhdf5 -lm -lz
...
```
## Turn on electro-coalescence scheme and Set charging rate(default is zero)
```ruby
vi scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm/run.conf
...
- sdm_elecol	 = 0, ! Flag for electro coalescence scheme, 0: No Charge(default), 1: Columb force, 2: Image force, 3: Khain04, 4: Conducting sphere
- sdm_elerate  = 0.0d0, ! Factor for droplets charging rate, used when sdm_elecol!=0 
+ sdm_elecol	 = 4, ! Flag for electro coalescence scheme, 0: No Charge(default), 1: Columb force, 2: Image force, 3: Khain04, 4: Conducting sphere
+ sdm_elerate  = 0.3d0, ! Factor for droplets charging rate (0-7), used when sdm_elecol!=0 
...
```
> [!NOTE]
> Use sdm_elecol to select electro coalescence scheme, 0: No Charge(default), 1: Columb force, 2: Image force, 3: Khain04, 4: Conducting sphere, note that if sdm_elecol not equal 0, must set sdm_elerate to make electro-coalescence work.

## Clean 
```ruby
cd scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm
make allclean ENABLE_SDM=T
sh clean.sh
```

## Compile
```ruby
module purge
module load intel/19.1.3 mpt hdf5/1.12.0 netcdf/4.7.4
cd scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm
make ENABLE_SDM=T
```

## Run the bacth job
```ruby
cd scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm
```
PBS system job script example
```ruby
vi iwaya_run.sh
```

```ruby
#!/bin/bash
#PBS -q S
#PBS -l select=4:ncpus=20:mpiprocs=20
#PBS -l walltime=4:00:00
#PBS -N scale-sdm

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}
module purge
module load intel/19.1.3 mpt hdf5/1.12.0 netcdf/4.7.4

# run
mpiexec_mpt dplace -s1 ./scale-les_init init.conf || exit
mpiexec_mpt dplace -s1 ./scale-les  run.conf || exit
```
check the simulation time step by
```ruby
tail -f LOG.pe000000
```
# Plot panel
## Plot collision-coalescence kernel
Trun on the collision-coalescence kernel output option
```ruby
vi contrib/SDM/sdm_coalescence.f90
```
```
...
- logical,parameter :: output_kernel = .false. 
+ logical,parameter :: output_kernel = .True. 
...
```
And repeat 
>simulation setup

```ruby
cd scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm/cc_output
python plot_cc_kernel.py
```
## Plot the Spatial Structure Using R
```ruby
cd QHYD_2Dplot
module purge
module load intel/19.1.3 mpt hdf5/1.12.0 netcdf/4.7.4 R/4.1.0
Rscript QHYD_2Dplot.R
for i in *.pdf ; do convert -density 400 $i ${i/pdf/png} ; done
for i in QHYD_overlay ; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
cd ..
```
## Plot the Time Series of Water Path and Accumulated Precipitation
```ruby
cd time_series_0D
Rscript time_series_0D.R
for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
for i in *.png; do display $i; done
cd ..
```
## Droplet Size Distribution (Mass Density).
```ruby
cd ptl_dist_1D_R
Rscript ptl_dist_1D.netcdf.R
for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
for i in *.png; do display $i; done
for i in mass_dens_drop; do convert  -delay 20 -loop 0 ${i}.*.png ${i}.gif; done
cd ..
```
