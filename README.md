Introduction of electro-coalescence treatment of droplets SCALE-SDM 0.0-2.2.2 for reproduce the results of Zhang et al. (2023, GMD)

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
+NETCDF_INCLUDE ?= $(NETCDF_INC)
+NETCDF_LIBS    ?= $(NETCDF_LIB) $(HDF5_LIB) -lm -lz
...
```
## Turn on the conducting spheres treatment option
```ruby
vi contrib/SDM/sdm_coalescence.f90
...
- real(RP) :: alpha=0.0D0 ! coefficient to decide the charge amount (Andronache 2004)
+ real(RP) :: alpha=0.2D0 ! coefficient to decide the charge amount (Andronache 2004)
...
```
```ruby
...
- !if ( (sd_r1 > (sd_r2*1.0d2)) .or. (sd_r2 > (sd_r1*1.0d2)) ) then
- !    call electro_coalescence_efficiency_image_charge(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
- !else
- !    call electro_coalescence_efficiency_conducting_sphere(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
- !end if
- !if ((sd_r1< 1.0e-7) .or. (sd_r2<1.0e-7)) then
- !    eff_elc = 0.0
- !end if
-  eff_elc = 0.0
+ if ( (sd_r1 > (sd_r2*1.0d2)) .or. (sd_r2 > (sd_r1*1.0d2)) ) then
+     call electro_coalescence_efficiency_image_charge(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
+ else
+     call electro_coalescence_efficiency_conducting_sphere(eff_elc,sd_r1,sd_r2,sd_vz1,sd_vz2,lmd_crs,vis_crs)
+ end if
+ if ((sd_r1< 1.0e-7) .or. (sd_r2<1.0e-7)) then
+     eff_elc = 0.0
+ end if
+ ! eff_elc = 0.0
...
```

## Clean 
```ruby
cd scale-les/test/case/warmbubble/2D_Lasher-trapp05_mod_sdm
make allclean ENABLE_SDM=T
sh clean.sh
```

## Compile
```ruby
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
module load intel/17.0.0 mpt hdf5/1.8.12 netcdf/4.4.1

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
module load intel/2015.6.233 zlib/1.2.8 xz/5.2.4 pcre/8.40 bzip2/1.0.6 openssl/1.1.1a curl/7.63.0 R/3.4.3
Rscript QHYD_2Dplot.R
for i in *.pdf ; do convert -density 400 $i ${i/pdf/png} ; done
for i in QHYD_overlay ; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
cd ..
```
## Plot the Time Series of Water Path and Accumulated Precipitation
```ruby
cd time_series_0D
module purge
module load intel/2015.6.233 zlib/1.2.8 xz/5.2.4 pcre/8.40 bzip2/1.0.6 openssl/1.1.1a curl/7.63.0 R/3.4.3
Rscript time_series_0D.R
for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
for i in *.png; do display $i; done
cd ..
```
## Make a 2D Color Map of Ice Particle Distribution
```ruby
cd ptl_dist_2Dmap_R
module purge
module load intel/2015.6.233 zlib/1.2.8 xz/5.2.4 pcre/8.40 bzip2/1.0.6 openssl/1.1.1a curl/7.63.0 R/3.4.3
Rscript ptl_dist_2Dmap.netcdf.overlay.R
for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
for i in *.png; do display $i; done
animate -delay 20 *.png
for i in aspect density mass term_vel_ice ; do convert  -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
cd ..
```
## Aerosol Size Distribution (Number Density), Freezing Temperature Distribution, Droplet Size Distribution (Mass Density), and Terminal Velocities of Droplets and Ice.
```ruby
cd ptl_dist_1D_R
module purge
module load intel/2015.6.233 zlib/1.2.8 xz/5.2.4 pcre/8.40 bzip2/1.0.6 openssl/1.1.1a curl/7.63.0 R/3.4.3
Rscript ptl_dist_1D.netcdf.R
less numbers.txt
for i in *.pdf; do convert -density 400 $i ${i/pdf/png}; done
for i in *.png; do display $i; done
for i in mass_dens_drop num_dens_amsul prob_dens_freezing_temp term_vel_drop term_vel_ice xz_SD; do convert  -delay 20 -loop 0 ${i}.*.png ${i}.gif; done
cd ..
```
