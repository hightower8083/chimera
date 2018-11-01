# CHIMERA: a code for FEL and laser plasma simulations

[![Build Status master](https://img.shields.io/travis/hightower8083/chimera/master.svg?label=master)](https://travis-ci.org/hightower8083/chimera/branches)
[![Build Status dev](https://img.shields.io/travis/hightower8083/chimera/dev.svg?label=dev)](https://travis-ci.org/hightower8083/chimera/branches)
[![CRAN](https://img.shields.io/cran/l/devtools.svg)](LICENSE)

by Igor A Andriyash (<igor.andriyash@gmail.com>)

CHIMERA is a relativistic electromagnetic particle-in-cell code, based on a quasi-cylindric pseudo-spectral analytical time domain (PSATD) Maxwell solver. More details on this method can be found in the original publications [<cite>[1]</cite>,<cite>[2]</cite>,<cite>[3]</cite>]. 

## Installation

The code was tested only under Linux and MacOS. The recommended installation is through the [Anaconda](https://www.continuum.io/why-anaconda) distribution.
If Anaconda is not your default Python installation, download and install it from [here](https://www.continuum.io/downloads).

- install the dependances:
```
conda install gcc scipy h5py
conda install -c conda-forge fftw
```
- clone the code source and `cd` into the folder (e.g. into `~/src/`)
```
git clone https://github.com/hightower8083/chimera.git
cd chimera
```
- compile the Fortran 90 subroutines and install the code:
```
make
make install
```
- code may be uninstalled using:
```
make uninstall
make clean
```

**NB**: 
- you may need to adapt Makefile for your installation of Anaconda and FFTW. 
- you may need to heck the correct address to FFTW3 /lib/ and /bin/


## Running

To run CHIMERA in multiprocessor mode specify the OMP_NUM_THREADS variable. 

For more information see demos in ./doc/


## References

\[[1]\] Igor A. Andriyash, Rémi Lehe and Agustin Lifschitz, *Laser-plasma interactions with a Fourier-Bessel particle-in-cell method*, Physics of Plasmas **23**, 033110 
(2016)

\[[2]\] Rémi Lehe, Manuel Kirchen, Igor A. Andriyash, Brendan B. Godfrey and Jean-Luc Vay, *A spectral, quasi-cylindrical and dispersion-free Particle-In-Cell algorithm*, 
Computer Physics Communications **203**, 66 (2016)

\[[3]\] Igor A. Andriyash, Rémi Lehe and Victor Malka, *A spectral unaveraged algorithm for free electron laser simulations*, Journal of Computational Physics **282**, 397 (2015)

[1]:http://dx.doi.org/10.1063/1.4943281
[2]:http://dx.doi.org/10.1016/j.cpc.2016.02.007
[3]:http://dx.doi.org/10.1016/j.jcp.2014.11.026

