# CHIMERA: a code for FEL and laser plasma simulations

[![Build Status master](https://img.shields.io/travis/hightower8083/chimera/master.svg?label=master)](https://travis-ci.org/hightower8083/chimera/branches)
[![Build Status dev](https://img.shields.io/travis/hightower8083/chimera/dev.svg?label=dev)](https://travis-ci.org/hightower8083/chimera/branches)

by Igor A Andriyash (<igor.andriyash@gmail.com>)

CHIMERA is a relativistic electromagnetic particle-in-cell code, based on a quasi-cylindric pseudo-spectral analytical time domain (PSATD) Maxwell solver. More details on 
this method can be found in the original publications [<cite>[1]</cite>,<cite>[2]</cite>,<cite>[3]</cite>]. 

System requirements
- code runs under Linux or MacOS
- Fortran 90/95 compiler with OpenMP support (e.g. recent gfortran or ifort)
- Python with NumPy, SciPy; Additionally Matplotlib, Ipython and Jypyter are recommended
- FFTW3 (http://www.fftw.org), better to be compiled with the same compiler and "-fPIC" option enabled

To install CHIMERA
- clone the code folder into your working directory and add it to the PYTHONPATH
- Check that Makefile contains a correct address to FFTW3 /lib/ and /bin/
- compile the Fortran modules using *make*. 

To run CHIMERA in multiprocessor mode specify the OMP_NUM_THREADS variable. For more information see demo in ./doc/

\[[1]\] Igor A. Andriyash, Rémi Lehe and Agustin Lifschitz, *Laser-plasma interactions with a Fourier-Bessel particle-in-cell method*, Physics of Plasmas **23**, 033110 
(2016)

\[[2]\] Rémi Lehe, Manuel Kirchen, Igor A. Andriyash, Brendan B. Godfrey and Jean-Luc Vay, *A spectral, quasi-cylindrical and dispersion-free Particle-In-Cell algorithm*, 
Computer Physics Communications **203**, 66 (2016)

\[[3]\] Igor A. Andriyash, Rémi Lehe and Victor Malka, *A spectral unaveraged algorithm for free electron laser simulations*, Journal of Computational Physics **282**, 397 (2015)

[1]:http://dx.doi.org/10.1063/1.4943281
[2]:http://dx.doi.org/10.1016/j.cpc.2016.02.007
[3]:http://dx.doi.org/10.1016/j.jcp.2014.11.026

