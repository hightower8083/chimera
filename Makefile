PGM = fimera

SRC_common = ./f90/fb_io.f90              ./f90/fb_math.f90         \
             ./f90/fb_math_env.f90        ./f90/grid_deps.f90       \
             ./f90/grid_deps_env.f90      ./f90/grid_deps_chnk.f90  \
             ./f90/grid_deps_env_chnk.f90 ./f90/maxwell_solvers.f90 \
             ./f90/particle_tools.f90     ./f90/devices.f90         \
             ./f90/utils.f90              ./f90/synchrad.f90

#FFTI = /home/sources/magins/andriyash/CODES/fftw/include
#FFTL = /home/sources/magins/andriyash/CODES/fftw/lib

#FFTI = /home/sources/magins/andriyash/CODES/fftw-gcc/include
#FFTL = /home/sources/magins/andriyash/CODES/fftw-gcc/lib

#FFTI = /Users/igor/CODES/fft/include
#FFTL = /Users/igor/CODES/fft/lib

FFTI = /usr/local/include
FFTL = /usr/local/lib

#FLAGS_G = -c -DF2PY_REPORT_ON_ARRAY_COPY=1 --opt='-Og -Wall -Wline-truncation  -Wcharacter-truncation \
# -Wextra -Wsurprising  -Waliasing -Wimplicit-interface  -Wunused-parameter  -fwhole-file -fcheck=all  \
# -std=f2008 -pedantic -fbacktrace -fopenmp -lm -lfftw3 -I$(FFTI)' -L$(FFTL) -lm -lfftw3 -lgomp
#FLAGS_I = -c --fcompiler=intelem --opt='-O3 -openmp -xHost -ipo -heap-arrays 24576 -I$(FFTI) -lfftw3' -L$(FFTL) -lm -lfftw3 -liomp5

FLAGS_G = -c --opt='-O3 -ffast-math -march=native -fopenmp -lm -lfftw3 -I$(FFTI)' -L$(FFTL) -lm -lfftw3 -lgomp
FLAGS_I = -c --fcompiler=intelem --opt='-O3 -openmp -xHost -ipo -I$(FFTI) -lfftw3' -L$(FFTL) -lm -lfftw3 -liomp5

F90 = f2py

# This is done when typing 'make gfortran' or simply 'make' (default behavior)
gfortran :
	$(F90) $(FLAGS_G) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

# This is done when typing 'make ifort'
ifort :
	$(F90) $(FLAGS_I) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

clean :
	rm -rf  ./moduls/*.so ./moduls/*.pyc ./*.pyc ./moduls/.nfs0* \
           ./.nfs0* ./f90/.nfs0*  ./doc/.ipynb_checkpoints/* ./*.so.*
