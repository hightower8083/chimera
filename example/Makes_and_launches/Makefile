PGM = fimera

SRC_common = ./src/fb_math.f90  ./src/grid_deps.f90  \
             ./src/devices.f90  ./src/particle_tools.f90 \
             ./src/utils.f90    ./src/synchrad.f90

#FFTI = /home/sources/magins/andriyash/CODES/fftw/include
#FFTL = /home/sources/magins/andriyash/CODES/fftw/lib

#FFTI = /home/sources/magins/andriyash/CODES/fftw-gcc/include
#FFTL = /home/sources/magins/andriyash/CODES/fftw-gcc/lib

#FFTI = /Users/igor/CODES/fft/include
#FFTL = /Users/igor/CODES/fft/lib

FFTI = /usr/local/include
FFTL = /usr/local/lib

#FLAGS_G = -c -DF2PY_REPORT_ON_ARRAY_COPY=1 --opt='-O3 -ffast-math -lm -lfftw3 -I/usr/local/include/' -lm -lfftw3
#FLAGS_I = -c --fcompiler=intelem --opt='-O3 -openmp -xHost -ipo -heap-arrays 24576 -I$(FFTI) -lfftw3' -L$(FFTL) -lm -lfftw3 -liomp5

FLAGS_G = -c --opt='-O3 -ffast-math -fopenmp -lm -lfftw3 -I$(FFTI)' -L$(FFTL) -lm -lfftw3 -lgomp
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
	rm -f  ./moduls/*.so ./moduls/*.pyc ./src/*.pyc \
          ./*.pyc ./moduls/.nfs0* ./.nfs0* ./src/.nfs0* \
          ./.ipynb_checkpoints/*
