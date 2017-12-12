F90 = f2py
PGM = fimera
FFTPATH = $(HOME)/CODES/fftw

SRC_common = ./f90/fb_io.f90              ./f90/fb_math.f90         \
             ./f90/fb_math_env.f90        ./f90/grid_deps.f90       \
             ./f90/grid_deps_env.f90      ./f90/grid_deps_chnk.f90  \
             ./f90/grid_deps_env_chnk.f90 ./f90/maxwell_solvers.f90 \
             ./f90/particle_tools.f90     ./f90/devices.f90         \
             ./f90/utils.f90              ./f90/SR.f90

FLAGS_G = -c --opt='-O3 -ffast-math -march=native -fopenmp          \
          -lm -lfftw3 -I$(FFTPATH)/include' -L$(FFTPATH)/lib -lm -lfftw3 -lgomp

FLAGS_I = -c --fcompiler=intelem --opt='-O3 -openmp -xHost          \
          -ipo -I$(FFTPATH)/include -lfftw3' -L$(FFTPATH)/lib -lm -lfftw3 -liomp5

FLAGS_GD = -c -DF2PY_REPORT_ON_ARRAY_COPY=1 --opt='-Og -Wall        \
           -Wline-truncation  -Wcharacter-truncation -Wextra        \
           -Wsurprising  -Waliasing -Wimplicit-interface            \
           -Wunused-parameter  -fwhole-file -fcheck=all             \
           -std=f2008 -pedantic -fbacktrace -fopenmp -lm -lfftw3    \
           -I$(FFTPATH)/include' -L$(FFTPATH)/lib -lm -lfftw3 -lgomp

FLAGS_ID = -c --fcompiler=intelem --opt='-O3 -openmp -xHost -ipo    \
           -heap-arrays 24576 -I$(FFTPATH)/include -lfftw3' -L$(FFTPATH)/lib -lm   \
           -lfftw3 -liomp5

# This is done when typing 'make gfortran' or simply 'make' (default behavior)
gfortran :
	$(F90) $(FLAGS_G) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

# This is done when typing 'make ifort'
ifort :
	$(F90) $(FLAGS_I) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

# This is done when typing 'make debug' or 'make idebug'
debug :
	$(F90) $(FLAGS_GD) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

idebug :
	$(F90) $(FLAGS_ID) -m $(PGM) $(SRC_common)
	mv *.so ./moduls/

clean :
	rm -rf  ./*.so* ./moduls/*.so*                \
           ./*.pyc ./moduls/*.pyc ./utils/*.pyc  \
           ./.nfs0* ./moduls/.nfs0* ./f90/.nfs0* \
           ./doc/.ipynb_checkpoints ./.DS_Store
