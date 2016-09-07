#! /usr/bin/env python
import numpy as np
import sys, os
sys.path.append('./moduls/')
sys.path.append('./src/')
from time import time
from mpi4py import MPI
from input_proc import SimpleData
from synchrotron import Calculator

# script args: regime ('coh','incoh','both'), component ('x','y','z','all'), number of particles to take

comm = MPI.COMM_WORLD
T0 = time()
data = SimpleData()
if comm.rank == 0: os.system('cp -f input.py '+data.out_folder)
calc = Calculator(comm,data,sys.argv)
calc.main_calc()
timing = time()-T0
calc.writeout(timing)
