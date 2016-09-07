import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from mpi4py.MPI import COMM_WORLD as comm
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from chimera_tools import *

def stack_overlapped(x):
	loc_arrays = []
	for i in xrange(len(x)-1):
		dat = x[i].copy()
		dat[:,-1] += x[i+1][:,0]
		dat = dat[:,1:]
		loc_arrays.append(dat)
	return np.concatenate(loc_arrays,axis=1)

out_folder = '/work/sources/magins/andriyash/code/tests.out/lpa_andreas_01/'
dt = 0.04
Steps2Do = 11
BoxGrid = (-140.,0.0,60.,0.04,0.3)

densprof = np.vectorize( lambda x,y,z:1.)

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('Rsliced','SpaceCharge','StillAsBackground')}

seed_in = {'a0':2.,'k0':1.0,'x0':-42.0,'x_foc':1250.,'Lx':11.25,'LR':12.5}

specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.00465, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt, 
'Features':('Rsliced',)}

specie2_in = {'FixedCell':(2,2,4), 'Charge': 1., 'Density':0.00465, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt, 
'Features':('Rsliced', 'Still')}

solver = Solver(comm,solver_in)
solver.add_gauss_beam(seed_in)      
specie1 = Specie(comm,specie1_in)
specie1.add_particles( specie1.gen_parts((BoxGrid[0]+10.,BoxGrid[1],0.0,BoxGrid[2]),) )
specie2 = Specie(comm,specie2_in)
specie2.add_particles( specie2.gen_parts((BoxGrid[0]+10.,BoxGrid[1],0.0,BoxGrid[2]),) )
chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':()}
Chimera = ChimeraRun(comm,chimera_in)

Chimera.make_halfstep()

if comm.rank==0: timestart = time.time()
for i in xrange(Steps2Do):
   Chimera.make_step(i*dt)
   if comm.rank==0 and np.mod(i, (Steps2Do-1)//10)==0: print i
if comm.rank==0: print 'done in %f minutes' % ((time.time()-timestart)/60.,)
