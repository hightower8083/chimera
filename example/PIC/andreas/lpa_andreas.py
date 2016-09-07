import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from mpi4py.MPI import COMM_WORLD as comm
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from chimera_tools import *

from RunRestart import StartFromSave
import cPickle as pickle
runsave = StartFromSave()

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
Steps2Do = 46876
BoxGrid = (-140.,0.0,63.,0.04,0.3)
fld_out_step = 1875
dns_out_step = 625
phs_out_step = 9375

densprof = np.vectorize( lambda x,y,z: (x>=0.)*(x<=1250.)*x/1250. + (x>1250.)*(x<=2500.)*(1-0.5*(x-1250.)/1250.))
solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('Rsliced','SpaceCharge','StillAsBackground')}
seed_in = {'a0':2.,'k0':1.0,'x0':-42.0,'x_foc':1250.,'Lx':11.25,'LR':12.5}
specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.00465, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt, 
'Features':('Rsliced',)}
specie2_in = {'FixedCell':(2,2,4), 'Charge': 1., 'Density':0.00465, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt, 
'Features':('Rsliced', 'Still')}
MovingFrame = {'TimeStep':dt,'Steps':20,'AbsorbLayer':10.,'AddPlasma':densprof}

if sys.argv[-1]=='sim':
	solver = Solver(comm,solver_in)
	solver.add_gauss_beam(seed_in)		# check if need seed 

	specie1 = Specie(comm,specie1_in)
	specie2 = Specie(comm,specie2_in)	# check used number of species
	chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':(MovingFrame,)}
	Chimera = ChimeraRun(comm,chimera_in)
	Chimera.make_halfstep()

	timestart = time.time()
	if comm.rank==0:
		os.system('rm -rf '+out_folder)
		os.system('mkdir ' +out_folder)
	for i in xrange(Steps2Do):
		Chimera.make_step(i*dt)
		if comm.rank==0: print i
		if np.mod(i,fld_out_step)==0:
			ee = comm.gather(solver.EB[:,1:])
			if comm.rank==0:
				istr = str(i)
				while len(istr)<7: istr='0'+istr
				np.save(out_folder+'ee_'+istr+'.npy',np.concatenate(ee,axis=1))
		if np.mod(i,dns_out_step)==0:
			dens = comm.gather(specie1.get_dens_on_grid(solver_in['MaxAzimuthMode']))
			if comm.rank==0:
				istr = str(i)
				while len(istr)<7: istr='0'+istr
				if comm.size!=1:
					np.save(out_folder+'edens_'+istr+'.npy',stack_overlapped(dens))
				else:
					np.save(out_folder+'edens_'+istr+'.npy',np.concatenate(dens,axis=1))
		if np.mod(i,phs_out_step)==0:
			phs = comm.gather(specie1.particles)
			if comm.rank==0:
				istr = str(i)
				while len(istr)<7: istr='0'+istr
				np.save(out_folder+'phs'+istr+'.npy', np.concatenate(phs,axis=1))

	runsave.SaveRun(Chimera)
	runsave.Step2Start = i
	file2save = open(out_folder+'savedrun_'+str(comm.rank)+'.p','wb')
	pickle.dump(runsave, file2save, pickle.HIGHEST_PROTOCOL)

	if comm.rank==0: print 'done in %f minutes' % ((time.time()-timestart)/60.,)
