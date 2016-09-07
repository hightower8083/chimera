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

fld_out_step = 100
dns_out_step = 100
phs_out_step = 1000

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
	if comm.rank==0: print 'done in %f minutes' % ((time.time()-timestart)/60.,)
