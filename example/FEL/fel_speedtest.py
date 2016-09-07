import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from mpi4py.MPI import COMM_WORLD as comm
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from chimera_tools import *

out_folder = './../fel_test/'
dt = 0.015
Steps2Do = 500
BoxGrid     = (-0.01 ,0.01 ,0.08  ,0.0002,0.001  )
BoxGridBeam = (-0.006,0.006,0.0036,0.0001,0.00012)
g0 = 391.397366
K0 = 1.95

k_res = 2*g0**2/(1+K0**2/2)
gg = g0/(1.+K0**2/2)
vb = (1.-gg**-2)**0.5

solver_in = {'MaxAzimuthMode':2, 'Grid':BoxGrid, 'KxShift':k_res,'TimeStep':dt, 'Rcut':0.05, 'Features':()}
seed_in = {'a0':0.107,'k0':k_res,'x0':0.0,'x_foc':25.0,'Lx':0.00107,'LR':0.00656215} 
Devices = (['MagneticUndulator1D',1.95, np.array([1.,0.,1e9,1.]),np.array([0.,-0.01,0.01])],)
specie_in = {'RandCell':30, 'Charge':-1., 'Density':56.6, 'Mass':1., 'Grid':BoxGridBeam, 'TimeStep':dt, 
'MomentaMeans':(g0,0.,0.), 'MomentaSpreads':(1e-4*g0,0.008,0.008), 'NoiseRemv':(k_res,), 
'Features':(),'Devices':Devices}
MovingFrame = {'TimeStep':dt,'Steps':1,'Velocity':vb}

solver = Solver(comm,solver_in)
solver.add_gauss_beam(seed_in)
specie = Specie(comm,specie_in)
specie.add_particles(specie.gen_parts((-0.00535,0.00535,0.0,0.0034)))
specie.correct_fel()
chimera_in = {'Solvers':(solver,),'Particles':(specie,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(comm,chimera_in)
Chimera.make_halfstep()

if comm.rank==0: timestart = time.time()
for i in xrange(Steps2Do): Chimera.make_step(i*dt)
if comm.rank==0: print 'done in %f minutes' % ((time.time()-timestart)/60.,)
