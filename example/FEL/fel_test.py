import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from diagnostics import Diagnostics

out_folder = './../fel_test/'

g0 = 391.397366
K0 = 1.95
dt = 0.015
Steps2Do = int(450./dt)+1

BoxGrid     = (-0.01,0.01,0.08,0.0001,0.001 )
BoxGridBeam = (-0.01,0.01,0.04,0.0001,0.0001)

k_res = 2*g0**2/(1+K0**2/2)
gg = g0/(1.+K0**2/2)**0.5
vb = (1.-gg**-2)**0.5

solver_in = {'MaxAzimuthMode':0, 'Grid':BoxGrid, 'KxShift':  k_res,'TimeStep':dt, 'Rcut':0.04, 'Features':('Xchunked',)}
#solver_in = {'MaxAzimuthMode':2, 'Grid':BoxGrid, 'KxShift':  k_res,'TimeStep':dt, 'Rcut':0.04, 'Features':('Xchunked',)}
#solver_in2 ={'MaxAzimuthMode':2, 'Grid':BoxGrid, 'KxShift':2*k_res,'TimeStep':dt, 'Rcut':0.04, 'Features':('Xchunked',)}
#solver_in3 ={'MaxAzimuthMode':2, 'Grid':BoxGrid, 'KxShift':3*k_res,'TimeStep':dt, 'Rcut':0.04, 'Features':('Xchunked',)}

seed_in = {'a0':0.107,'k0':k_res,'x0':0.0,'x_foc':25.0,'Lx':0.00107,'LR':0.00656215}
devices_in = (
{'DEV':['MagneticUndulator1D',1.95, np.array([1.,0.,1e9,1.]),np.array([0.,-0.01,0.01])]},)

specie_in = {'RandCell':60, 'Charge':-1., 'Density':56.6, 'Mass':1., 'Grid':BoxGridBeam, 'TimeStep':dt,
'MomentaMeans':(g0,0.,0.), 'MomentaSpreads':(1e-4*g0,0.008,0.008),'Devices':devices_in,'Xchunked':(4,6)}

MovingFrame = {'TimeStep':dt,'Steps':1,'Velocity':vb,'Features':('Staged','NoSorting')}

diags_in = (
{'Type':'EnergyEM','Step':int(1./dt)},
{'Type':'Power'   ,'Step':int(1./dt)},
{'Type':'Density' ,'Step':int(10./dt),'Features':{'MaxMode':2,},},
)

solver = Solver(solver_in)
solver.add_gauss_beam(seed_in)
#solver3 = Solver(solver_in3)

specie = Specie(specie_in)
specie.add_particles(specie.gen_parts(Domain=(-0.00535,0.00535,0.0,0.0034)))
specie.denoise((k_res,3*k_res,))
specie.correct_fel()

chimera_in = {'Solvers':(solver,),'Particles':(specie,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(chimera_in)

if sys.argv[-1]=='sim':
	Chimera.make_halfstep()
	Diags = Diagnostics(Chimera,diags_in,out_folder)
	timestart = time.time()
	for i in xrange(Steps2Do):
		Chimera.make_step(i)
		Diags.do_diags(i)
		if np.mod(i,Steps2Do/20)==0: print i
	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
