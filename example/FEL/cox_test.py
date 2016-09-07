import sys, time
sys.path.append('./src/')
sys.path.append('./moduls/')

import fimera as chimera
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from diagnostics import Diagnostics

import numpy as np

K0 = 1.71856
g0 = 352.2505
Periods = 102
StepsPerPeriod = 50

k_res = 2*g0**2/(1+K0**2/2)
gg = g0/(1.+K0**2/2)**.5
vb = (1.-gg**-2)**0.5

out_folder = './../cox_01/'
dt = 1./StepsPerPeriod
Steps2Do = int(Periods/dt)+1

densprof = lambda x,y,z: np.exp(-0.5*x**2/5e-5**2)*np.exp(-0.5*(y**2+z**2)/1.42e-4**2)

BoxGrid     = (-0.005 ,0.005,0.06,0.000025, 0.0003  )
BoxGridBeam = (-0.005 ,0.005,0.01,0.000025, 0.000005)

solver_in = {'MaxAzimuthMode':0, 'Grid':BoxGrid, 'KxShift':  k_res,'TimeStep':dt, 'Rcut':0.03, 'Features':('Xchunked',)}

seed_in = {'a0':0.074,'k0':k_res,'x0':0.0,'x_foc':10.0,'Lx':0.0016,'LR':0.00892}
Devices = ({'DEV':[chimera.undul_analytic_taper,np.array([K0, 1.,1.,100.,0.027])]},)

specie_in = {'RandCell':50, 'Charge':-1., 'Density':6e5, 'Mass':1., 'Grid':BoxGridBeam, 'TimeStep':dt,
'MomentaMeans':(g0,0.,0.), 'MomentaSpreads':(3.5225,0.35225,0.35225),'Devices':Devices,'Xchunked':(4,6)}

diags_in = (
{'Type':'EnergyEM' ,'Step':int(1./dt)},
{'Type':'Power'    ,'Step':int(1./dt)},
{'Type':'Particles','Step':int(20./dt)},
)

solver = Solver(solver_in)
specie = Specie(specie_in)
MovingFrame = {'TimeStep':dt,'Steps':1,'Velocity':vb,'Features':('Staged','NoSorting')}
chimera_in = {'Solvers':(solver,),'Particles':(specie,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(chimera_in)

if sys.argv[-1]=='sim':
	Diags = Diagnostics(Chimera,diags_in,out_folder)
	solver.add_gauss_beam(seed_in)
	specie.add_particles(specie.gen_parts(Domain=(-2e-4,2e-4,0.0,0.0005),ProfileFunc=densprof))
	specie.coxinel_line()
	specie.denoise((k_res,2*k_res,3*k_res))

	Chimera.make_halfstep()

	timestart = time.time()
	for i in xrange(Steps2Do):
		Diags.do_diags(i)
		Chimera.make_step(i)
		if np.mod(i,Steps2Do/20)==0: print i
	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
