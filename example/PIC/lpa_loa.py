import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from diagnostics import Diagnostics

dt = 0.04
Steps2Do = 140000
BoxGrid = (-90.,0.0,55.,0.04,0.35)
out_folder = '../tests.out/lpa_loa/'

densprof = lambda x: np.interp(x*0.8e-3, [0.,0.8,0.95,4.25],[0.,1.,0.35,0.0])

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('SpaceCharge','StillAsBackground','Xchunked')}
seed_in = {'a0':1.53,'k0':1.0,'x0':-35.,'x_foc':1000.,'Lx':11.25,'LR':17.5}

specie_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.0116, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(4,21),'Features':('NoSorting',)}
specie2_in = {'FixedCell':(2,2,4), 'Charge':1., 'Density':0.0116, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(4,21), 'Features':('Still','NoSorting')}

MovingFrame = {'TimeStep':dt,'Steps':20,'AbsorbLayer':10.,'AddPlasma':densprof,'Features':()}

diags_in = (
{'Type':'Fields','Step':500},
{'Type':'Density','Step':500},
{'Type':'Particles','Step':2000},
)

solver = Solver(solver_in)

if sys.argv[-1]=='sim':
	solver.add_gauss_beam(seed_in)
	specie = Specie(specie_in)
	specie2 = Specie(specie2_in)
	chimera_in = {'Solvers':(solver,),'Particles':(specie,specie2,),'MovingFrames':(MovingFrame,)}
	Chimera = ChimeraRun(chimera_in)
	Diags = Diagnostics(Chimera,diags_in,out_folder)
	Chimera.make_halfstep()
	timestart = time.time()
	for i in xrange(Steps2Do):
		Chimera.make_step(i)
		Diags.do_diags(i)
		ff = open('out',mode='a');	ff.write(str(i)+'\n')

	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
