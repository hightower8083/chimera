import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from chimera_tools import *
import fimera as chimera

dt = 0.04
Steps2Do = 120
BoxGrid = (0.,70.0,50.,0.04,0.25)

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('SpaceCharge','StillAsBackground')}
seed_in = {'a0':4.0,'k0':1.0,'x0':-45.0,'x_foc':-45.0,'Lx':12.5,'LR':20.}

#specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.0005734, 'Mass':1.,'MomentaSpreads':(0.01,0.01,0.01), 'Grid':BoxGrid, 'TimeStep':dt}
#specie2_in = {'FixedCell':(2,2,4), 'Charge': 1., 'Density':0.0005734, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt, 'Features':('Still',)}

specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.0005734, 'Mass':1.,'MomentaSpreads':(0.01,0.01,0.01), 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(8,20), 
'Features':('NoSorting',)}
specie2_in = {'FixedCell':(2,2,4), 'Charge': 1., 'Density':0.0005734, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt, 'Features':('Still','NoSorting',),'Xchunked':(8,20)}
MovingFrame = {'TimeStep':dt,'Steps':5,'AbsorbLayer':10.,'AddPlasma':None, 'Features':()}

solver = Solver(solver_in)
specie1 = Specie(specie1_in)
specie1.add_particles( specie1.gen_parts(Domain=(BoxGrid[0]+10.,BoxGrid[1],0.0,BoxGrid[2]),) )
chimera_in = {'Solvers':(solver,),'Particles':(specie1,)}
Chimera = ChimeraRun(chimera_in)

if sys.argv[-1]=='go':
	solver.add_gauss_beam(seed_in)
	specie2 = Specie(specie2_in)
	specie2.add_particles( specie2.gen_parts(Domain=(BoxGrid[0]+10.,BoxGrid[1],0.0,BoxGrid[2]),) )
	chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':()}
	Chimera = ChimeraRun(chimera_in)

	Chimera.make_halfstep()

	timestart = time.time()
	for i in xrange(Steps2Do):
		Chimera.make_step(i)
		if np.mod(i, (Steps2Do-1)//10)==0: print i
	print 'done with %f seconds per step' % ((time.time()-timestart)/Steps2Do,)
