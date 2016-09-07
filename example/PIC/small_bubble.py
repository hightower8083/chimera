import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from diagnostics import Diagnostics

dt = 0.04
Steps2Do = 601
BoxGrid = (-50.,0.0,20.,0.04,0.25)
out_folder = '../small_bubble1/'

densprof = lambda x: np.interp(x, [0.,40.,50.,300],[0.,1.,0.5,0.5])

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('SpaceCharge','StillAsBackground','Xchunked')}
seed_in = {'a0':3.0,'k0':1.0,'x0':-16.0,'x_foc':45.0,'Lx':4,'LR':4}

specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.005, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(4,41),'Features':('NoSorting',)}
specie2_in = {'FixedCell':(2,2,4), 'Charge':1., 'Density':0.005, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(4,41), 'Features':('Still','NoSorting')}

MovingFrame = {'TimeStep':dt,'Steps':40,'AbsorbLayer':10.,'AddPlasma':densprof,'Features':()}

diags_in = (
{'Type':'Fields','Step':200},
{'Type':'Density','Step':200,'ModeMax':3},
)

solver = Solver(solver_in)
solver.add_gauss_beam(seed_in)		# check if need seed 
specie1 = Specie(specie1_in)
specie2 = Specie(specie2_in)	# check used number of species
chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(chimera_in)

if sys.argv[-1]=='sim':
	Diags = Diagnostics(Chimera,diags_in,'../small_bubble1/')
	Chimera.make_halfstep()
	timestart = time.time()
	for i in xrange(Steps2Do):
		Chimera.make_step(i)
		Diags.do_diags(i)
		istr = str(i)
		while len(istr)<7: istr='0'+istr
		print istr
	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
