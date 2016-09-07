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
Steps2Do = 40
BoxGrid = (-10.,10.0,10.,0.05,0.05)

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('SpaceCharge',)}

specie_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.01, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(2,1)}
#specie_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.01, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt}

solver = Solver(solver_in)
specie = Specie(specie_in)

specie.add_particles(specie.gen_parts(Domain=(-3.,3.,0.,3.),) )

chimera_in = {'Solvers':(solver,),'Particles':(specie,)}
Chimera = ChimeraRun(chimera_in)
Chimera.make_halfstep()

t=time.time()
for i in range(200):
	Chimera.make_step(i)
print time.time()-t
