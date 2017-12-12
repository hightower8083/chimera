from __future__ import print_function,division
from copy import deepcopy
import numpy as np
import sys, time

from chimera.moduls.species import Specie
from chimera.moduls.solvers import Solver
from chimera.moduls.chimera_main import ChimeraRun
from chimera.moduls.diagnostics import Diagnostics
import chimera.moduls.fimera as chimera

xgmin, xgmax = -20., 1.
Rg = 16.0

dx = 0.04
dr = 0.25
dt = dx

SimulationSteps = int(22/dt)

a0 = 3.
Lx = 4.
w0 = 4.
x_foc = 45.
x0 = -14.

dens = 0.005
densprof = lambda x: np.interp(x, [0.,40.,50.,300],[0.,1.,0.5,0.5])

solver_in = {
    'Grid':(xgmin, xgmax,Rg,dx,dr),'TimeStep':dt, 'MaxAzimuthMode':1,
    'Xchunked':(4,10),'Features':('SpaceCharge','StillAsBackground',)
}

laser_in = {
    'a0':a0,'k0':1.0,'x0':x0,'x_foc':x_foc,'Lx':Lx,'LR':w0
}

solver = Solver(solver_in)
solver.add_gauss_beam(laser_in)

electrons_in = {
    'Grid':(xgmin, xgmax,Rg,dx,dr), 'TimeStep':dt,
    'Density':dens, 'FixedCell':(2,2,4),
    'Xchunked':(4,10),'Features':('NoSorting',)
}

ions_in = deepcopy(electrons_in)
ions_in['Charge'] = 1
ions_in['Mass'] = 1886
ions_in['Features'] += ('Still',)

electrons = Specie(electrons_in)
ions = Specie(ions_in)

MovingFrame = {
    'TimeStep':dt,'Steps':10,'AbsorbLayer':175,
    'AddPlasma':densprof,'Features':('IonsOnTop',)}

chimera_in = {
    'Solvers':(solver, ), 'Particles':(electrons, ions),
    'MovingFrames':(MovingFrame, )}

Chimera = ChimeraRun(chimera_in)
Diags = Diagnostics(Chimera, (), out_folder=None)

ti = time.time()
for i in range(1,SimulationSteps+1):
    Chimera.make_step(i)
    sys.stdout.write('\r'+str(i)+'/'+str(SimulationSteps))
    sys.stdout.flush()
print('\nDone in {:g} seconds'.format(time.time()-ti) )
