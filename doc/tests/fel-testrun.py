from __future__ import print_function,division
from copy import deepcopy
import numpy as np
import sys, time

from chimera.moduls.species import Specie
from chimera.moduls.solvers import Solver
from chimera.moduls.chimera_main import ChimeraRun
from chimera.moduls.diagnostics import Diagnostics
import chimera.moduls.fimera as chimera

K0 = 1.95
Periods = 10
lam0 = 2.8

g0 = 200/0.511
lbx  = 80e-4/lam0
lbr  = 80e-4/lam0
dens = 20e-12/1.6022e-19/(np.pi*80e-4**3)/(1.1e21/2.8e4**2)

gg = g0/(1.+K0**2/2)**.5
k_res = 2*gg**2
vb = (1.-gg**-2)**0.5

Lgx  = 200e-4  /lam0
Rg   = 1000e-4  /lam0
Rg_cut = 700e-4/lam0
Nx,Nr = 120,120

dt = 1./30
SimulationSteps = int(Periods/dt)

solver_in = {
    'Grid':(-0.5*Lgx, 0.5*Lgx, Rg, Lgx/Nx, Rg/Nr), 
    'TimeStep':dt, 'MaxAzimuthMode':0, 
    'KxShift':k_res,'Rcut':Rg_cut ,
    'CoPropagative':vb, 'Xchunked':(4,6),
    'Features':{'NoPoissonCorrection':True,},
}

seed_in = {
    'a0':0.15,'k0':k_res,'x0':-30e-4/lam0,
    'x_foc':70.0/lam0,'Lx':15e-4/lam0,'LR':180e-4/lam0
}

solver = Solver(solver_in)
solver.add_gauss_beam(seed_in)

beam_in = {
    'Grid':(-0.5*Lgx, 0.5*Lgx, Rg, Lgx/Nx, lbr/Nr),'TimeStep':dt,'Charge':-1.,'Mass':1.,
    'Density':dens, 'RandCell':50,'MomentaMeans':(g0,0.,0.), 'MomentaSpreads':(1e-4*g0,2e-5*g0,2e-5*g0),
    'Xchunked':(4,6),
}

beam_in['Devices'] = ([chimera.undul_analytic,np.array([K0, 1., 1., Periods])],)

beam = Specie(beam_in)
beam.add_particles(*beam.gen_parts(Domain=(-0.5*lbx, 0.5*lbx, 0.0, lbr)))
beam.denoise((k_res,))

MovingFrame = {
    'TimeStep':dt,'Steps':1,'Velocity':vb,'Features':('Staged','NoSorting')
}

chimera_in = {'Solvers':(solver,),'Particles':(beam,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(chimera_in)

Diags = Diagnostics(Chimera,(),out_folder=None)

diag_step = int(1./dt)
nrg = Diags.nrg_out({'Features':('Return',)})[0]
pwr = Diags.pwr_out({'Features':('Return',)})[0]

diag_nrg = np.zeros((int(SimulationSteps/diag_step)+1,nrg.shape[0]))
diag_pwr = np.zeros((int(SimulationSteps/diag_step)+1,pwr.shape[0]))

diag_nrg[0] = nrg
diag_pwr[0] = pwr

ti = time.time()
for i in range(1,SimulationSteps+1):
    Chimera.make_step(i)
    sys.stdout.write('\r'+str(i)+'/'+str(SimulationSteps))
    sys.stdout.flush()
    if np.mod(i,diag_step)==0:
        diag_nrg[int(i/diag_step)] = Diags.nrg_out({'Features':('Return',)})[0]
        diag_pwr[int(i/diag_step)] = Diags.pwr_out({'Features':('Return',)})[0]
print('\nDone in {:g} seconds'.format(time.time()-ti))
