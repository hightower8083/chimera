import numpy as np
from noinput import *

GridFEL = (0.005,0.06,0.03,0.0003)
GridFEL0= (0.005,0.0005,0.0005,0.000005)
k_res = 98811.0

Devices = (
['MagneticUndulator1D',1.71856, (1.,0.,1e10,4.),0.],
)

Seeds_in = (
['Gauss',0.076,k_res, (0.0,0.,0.,10.),(0.001,0.00892)],
)

Solvers_in = (
[ 'FEL.0',GridFEL,(0.,1.e10),(0.9*k_res,1.1*k_res),],
)

Particles_in = (
[(1e3,6e5,-1.,1.),('FEL',GridFEL0,(1, k_res)),
(-2e-4,2e-4,0.0005), (352.2505,0.,0.), (3.5225,0.35225,0.35225)],
)

PlasmaProfiles_in = (
('gauss', 'x', (0,    5e-5)),
('gauss', 'r', (0, 1.42e-4)),
)

out_folder = './../tests.out/cox_fel_01/'
simul_time = 150.
dt = 0.04
diag_nrg_step = int(1/dt)
diag_fld_step = int(10/dt)
diag_phs_step = int(10/dt)
