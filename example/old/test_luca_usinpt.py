import numpy as np
from noinput import *

GridFEL = (0.04,0.08,0.05,0.001)
GridFEL0= (0.0,0.0,0.0036,0.00007)
k_res = 105600.205166

Solvers_in = (
['FEL.0',GridFEL,(0.,1.e10),(0.989*k_res,1.007*k_res),],
['FEL.0',GridFEL,(0.,1.e10),(2.985*k_res,3.009*k_res),],
)

Devices = [
['MagneticUndulator1D',1.95, (1.,0.,1e9,1.),[0.,-0.01,0.01]],
]

Seeds_in = (
['Gauss',0.107,k_res, (0.,0.,0.,25.),(0.00107,0.00656215)],
)

Particles_in = (
[(2400,56.6,-1.,1.),('FEL',GridFEL0,(k_res,2*k_res,3*k_res)),
(-0.00535,0.00535,0.0034),(391.397366,0.,0.),(0.03914,0.008,0.008)],
)

#PlasmaProfiles_in = (
#('gauss', 'r', (0,0.0034)),
#)

dopoiss = 0
out_folder = '../fel_luca_flattop/'
simul_time = 400.
dt = 0.01
#track_diag = True
#diag_int =   10.
#track_step = 2#int(1./dt)
diag_nrg_step = int(1/dt)
diag_pwr_step = int(200/dt)
diag_fld_step = int(100/dt)
diag_phs_step = int(100/dt)
