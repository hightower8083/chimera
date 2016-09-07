import numpy as np
from noinput import *

GridFEL = (0.05, 0.02,0.006,0.0003)
k_res = 105599.78482459845

Solvers_in = ( [ 'FEL.1',GridFEL,(0.,1.e10),(0.985*k_res,1.003*k_res),], )
Particles_in = ( [(2e3,56.59688,-1.,1.),('FEL',GridFEL,(1, k_res)), (-0.005,0.005,0.0034), (391.388155,0.,0.), (0.03914,0.00105,0.00105)], )

Devices = ( ['MagneticUndulator1D',1.95, (1.,0.,1e10,4.),0.], )
Seeds_in = ( ['Gauss',0.1,k_res, (0.,0.,0.,29.),(0.00214,0.00656)], )

out_folder = './../tests.out/test_fel_01/'
simul_time = 600.
dt = 0.03
diag_nrg_step = int(2/dt)
diag_fld_step = int(50/dt)
diag_phs_step = int(50/dt)
