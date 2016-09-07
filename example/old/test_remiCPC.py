import numpy as np
from noinput import *

GridPIC = (-15.,0.,60.,0.04,0.36)

Solvers_in = ( [ 'PIC.1',GridPIC, (0.,1.e10),(1,((300.,1.),(30.,0.)),1)], )
Seeds_in = ( ['Gauss', 4., 1., (-40.,0.,0.,-87.5),(12.5,20.)], )
RunWindows_in = ( (0.,1e9, 1.,10,(10., 0.)), )

PlasmaProfiles_in = (('ramp_up','x',(0, 125.)), ('plateau','x',(125.,1e9 )),('ramp_up','x',(125., 250.)), ('plateau','x',(250., 500.)),('ramp_down','x',(500.,750.)),)
pf=2.
Particles_in = (
[((16,4),5.818e-4/pf,-1.,1.), ('PIC',GridPIC,0), (1.e9,1.e9,1.e9), (0.,0.,0.), (0.,0.,0.)],
[((16,4),5.818e-4/pf,1.,1886.), ('PIC',GridPIC,0), (1.e9,1.e9,1.e9), (0.,0.,0.), (0.,0.,0.)],
)

out_folder = './../tests.out/test_remiCPC/'
simul_time = 1000.
dt = 0.036
diag_eb_step  = int(100./dt)
diag_dns_step = int(100./dt)
diag_phs_step = int(100./dt)
diag_mov_step = 100
