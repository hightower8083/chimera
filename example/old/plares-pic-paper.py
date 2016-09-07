import numpy as np
from noinput import *

GridPIC = (-55.,3.,32.,0.02,0.2)

Solvers_in = (
[ 'PIC.1',GridPIC, (0.,1.e10),(1,((250.,1.),(30.,0.)),1)],
)

Seeds_in = (
['Gauss', 3., 1., (-15.,0.,0.,300.),(5.,4.,1.)],
)

pf = 2. # profile factor

Particles_in = (
[((16,4),0.005/pf,-1.,1.),
('PIC',GridPIC,0),
(1.e9,1.e9,1.e9),
(0.,0.,0.),
(0.,0.,0.)],
[((16,4),0.005/pf,1.,1886.),
('PIC',GridPIC,0),
(1.e9,1.e9,1.e9),
(0.,0.,0.),
(0.,0.,0.)],
)

PlasmaProfiles_in = (
('ramp_up'  ,'x', (0   , 150.)),
('plateau'  ,'x', (150., 1e9 )),
('ramp_up'  ,'x', (150., 300.)),
('ramp_down','x', (300., 315.)),
)

RunWindows_in = (
(0.,1e9, 1.,10,(7., 0.)),
)

out_folder = '/work/sources/magins/andriyash/code/tests.out/4paper_02_01e/'
simul_time = 500.
dt = 0.02
diag_eb_step  = int(50./dt)
diag_dns_step = int(50./dt)
diag_phs_step = int(100./dt)
diag_mov_step = 30

