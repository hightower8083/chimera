import numpy as np
from noinput import *

GridPIC = (-50.,0.,20.,0.04,0.25)
Solvers_in = ( [ 'PIC.1',GridPIC, (0.,1.e10),(1,((250.,1.),(30.,0.)))], )
RunWindows_in = ( (0.,1e9, 1.,10,(10., 0.)), )

Seeds_in = ( ['Gauss', 3., 1., (-15.,0.,0.,45.),(4.,4.)], )

PlasmaProfiles_in = (
('ramp_up' ,'x', (0 , 20.)),
('plateau' ,'x', (20., 1e9 )), 
('ramp_up' ,'x', (20., 40.)), 
('ramp_down','x', (40., 50.)), 
)

pf = 2. # prof fact
Particles_in = (
[((2,2,4),0.005/pf,-1.,1.), ('PIC',GridPIC,0), (1.e9,1.e9,1.e9), (0.,0.,0.), (0.,0.,0.)],
[((2,2,4),0.005/pf,1.,1886.), ('PIC',GridPIC,0), (1.e9,1.e9,1.e9), (0.,0.,0.), (0.,0.,0.)], 
)

out_folder = './../tests.out/small_bubble_test/'
simul_time = 100. 
dt = 0.04
diag_eb_step = int(20./dt) 
diag_dns_step = int(20./dt) 
diag_mov_step = 50
