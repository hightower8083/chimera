import numpy as np

EE = 180.
lambd_r = 200e-9
lambd_u = 20e-3

def cox_transport(particles,particles_cntr):
	parts1 = np.zeros_like(particles[:6,:])
	parts1[:3] = (particles[:3].T - particles[:3].mean(-1)).T*lambd_u

	pz_mean = ((particles[[3,4,5]]**2).sum(0)**0.5).mean()
	parts1[3] = ((particles[[3,4,5]]**2).sum(0)**0.5 - pz_mean)/(particles[[3,4,5]]**2).sum(0)**0.5
	parts1[[4,5]] = particles[[4,5]]/particles[3]
	
	r11  =  10.
	r126 = -0.63*r11
	r522 = -0.3
	r56  =  0.5e-3
	L_bw = -1.

	g_loc = (EE+0.511)/0.511; brho = EE*1e6/3e8; 
	aw = (2*g_loc**2*lambd_r/lambd_u - 1.)**0.5; 
	B_loc = 2**.5*aw/93.4/lambd_u; K_loc = brho/B_loc; 
	K_loc = 1/(2.**.5*K_loc)

	parts1[0] = parts1[0] + r56*parts1[3] + r522*(parts1[4]**2+parts1[5]**2)
	parts1[[1,2]] = r11*parts1[[1,2]] + r126*parts1[3]*parts1[[4,5]]
	parts1[[4,5]] = parts1[[4,5]]/r11

	parts1[0] = parts1[0] + L_bw*parts1[3]/g_loc**2
	parts1[1] = parts1[1] + L_bw*parts1[4]
	tmp = +parts1[2]
	parts1[2] = np.cos(K_loc*L_bw)*tmp + np.sin(K_loc*L_bw)/K_loc*parts1[5]
	parts1[5] = -K_loc*np.sin(K_loc*L_bw)*tmp + np.cos(K_loc*L_bw)*parts1[5]

	particles[3] = pz_mean/(1.0 - parts1[3])
	particles[4] = parts1[4]*particles[3]
	particles[5] = parts1[5]*particles[3]
	particles[:3] = (parts1[:3].T/lambd_u+particles[:3].mean(-1)).T
#	particles[3] -= particles[0]/1e-3*7.2 # chirp remove r11  =  20. r56  = 1e-3
	particles[3] -= particles[0]*1.4e4 # chirp remove r11  =  10. r56  =  0.5e-3
	particles[6] = (1.+(particles[[3,4,5]]**2).sum(0))**0.5
	particles_cntr[:3] = particles[:3]

	return particles,particles_cntr

