import numpy as np
from scipy.constants import c

def cox_transport(particles,particles_cntr, Tr):

	parts1 = np.zeros_like(particles[:6,:])
	ww = np.abs(particles[-1])
	parts1[:3] = particles[:3]*Tr['lambda_u']

	pz_mean = (particles[3]*ww).sum()/ww.sum()
	parts1[3] = (particles[3]-pz_mean)/pz_mean
	parts1[4:] = particles[4:6]/particles[3]
	
	if 'r22' not in Tr: Tr['r22']=Tr['r11']
	r116 = -0.63*Tr['r11']
	r226 = -0.63*Tr['r22']
	r522 = -0.3

	parts1[0] = parts1[0] + Tr['r56']*parts1[3] + r522*(parts1[4]**2+parts1[5]**2)
	parts1[1] = Tr['r11']*parts1[1] + r116*parts1[3]*parts1[4]
	parts1[2] = Tr['r22']*parts1[2] + r226*parts1[3]*parts1[5]
	parts1[[4,5]] = parts1[[4,5]]/Tr['r11']

	g_loc = (Tr['EE']+0.511)/0.511
	B_loc = Tr['K0']/93.4/Tr['lambda_u']
	K_loc = c*B_loc/(Tr['EE']*1e6*2.**.5)

	parts1[0] = parts1[0] - Tr['foc']*parts1[3]/g_loc**2
	parts1[2] = parts1[2] - Tr['foc']*parts1[5]
 	tmp = +parts1[1]
	parts1[1] = np.cos(K_loc*Tr['foc'])*tmp - np.sin(K_loc*Tr['foc'])/K_loc*parts1[4]
	parts1[4] = K_loc*np.sin(K_loc*Tr['foc'])*tmp + np.cos(K_loc*Tr['foc'])*parts1[4]

	if 'x0' in Tr:
		parts1[0  ] = parts1[0  ] + Tr['x0']*parts1[3  ]/g_loc**2
		parts1[1:3] = parts1[1:3] + Tr['x0']*parts1[4:6]

	particles[:3] = parts1[:3]/Tr['lambda_u']
	particles[3] = pz_mean*(1.0+parts1[3])
	particles[4] = parts1[4]*particles[3]
	particles[5] = parts1[5]*particles[3]
	particles[6] = (1.+(particles[[3,4,5]]**2).sum(0))**0.5

	particles_cntr[:3] = particles[:3]

	return particles,particles_cntr
