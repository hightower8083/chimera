import numpy as np
import sys
import matplotlib.pyplot as plt

inpt = __import__(sys.argv[1])
dt = inpt.dt

a0 = [inpt.seed_in['a0'],]
Qe = [0,]

TotSteps = int(sys.argv[2])

for diag in inpt.diags_in:
	if diag['Type']=='Fields': diag_freq_lasr = diag['Step']
	if diag['Type']=='Density': diag_freq_plas = diag['Step']
	if diag['Type']=='Particles': diag_freq_eons = diag['Step']
	
for i in np.arange(diag_freq_lasr,TotSteps+1,diag_freq_lasr):
	istr = str(i)
	while len(istr)<7: istr='0'+istr
	a0.append(np.abs(np.load('./data/ee_0_'+istr+'.npy')[:,:5,0,2]).max())

for i in np.arange(diag_freq_eons,TotSteps+1,diag_freq_eons):
	istr = str(i)
	while len(istr)<7: istr='0'+istr
	phs = np.load('./data/phs0_'+istr+'.npy')
	if phs.shape[-1]!=0:
		Qe.append((phs[-1]*(phs[3]>150)*( ( (0.8e-3*phs[1:3])**2).sum(0)<0.0025 )).sum()*0.8*178.62)
	else:
		Qe.append(0.0)

phs = None
a0 = np.array(a0)
Qe = np.array(Qe)

x_fld = ((diag_freq_lasr*np.arange(a0.shape[0]))*dt+inpt.seed_in['x0'])*0.8e-3
x_eons= ((diag_freq_eons*np.arange(Qe.shape[0]))*dt+inpt.seed_in['x0']-inpt.electrons_in['Density']**-0.5)*0.8e-3

xx = np.r_[0:TotSteps*dt*0.8e-3:200j]
aa = inpt.seed_in['a0']/(1+((xx-inpt.seed_in['x_foc'])**2/(np.pi*inpt.seed_in['LR']**2*0.8e-3)))**.5

plt.plot(xx,aa,'--',c='k')
plt.plot(x_fld,a0,lw=1.5,c='r')
plt.fill_between(xx,inpt.densprof(xx/0.8e-3)*(inpt.electrons_in['Density']*1.1e21/0.8**2*1e-19),alpha=0.3,color='b',interpolate=True)
plt.fill_between(0.5*(x_eons[1:]+x_eons[:-1]),-(Qe[1:]-Qe[:-1])/(x_eons[1]-x_eons[0])*1e-3,alpha=0.6,color='r',interpolate=True)
plt.ylim(0,plt.ylim()[1])
plt.xlim(0,xx[-1])
plt.legend(('$a_0$ (w/o plasma)','$a_0$','$n_e$ [$10^{19}$ cm$^{-3}$]','$dQ_e/dz|_{\gamma>150}$ [pC/$\mu$m]'),loc=2)
plt.xlabel('Propagation distance (mm)',fontsize=18)
plt.show()
