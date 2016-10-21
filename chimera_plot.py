import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import medfilt,medfilt2d

mymap = plt.cm.Reds;mymap._init();mymap._lut[:,-1] = abs(np.sin(np.r_[0:.5*np.pi:259j]))
mymap1 = plt.cm.seismic;mymap1._init();mymap1._lut[:,-1] = abs(np.sin(np.r_[-0.5*np.pi:0.5*np.pi:259j]))
mymap2 = plt.cm.PiYG;mymap2._init();mymap2._lut[:,-1] = abs(np.sin(np.r_[-0.5*np.pi:0.5*np.pi:259j]))**2
Ntheta = 60

inpt = __import__(sys.argv[1])
specie_in = inpt.electrons_in
BoxGrid = inpt.BoxGrid
PlasmaGrid = inpt.PlasmaGrid
dt =inpt.dt
MovingFrame = inpt.MovingFrame
solver = inpt.solver
out_folder = inpt.out_folder
if 'Velocity' in MovingFrame:
	vb = MovingFrame['Velocity']
else:
	vb=1.
comp = sys.argv[2]

def dens_plot(t,th=0,vmax = 3,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	n_max = vmax*specie_in['Density']
	extnt = np.array([PlasmaGrid[0],PlasmaGrid[1],-PlasmaGrid[2],PlasmaGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8e-3
	xlim *= 0.8e-3
	plt.clf()
	dat = np.load(out_folder+'dens_'+comp+'_'+tstr+'.npy')
	phase = np.exp(-0.5j*np.pi*th*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]-1))
	for i in xrange(len(phase)):
		dat0 += 0.5*(3-(-1)**i)*np.real(phase[i]*np.hstack((dat[:,::-1,i],(-1)**i*dat[:,1:,i])))
	plt.imshow(dat0.T, origin='lower',aspect='auto',vmin=-n_max,vmax=0.,extent=extnt,**auxargs)
	plt.xlim(xlim)

def fld_plot(t,th=0,s=2,vmax=None,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	extnt = np.array([BoxGrid[0],BoxGrid[1],-BoxGrid[2],BoxGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8e-3
	xlim *= 0.8e-3
#	plt.clf()
	dat = np.load(out_folder+'ee_'+comp+'_'+tstr+'.npy')
	phase = np.exp(-0.5j*np.pi*th*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(len(phase)):
		dat0 += np.real(phase[i]*np.hstack((dat[:,::-1,i,s],(-1)**i*dat[:,:,i,s])))
	if vmax==None:
		vmax=np.abs(dat0).max()
		vmin=-vmax
	else:
		vmin=-vmax
	plt.imshow(dat0.T, origin='lower',aspect='auto',vmax=vmax,vmin=vmin,extent=extnt,**auxargs)
	plt.xlim(xlim)

def envel_plot(t,th=0,s=2,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	extnt = np.array([BoxGrid[0],BoxGrid[1],-BoxGrid[2],BoxGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8e-3
	xlim *= 0.8e-3

	dat = np.load(out_folder+'ee_'+comp+'_'+tstr+'.npy')
	phase = np.exp(-0.5j*np.pi*th*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(len(phase)):
		dat0 += np.real(phase[i]*np.hstack((dat[:,::-1,i,s],(-1)**i*dat[:,:,i,s])))
	for i in np.arange(dat0.shape[1]): dat0[:,i] = medfilt(np.abs(dat0[:,i]), 4*int(0.5/inpt.dt)+1 )
	plt.imshow(dat0.T, origin='lower',aspect='auto',extent=extnt,cmap=mymap,**auxargs)
	plt.xlim(xlim)

def envel_contour(t,th=0,s=2,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	extnt = np.array([BoxGrid[0],BoxGrid[1],-BoxGrid[2],BoxGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8e-3
	xlim *= 0.8e-3

	dat = np.load(out_folder+'ee_'+comp+'_'+tstr+'.npy')
	phase = np.exp(-0.5j*np.pi*th*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(len(phase)):
		dat0 += np.real(phase[i]*np.hstack((dat[:,::-1,i,s],(-1)**i*dat[:,:,i,s])))
	#for i in np.arange(dat0.shape[1]): dat0[:,i] = medfilt(np.abs(dat0[:,i]), 4*int(0.5/inpt.dt)+1 )
	#x,y = np.mgrid[extnt[0]:extnt[1]:1j*dat0.shape[1],extnt[0]:extnt[1]:1j*dat0.shape[0]]
	dat0 = medfilt(np.abs(dat0),(3*int(0.5/inpt.dt)+1,1))
	plt.contour(dat0.T, origin='lower',aspect='auto',extent=extnt,cmap=mymap,**auxargs)
	plt.xlim(xlim)

def get_pwrO():
	pwrO = np.loadtxt(out_folder+'pwrO_'+comp+'.txt')
	return pwrO.reshape(pwrO.shape[0],Ntheta,pwrO.shape[1]/Ntheta)

def pwrO_plot(pwrO):
	Rgrid = solver.Args['RgridFull'][1:]
	Rgrid = Rgrid - Rgrid[0]
	xO = Rgrid[None,:]*np.sin(np.r_[0:2*np.pi:Ntheta*1j])[:,None]
	yO = Rgrid[None,:]*np.cos(np.r_[0:2*np.pi:Ntheta*1j])[:,None]
	plt.pcolormesh(xO,yO,pwrO,cmap='hot',shading='gouraud')
