import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.signal import medfilt,medfilt2d

#mymap = plt.cm.Reds;mymap._init();mymap._lut[:,-1] = abs(np.sin(np.r_[0:.5*np.pi:259j]))
#mymap1 = plt.cm.seismic;mymap1._init();mymap1._lut[:,-1] = abs(np.sin(np.r_[-0.5*np.pi:0.5*np.pi:259j]))
#mymap2 = plt.cm.PiYG;mymap2._init();mymap2._lut[:,-1] = abs(np.sin(np.r_[-0.5*np.pi:0.5*np.pi:259j]))**2
Ntheta = 60

inpt = __import__(sys.argv[1])
solver=inpt.solver
electrons = inpt.electrons
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
	extnt = np.array([electrons.Args['leftX'],electrons.Args['rightX'],-electrons.Args['upperR'],electrons.Args['upperR']])
	#extnt = np.array([PlasmaGrid[0],PlasmaGrid[1],-PlasmaGrid[2],PlasmaGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8 #e-3
	xlim *= 0.8 #e-3
	plt.clf()
	dat = np.load(out_folder+'dens_'+comp+'_'+tstr+'.npy')
	phase_p = np.exp(1.j*th*np.arange(dat.shape[2]))
	phase_m = np.exp(1.j*(th+np.pi)*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(dat.shape[2]):
		dat0 += np.real(np.hstack((dat[:,::-1,i]*phase_m[i],dat[:,:,i]*phase_p[i] )))
	plt.imshow(dat0.T, origin='lower',aspect='auto',vmin=-n_max,vmax=0.,extent=extnt,**auxargs)
	plt.xlim(xlim)
	return dat0, extnt

def fld_plot(t,th=0,s=2,vmax=None,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	extnt = np.array([solver.Args['leftX'],solver.Args['rightX'],-solver.Args['upperR'],solver.Args['upperR']])
#	extnt = np.array([BoxGrid[0],BoxGrid[1],-BoxGrid[2],BoxGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
	extnt *= 0.8 #e-3
	xlim *= 0.8 #e-3
	dat = np.load(out_folder+'ee_'+comp+'_'+tstr+'.npy')
	phase_p = np.exp(1.j*th*np.arange(dat.shape[2]))
	phase_m = np.exp(1.j*(th+np.pi)*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(dat.shape[2]):
		 dat0 += np.real(np.hstack((dat[:,::-1,i,s]*phase_m[i],dat[:,:,i,s]*phase_p[i] )))
	if vmax==None:
		vmax=np.abs(dat0).max()
		vmin=-vmax
	else:
		vmin=-vmax
	plt.imshow(dat0.T, origin='lower',aspect='auto',vmax=vmax,vmin=vmin,extent=extnt,**auxargs)
	plt.xlim(xlim)
	return dat0, extnt

def envel_plot(t,th=0,s=2,**auxargs):
	tstr = str(t)
	while len(tstr)<7: tstr='0'+tstr
	extnt = np.array([BoxGrid[0],BoxGrid[1],-BoxGrid[2],BoxGrid[2]])
	extnt[:2] += t*dt*vb
	if 'AbsorbLayer' in MovingFrame:
		xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
	else:
		xlim = np.array([extnt[0], extnt[1]])
#	extnt *= 0.8e-3
#	xlim *= 0.8e-3

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
#	extnt *= 0.8e-3
#	xlim *= 0.8e-3
	dat = np.load(out_folder+'ee_'+comp+'_'+tstr+'.npy')
	phase = np.exp(-0.5j*np.pi*th*np.arange(dat.shape[2]))
	dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
	for i in xrange(len(phase)):
		dat0 += np.real(phase[i]*np.hstack((dat[:,::-1,i,s],(-1)**i*dat[:,:,i,s])))
	dat0 = medfilt(np.abs(dat0),(4*int(0.5/inpt.dt)+1,1))
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

class iplots:
	def __init__(self,t,**auxargs):
		self.press  = False
		self.press_hold = False
		self.dens_plot_cont(t,**auxargs)
		self.connect()

	def connect(self):
		'connect to all the events we need'
		self.cidpress = self.fig.canvas.mpl_connect(
		  'button_press_event', self.on_press)
		self.cidrelease = self.fig.canvas.mpl_connect(
		  'button_release_event', self.on_release)
		self.cidmotion = self.fig.canvas.mpl_connect(
		  'motion_notify_event', self.on_motion)

	def dens_plot_cont(self,t,th=0,vmax=3,**auxargs):
		tstr = str(t)
		while len(tstr)<7: tstr='0'+tstr
		n_max = vmax*specie_in['Density']
		self.val = n_max
		self.val_on_motion = n_max

		extnt = np.array([PlasmaGrid[0],PlasmaGrid[1],-PlasmaGrid[2],PlasmaGrid[2]])
		extnt[:2] += t*dt*vb
		if 'AbsorbLayer' in MovingFrame:
			xlim = np.array([extnt[0]+MovingFrame['AbsorbLayer'], extnt[1]])
		else:
			xlim = np.array([extnt[0], extnt[1]])
		extnt *= 0.8e-3
		xlim *= 0.8e-3
		dat = np.load(out_folder+'dens_'+comp+'_'+tstr+'.npy')

		phase_p = np.exp(1.j*th*np.arange(dat.shape[2]))
		phase_m = np.exp(1.j*(th+np.pi)*np.arange(dat.shape[2]))
		self.dat0 = np.zeros((dat.shape[0], 2*dat.shape[1]))
		for i in xrange(dat.shape[2]):
			self.dat0 += np.real(np.hstack((dat[:,::-1,i]*phase_m[i],dat[:,:,i]*phase_p[i] )))	

		self.fig, ax = plt.subplots(1,1)
		self.pl = ax.imshow(np.abs(self.dat0).T, origin='lower',aspect='auto',vmin=0,vmax=self.val,**auxargs)
		self.datshape = self.dat0.shape

	def on_press(self,event):
		self.ix_press, self.iy_press=int(event.xdata), int(event.ydata)
		if event.button ==3:
			self.press_hold = True
		else:
			self.val = np.abs(self.dat0[self.ix_press, self.iy_press])
			self.press = True
			self.pl.set_clim(0,self.val)
			plt.draw()

	def on_release(self,event):
		self.press  = False
		self.press_hold = False
		self.val = self.val_on_motion

	def on_motion(self,event):
		if self.press_hold:
			self.val_on_motion = self.val*(1-(self.iy_press-int(event.ydata))/self.datshape[1])
			self.pl.set_clim(0,self.val_on_motion)
			plt.draw()

#print self.ix_press-int(event.xdata), self.iy_press-int(event.ydata)
