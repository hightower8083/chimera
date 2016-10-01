import numpy as np
from scipy.special import jn_zeros,jn,j1
import chimera.moduls.fimera as chimera

class Solver:
	def __init__(self,solver_in):
		self.Configs = solver_in
		if 'Features' not in self.Configs: self.Configs['Features'] = ()

		Nko = self.Configs['MaxAzimuthMode']
		leftX,rightX, lengthR, dx,dr = self.Configs['Grid']
		dt = self.Configs['TimeStep']

		print 'Constructing a solver in the cylindric box', (leftX,rightX, lengthR)
		print 'spatial and temporal resolutions are', (dx,dr,dt)

		if 'TimeActive' in self.Configs:
			self.TimeInterval = self.Configs['TimeActive']
			print 'Active interval is limited to', self.TimeInterval
		else:
			self.TimeInterval = (0.0,np.inf)

		if 'KxShift' in self.Configs:
			kx0 = 2*np.pi*self.Configs['KxShift']
			print 'Spectral domain is shifted to', self.Configs['KxShift']
		else:
			kx0 = 0.0

		if 'Xchunked' in self.Configs:
			nthrds = self.Configs['Xchunked'][0]
			Nxchunk = 2*int(np.round(0.5/nthrds/dx*(rightX - leftX)))
			Nx = Nxchunk*nthrds
		else:
			Nx = int(np.round(0.5/dx*(rightX - leftX))*2)


		kx_env = 2*np.pi*np.fft.fftfreq(Nx,dx)
		kx     = kx0 + kx_env
		Xgrid  = rightX - dx*np.arange(Nx)[::-1]
		leftX = Xgrid[0]
		dkx = (kx[1]-kx[0])/(2.*np.pi)

		Nkr = int(np.round(lengthR/dr))
		Nr = Nkr+1
		RgridFull = dr*(np.arange(Nr)-0.5)
		lengthR = RgridFull[-1] + dr

		if 'KxShift' in self.Configs:
			Mmin    ,Mmax    ,Mtot     = -Nko  ,Nko,2*Nko+1
			Mmin_ext,Mmax_ext,Mtot_ext = -Nko-1,Nko+1,2*Nko+3
		else:
			Mmin    ,Mmax    ,Mtot     = 0,Nko,Nko+1
			Mmin_ext,Mmax_ext,Mtot_ext = 0,Nko+1,Nko+2

		print 'Grid resolutions are ', (Nx,Nkr,Mtot)

		kr   = np.zeros((Nkr,Mtot_ext))
		kr_g = np.zeros((Nx,Nkr,Mtot ))
		w    = np.zeros((Nx,Nkr,Mtot ))

		for jm in np.arange(Mmin_ext,Mmax_ext+1):
			kr[:,jm] = jn_zeros(jm,Nkr)/lengthR
		for jm in np.arange(Mmin,Mmax+1):
			kr_g[:,:,jm], kx_g = np.meshgrid(kr[:,jm],kx)
			w[:,:,jm] = np.sqrt(kx_g**2 + kr_g[:,:,jm]**2)

		Out    = np.zeros((Nkr,Nkr,Mtot))
		In     = np.zeros((Nkr,Nkr,Mtot))
		DpS2R  = np.zeros((Nkr,Nkr,Mtot))
		DmS2R  = np.zeros((Nkr,Nkr,Mtot))
		DpS2S  = np.zeros((Nkr,Nkr,Mtot))
		DmS2S  = np.zeros((Nkr,Nkr,Mtot))

		for jm in np.arange(Mmin,Mmax+1):
			Out[:,:,jm] = jn(jm, RgridFull[1:,None]*kr[:,jm][None,:])
			In [:,:,jm] = np.linalg.inv(Out[:,:,jm])
			DpS2R[:,:,jm] = 0.5*kr[:,jm+1][None,:]*jn(jm,RgridFull[1:,None]*kr[:,jm+1][None,:])
			DmS2R[:,:,jm] = 0.5*kr[:,abs(jm-1)][None,:]*jn(jm,RgridFull[1:,None]*kr[:,abs(jm-1)][None,:])

		for jm in np.arange(Mmin,Mmax+1):
			DpS2S[:,:,jm] = In[:,:,jm].dot(DpS2R[:,:,jm])
			DmS2S[:,:,jm] = In[:,:,jm].dot(DmS2R[:,:,jm])

		VGrid = 2*np.pi*dx*dr*RgridFull
		VGrid = (VGrid+(RgridFull==0))**-1*(RgridFull>0.0)
		InCurr = In*VGrid[None,1:,None]

		if ('KxShift' in self.Configs) and (Nko>0):
			Out    = np.concatenate((  Out[:,:,Mmin:],  Out[:,:,:Mmax+1]),axis=-1)
			In     = np.concatenate((   In[:,:,Mmin:],   In[:,:,:Mmax+1]),axis=-1)
			DpS2S  = np.concatenate((DpS2S[:,:,Mmin:],DpS2S[:,:,:Mmax+1]),axis=-1)
			DmS2S  = np.concatenate((DmS2S[:,:,Mmin:],DmS2S[:,:,:Mmax+1]),axis=-1)
			w      = np.concatenate((    w[:,:,Mmin:],    w[:,:,:Mmax+1]),axis=-1)
			kr_g   = np.concatenate(( kr_g[:,:,Mmin:], kr_g[:,:,:Mmax+1]),axis=-1)
			InCurr = np.concatenate((InCurr[:,:,Mmin:],InCurr[:,:,:Mmax+1]),axis=-1)

		if 'Rcut' in self.Configs:
			indRcut = (RgridFull<self.Configs['Rcut']).sum()
			Rgrid  = RgridFull[:indRcut].copy()
			InCurr = InCurr[:,:indRcut-1]
			OutFull=  Out.copy()
			Out    =  Out[:indRcut-1]
			print 'Rgrid is cut after ',self.Configs['Rcut']
		else:
			OutFull =  Out
			Rgrid  = RgridFull

		Nr = Rgrid.shape[0]
		kx_g = np.asfortranarray(kx_g)
		kr_g = np.asfortranarray(kr_g)
		w = np.asfortranarray(w)
		DpS2S = np.asfortranarray(np.swapaxes(DpS2S,0,1))
		DmS2S = np.asfortranarray(np.swapaxes(DmS2S,0,1))
		InCurr = np.asfortranarray(np.swapaxes(InCurr,0,1))
		In = np.asfortranarray(np.swapaxes(In,0,1))
		Out = np.asfortranarray(np.swapaxes(Out,0,1))

		self.Args = {'leftX':leftX,'rightX':rightX,'lowerR':(Rgrid*(Rgrid>=0)).min(),\
		  'upperR':Rgrid.max(),'kx_g':kx_g,'kr_g':kr_g,'w':w,'kx0':kx0,'kx_env':kx_env,\
		  'Nkr':Nkr,'Nx':Nx,'Nko':Nko, 'Xgrid':Xgrid,'Rgrid':Rgrid,'lengthR':lengthR,\
		  'dkx':dkx,'dt':dt,'dx':dx,'dr':dr,'kx':kx,'VGrid':VGrid,'RgridFull':RgridFull}

		if 'CoPropagative' in self.Configs:
			self.PSATD_coeffs(self.Configs['CoPropagative'])
		else:
			self.PSATD_coeffs()

		self.Args['FBDiff'] = (DpS2S,DmS2S,kx)
		self.Args['PoissFact'] = np.asfortranarray(w**-2)

		self.Args['EnergyFact'] = 0.5*0.511e6*1.6022e-19/2.818e-13*lengthR**2/dkx*\
		  jn(np.abs(np.arange(Mmin,Mmax+1)[None,None,:])+1,kr_g*lengthR)**2

		if 'KxShift' in self.Configs:
			cutafter = 0.8
			fu_bandpass = lambda x : (x<cutafter)+(x>=cutafter)*np.cos(np.pi/2*(x-cutafter)/(1-cutafter))**2
			filt_bandpass = fu_bandpass(np.abs(kx_env)/np.abs(kx_env.max()))[:,None,None]

			filt_antialias = np.ones_like(filt_bandpass)
			fu_antialias = lambda x,x0 :1-np.exp(-(x-x0)**2/((x0+1.)/Nx)**2)
			cell_echos = np.abs(kx_env).max()/kx0*np.arange(20)-1.
			full_band = np.array([kx.min()/kx0, kx.max()/kx0])-1.

			for cellecho in cell_echos:
				if cellecho>full_band[0] and cellecho<full_band[1]:
					print 'possible grid echo is detected at', cellecho/abs(full_band).max()
					if 'NoAntiEcho' in self.Configs['Features']:
						continue
					elif abs(cellecho)/abs(full_band).max()>0.4:
						print 'will correct', cellecho/abs(full_band).max()
						filt_antialias *= fu_antialias( kx/kx0-1, cellecho  )[:,None,None]
					else:
						print 'echo is close to resonance; no correction will be performed'

			self.Args['CurrFact'] = np.asfortranarray((2*np.pi)**2/Nx*\
			  np.cos(0.5*np.pi*kr_g/kr_g.max(0).max(0))**2*filt_bandpass*filt_antialias)
			self.Args['DepProj'] = (Rgrid,1./dx,1./dr,kx0)
			self.Args['FBCurrIn'] = (kx_env,InCurr)
			self.Args['FBIn']     = (kx_env,In)
			self.Args['FBout']    = (kx_env,Out)
			self.Args['FBoutFull']    = (kx_env,OutFull)
		else:
			self.Args['CurrFact'] = np.asfortranarray( (2*np.pi)**2/Nx*\
			  np.cos(0.5*np.pi*kr_g/kr_g.max(0).max(0))**2*\
			  np.cos(0.5*np.pi*np.abs(kx_g)/np.abs(kx_g.max()))[:,:,None]**2)
			self.Args['DepProj'] = (Rgrid,1./dx,1/dr)
			self.Args['FBCurrIn']  = (kx,InCurr)
			self.Args['FBIn']      = (kx,In)
			self.Args['FBout']     = (kx,Out)
			self.Args['FBoutFull']     = (kx,OutFull)

		self.EB      = np.zeros((Nx,Nr,Mtot,6),dtype='complex',order='F')
		self.vec_spc = np.zeros((Nx,Nr,Mtot,3),dtype='complex',order='F')
		self.scl_spc = np.zeros((Nx,Nr,Mtot  ),dtype='complex',order='F')

		self.EG_fb      = np.zeros((Nx,Nkr,Mtot,6),dtype='complex',order='F')
		self.vec_fb     = np.zeros((Nx,Nkr,Mtot,3),dtype='complex',order='F')
		self.vec_fb_aux = np.zeros((Nx,Nkr,Mtot,3),dtype='complex',order='F')
		self.scl_fb     = np.zeros((Nx,Nkr,Mtot  ),dtype='complex',order='F')

		if 'SpaceCharge' in self.Configs['Features'] or 'StaticKick' in self.Configs['Features']:
			print 'Space charge is added'
			self.vec_fb_aux0 = np.zeros((Nx,Nkr,Mtot,3),dtype='complex',order='F')
			self.vec_fb_aux1 = np.zeros((Nx,Nkr,Mtot,3),dtype='complex',order='F')
		if 'StillAsBackground' in self.Configs['Features']:
			print '"Still" species as are treated as background'
			self.bg_spc = np.zeros_like(self.scl_spc)

		if 'NoPoissonCorrection' in self.Configs['Features']: 
			print  'Poisson correction will not be performed'

	def PSATD_coeffs(self,beta=1):
		kx_g,w = self.Args['kx_g'],self.Args['w']
		kx_g = beta*kx_g[:,:,None]
		dt = self.Configs['TimeStep']
		if 'SpaceCharge' in self.Configs['Features']:
			self.CPSATD1 = np.zeros(w.shape+(5,),dtype='double',order='F')
			self.CPSATD2 = np.zeros(w.shape+(5,),dtype='double',order='F')
		elif 'KxShift' in self.Configs:
			self.CPSATD1 = np.zeros(w.shape+(3,),dtype='complex',order='F')
			self.CPSATD2 = np.zeros(w.shape+(3,),dtype='complex',order='F')
		else:
			self.CPSATD1 = np.zeros(w.shape+(3,),dtype='double',order='F')
			self.CPSATD2 = np.zeros(w.shape+(3,),dtype='double',order='F')

		self.CPSATD1[:,:,:,0] = np.cos(dt*w)
		self.CPSATD1[:,:,:,1] = np.sin(dt*w)/w
		self.CPSATD2[:,:,:,0] = -w*np.sin(dt*w)
		self.CPSATD2[:,:,:,1] = np.cos(dt*w)

		if 'KxShift' in self.Configs:
			self.CPSATD1[:,:,:,2] = np.exp(-0.5j*kx_g*dt)*1j*kx_g/(w**2-kx_g**2)*\
			  (1.0-np.exp(1j*kx_g*dt)*(w/(1j*kx_g)*np.sin(dt*w)+np.cos(dt*w)))
			self.CPSATD2[:,:,:,2]  = np.exp(-0.5j*kx_g*dt)*w**2/(w**2-kx_g**2)*\
			  (1.0+np.exp(1j*kx_g*dt)*(1j*kx_g/w*np.sin(dt*w)-np.cos(dt*w)))
		else:
			self.CPSATD1[:,:,:,2] = -np.sin(dt*w)/w
			self.CPSATD2[:,:,:,2] = 1-np.cos(dt*w)

		if 'SpaceCharge' in self.Configs['Features']:

#			self.CPSATD1[:,:,:,3] = 0.5*(np.cos(dt*w)-1)/w**2
#			self.CPSATD1[:,:,:,4] = 0.5*(np.cos(dt*w)-1)/w**2

#			self.CPSATD2[:,:,:,3] = -0.5*np.sin(dt*w)/w
#			self.CPSATD2[:,:,:,4] = -0.5*np.sin(dt*w)/w

			self.CPSATD1[:,:,:,3] = (dt*w*np.cos(dt*w)-np.sin(dt*w))/w**3/dt
			self.CPSATD1[:,:,:,4] = (np.sin(dt*w)-dt*w)/w**3/dt

			self.CPSATD2[:,:,:,3] = (1-np.cos(dt*w)-dt*w*np.sin(dt*w))/w**2/dt
			self.CPSATD2[:,:,:,4] = (np.cos(dt*w)-1)/w**2/dt

	def maxwell_solver_init(self,px0):
		if 'SpaceCharge' not in self.Configs['Features'] and 'StaticKick' not in self.Configs['Features']: return
		beta0 = px0/np.sqrt(1+px0**2)
		kx_g,w = self.Args['kx_g'],self.Args['w']
		CPSATD1 = np.zeros(w.shape+(2,),dtype='complex',order='F')
		CPSATD2 = np.zeros(w.shape+(2,),dtype='complex',order='F')

		kx_g = beta0*kx_g[:,:,None]
		CPSATD1[:,:,:,0] = 1.j*kx_g/(w**2-kx_g**2)
		CPSATD1[:,:,:,1] = -1./(w**2-kx_g**2)
		CPSATD2[:,:,:,0] = w**2/(w**2-kx_g**2)
		CPSATD2[:,:,:,1] = 1.j*kx_g/(w**2-kx_g**2)
		self.EG_fb = chimera.maxwell_init_push(self.EG_fb,self.vec_fb,self.vec_fb_aux1,CPSATD1,CPSATD2)

	def maxwell_solver(self):
		if 'SpaceCharge' in self.Configs['Features']:
			self.EG_fb = chimera.maxwell_push_with_spchrg(self.EG_fb,self.vec_fb,self.vec_fb_aux0,self.vec_fb_aux1,self.CPSATD1,self.CPSATD2)
		else:
			self.EG_fb = chimera.maxwell_push_wo_spchrg(self.EG_fb,self.vec_fb,self.CPSATD1,self.CPSATD2)

	def fb_curr_in(self):
		self.vec_fb = chimera.fb_vec_in(self.vec_fb,self.vec_spc,self.Args['leftX'],*self.Args['FBCurrIn'])
		self.vec_fb = chimera.omp_mult(self.vec_fb,self.Args['CurrFact'])

	def fb_dens_in(self):
		self.scl_fb = chimera.fb_scl_in(self.scl_fb,self.scl_spc,self.Args['leftX'],*self.Args['FBCurrIn'])
		self.FBGradDens()
		self.vec_fb_aux1 = chimera.omp_mult(self.vec_fb_aux1,self.Args['CurrFact'])

	def fb_scl_spc_in(self):
		self.scl_fb = chimera.fb_scl_in(self.scl_fb,self.scl_spc,self.Args['leftX'],*self.Args['FBIn'])

	def fb_fld_out(self):
		self.EB = chimera.fb_eb_out(self.EB,self.EG_fb[:,:,:,:3],self.vec_fb_aux,\
		  self.Args['leftX'],*self.Args['FBout'])

		self.EB /= 2*np.pi
		if 'KxShift' in self.Configs:
			self.EB *= 2.0
		else:
			self.EB[:,:,1:] *= 2.0

		if 'KxShift' in self.Configs:
			for mm in range(self.EB.shape[2]):
				if mm==self.Args['Nko']:
					self.EB[:,0,mm] = self.EB[:,1,mm]
				else:
					self.EB[:,0,mm] = -self.EB[:,1,mm]
		else:
			self.EB[:,0,0] = self.EB[:,1,0]
			self.EB[:,0,1:] = -self.EB[:,1,1:]

	def poiss_corr(self):
		if 'NoPoissonCorrection' in self.Configs['Features']: return
		self.vec_fb_aux[:] = self.vec_fb
		self.FBDiv()
		self.FBGrad()
		if 'SpaceCharge' in self.Configs['Features']:
			self.vec_fb = chimera.poiss_corr_with_spchrg(self.vec_fb,self.vec_fb_aux,\
			  self.vec_fb_aux0,self.vec_fb_aux1,self.Args['PoissFact'],1./self.Configs['TimeStep'])
		else:
			self.vec_fb = chimera.poiss_corr_wo_spchrg(self.vec_fb,self.vec_fb_aux,self.Args['PoissFact'])

	def FBGrad(self):
		if 'KxShift' in self.Configs:
			self.vec_fb_aux = chimera.fb_grad_env(self.vec_fb_aux,self.scl_fb,*self.Args['FBDiff'])
		else:
			self.vec_fb_aux = chimera.fb_grad(self.vec_fb_aux,self.scl_fb,*self.Args['FBDiff'])

	def FBDiv(self):
		if 'KxShift' in self.Configs:
			self.scl_fb = chimera.fb_div_env(self.scl_fb,self.vec_fb_aux,*self.Args['FBDiff'])
		else:
			self.scl_fb = chimera.fb_div(self.scl_fb,self.vec_fb_aux,*self.Args['FBDiff'])

	def FBGradDens(self):
		if 'KxShift' in self.Configs:
			self.vec_fb_aux1 = chimera.fb_grad_env(self.vec_fb_aux1,self.scl_fb,*self.Args['FBDiff'])
		else:
			self.vec_fb_aux1 = chimera.fb_grad(self.vec_fb_aux1,self.scl_fb,*self.Args['FBDiff'])

	def G2B_FBRot(self):
		if 'KxShift' in self.Configs:
			self.vec_fb_aux = chimera.fb_rot_env(self.vec_fb_aux,self.EG_fb[:,:,:,3:],*self.Args['FBDiff'])
		else:
			self.vec_fb_aux = chimera.fb_rot(self.vec_fb_aux,self.EG_fb[:,:,:,3:],*self.Args['FBDiff'])
		self.vec_fb_aux = chimera.omp_mult(self.vec_fb_aux, self.Args['PoissFact'])

	def B2G_FBRot(self):
		if 'KxShift' in self.Configs:
			self.EG_fb[:,:,:,3:] = chimera.fb_rot_env(self.EG_fb[:,:,:,3:],self.vec_fb_aux,*self.Args['FBDiff'])
		else:
			self.EG_fb[:,:,:,3:] = chimera.fb_rot(self.EG_fb[:,:,:,3:],self.vec_fb_aux,*self.Args['FBDiff'])

	def add_gauss_beam(self,S):
		k0 = 2*np.pi*S['k0']
		a0 = 2*np.pi*S['a0']
		X_focus = S['x0']-S['x_foc']
		kx_g, kr_g,w = self.Args['kx_g'],self.Args['kr_g'],self.Args['w']
		w = w[:,:,:,None]

		if 'KxShift' in self.Configs:
			e_s0 = a0*0.5*(np.pi)**0.5*S['Lx']*S['LR']**2*self.Args['dkx']/self.Args['lengthR']**2
			self.vec_fb_aux[:,:,self.Args['Nko'],2] = e_s0/j1(self.Args['lengthR']*kr_g[:,:,self.Args['Nko']])**2*\
			  np.exp(-1j*kx_g*S['x0'])*\
			  np.exp(-0.25*(kx_g-k0)**2*S['Lx']**2 -0.25*kr_g[:,:,self.Args['Nko']]**2*S['LR']**2 )
			DT = -1.j*w
		else:
			Xgrid,Rgrid = self.Args['Xgrid'],self.Args['Rgrid']	# sin laser phase
			self.scl_spc[:,:,0] = a0*np.sin(k0*(Xgrid[:,None]-S['x0']))*\
			  np.exp(-(Xgrid[:,None]-S['x0'])**2/S['Lx']**2-Rgrid[None,:]**2/S['LR']**2)*\
			  (abs(Xgrid[:,None]-S['x0'])<3.5*S['Lx'])*(abs(Rgrid[None,:])<3.5*S['LR'])
			self.scl_spc[:,0,0] = 0.0
			self.fb_scl_spc_in()
			self.vec_fb_aux[:,:,:,2] = self.scl_fb/np.float(self.Args['Nx'])
			DT = -1.j*w*np.sign(kx_g[:,:,None,None] + (kx_g[:,:,None,None]==0))

		EE = self.vec_fb_aux.copy()
		self.FBDiv()
		self.FBGrad()
		self.vec_fb_aux = chimera.omp_mult(self.vec_fb_aux, self.Args['PoissFact'])
		EE += self.vec_fb_aux
		GG  = DT*EE

		self.vec_fb_aux[:] = np.cos(w*X_focus)*EE + np.sin(w*X_focus)/w*GG
		GG = -w*np.sin(w*X_focus)*EE + np.cos(w*X_focus)*GG
		EE = self.vec_fb_aux.copy()
		EE *= np.exp(1.j*kx_g[:,:,None,None]*X_focus)
		GG *= np.exp(1.j*kx_g[:,:,None,None]*X_focus)

		self.EG_fb[:,:,:,:3] += EE
		self.EG_fb[:,:,:,3:] += GG
		self.vec_fb_aux[:] = 0.0
		self.scl_fb[:] = 0.0

	def absorb_field(self,Lf,config='left'):
		Nfilt = int(Lf/self.Args['dx'])
		filt_shape = (0.5-0.5*np.cos(np.r_[0:np.pi:Nfilt*1j]))**2
		filt = np.ones(self.Args['Nx'])
		if config=='left':
			filt[:Nfilt] = filt_shape
		elif config=='right':
			filt[-Nfilt:] = filt_shape[::-1]
		elif config=='both':
			filt[:Nfilt] = filt_shape
			filt[-Nfilt:] = filt_shape[::-1]
		self.EG_fb[:,:,:,:3] = chimera.fb_filtr(self.EG_fb[:,:,:,:3],self.Args['leftX'],\
		  self.Args['kx'],filt)
		self.vec_fb_aux = chimera.fb_filtr(self.vec_fb_aux,self.Args['leftX'],\
		  self.Args['kx'],filt)
		self.B2G_FBRot()


	def FBRot(self,vec_in):
		if 'KxShift' in self.Configs:
			return chimera.fb_rot_env(np.empty_like(self.vec_fb),vec_in,*self.Args['FBDiff'])
		else:
			return chimera.fb_rot(np.empty_like(self.vec_fb),vec_in,*self.Args['FBDiff'])
