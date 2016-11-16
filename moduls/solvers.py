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

####################
		if 'KxShift' in self.Configs:
			Mmin    , Mmax    , Mtot     = -Nko-1, Nko+1 , 2*Nko+3
			Mmin_ext, Mmax_ext, Mtot_ext = Mmin-1, Mmax+1, Mtot+2
		else:
			Mmin    , Mmax    , Mtot     = 0, Nko+1 , Nko+2
			Mmin_ext, Mmax_ext, Mtot_ext = 0, Mmax+1, Mtot+1

		kr   = np.zeros((Nkr,Mtot_ext))
		kr_g = np.zeros((Nx,Nkr,Mtot ))
		w    = np.zeros((Nx,Nkr,Mtot ))
		DpS2S  = np.zeros((Nkr,Nkr,Mtot))
		DmS2S  = np.zeros((Nkr,Nkr,Mtot))
		for jm in np.arange(Mmin_ext,Mmax_ext+1):
			kr[:,jm] = jn_zeros(jm,Nkr)/lengthR
		for jm in np.arange(Mmin,Mmax+1):
			kr_g[:,:,jm], kx_g = np.meshgrid(kr[:,jm],kx)
			w[:,:,jm] = np.sqrt(kx_g**2 + kr_g[:,:,jm]**2)
			In = np.linalg.pinv(jn(jm, RgridFull[1:,None]*kr[:,jm][None,:]))
			DpS2S[:,:,jm] = 0.5*kr[:,jm+1][None,:]*In.dot(jn(jm,RgridFull[1:,None]*kr[:,jm+1][None,:]))
			DmS2S[:,:,jm] = 0.5*kr[:,abs(jm-1)][None,:]*In.dot(jn(jm,RgridFull[1:,None]*kr[:,abs(jm-1)][None,:]))
		if ('KxShift' in self.Configs) and (Nko>0):
			DpS2S  = np.concatenate((DpS2S[:,:,Mmin:],DpS2S[:,:,:Mmax+1]),axis=-1)
			DmS2S  = np.concatenate((DmS2S[:,:,Mmin:],DmS2S[:,:,:Mmax+1]),axis=-1)
			w      = np.concatenate((    w[:,:,Mmin:],    w[:,:,:Mmax+1]),axis=-1)
		w_ext = w.copy()
###################

		if 'KxShift' in self.Configs:
			Mmin    , Mmax    , Mtot     = -Nko  , Nko   , 2*Nko+1
			Mmin_ext, Mmax_ext, Mtot_ext = Mmin-1, Mmax+1, Mtot+2
		else:
			Mmin    , Mmax    , Mtot     = 0, Nko   , Nko+1
			Mmin_ext, Mmax_ext, Mtot_ext = 0, Mmax+1, Mtot+1

		print 'Grid resolutions are ', (Nx,Nkr,Mtot)

		kr   = np.zeros((Nkr,Mtot_ext))
		kr_g = np.zeros((Nx,Nkr,Mtot ))
		w    = np.zeros((Nx,Nkr,Mtot ))
		Out    = np.zeros((Nkr,Nkr,Mtot))
		In     = np.zeros((Nkr,Nkr,Mtot))

		for jm in np.arange(Mmin_ext,Mmax_ext+1):
			kr[:,jm] = jn_zeros(jm,Nkr)/lengthR
		for jm in np.arange(Mmin,Mmax+1):
			kr_g[:,:,jm], kx_g = np.meshgrid(kr[:,jm],kx)
			w[:,:,jm] = np.sqrt(kx_g**2 + kr_g[:,:,jm]**2)

		for jm in np.arange(Mmin,Mmax+1):
			Out[:,:,jm] = jn(jm, RgridFull[1:,None]*kr[:,jm][None,:])
			In [:,:,jm] = np.linalg.pinv(Out[:,:,jm])

		if ('KxShift' in self.Configs) and (Nko>0):
			Out    = np.concatenate((  Out[:,:,Mmin:],  Out[:,:,:Mmax+1]),axis=-1)
			In     = np.concatenate((   In[:,:,Mmin:],   In[:,:,:Mmax+1]),axis=-1)
			w      = np.concatenate((    w[:,:,Mmin:],    w[:,:,:Mmax+1]),axis=-1)
			kr_g   = np.concatenate(( kr_g[:,:,Mmin:], kr_g[:,:,:Mmax+1]),axis=-1)
			kr     = np.concatenate(( kr  [:  ,Mmin:], kr  [:  ,:Mmax+1]),axis=-1)

		VGrid = 2*np.pi*dx*dr*RgridFull
		VGrid = (VGrid+(RgridFull==0))**-1*(RgridFull>0.0)
		InCurr = In*VGrid[None,1:,None]

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

		self.Args['FBDiff'] = (DpS2S,DmS2S,kx)
		self.Args['PoissFact'] = np.asfortranarray(w**-2)
		self.Args['PoissFact_ext'] = np.asfortranarray(w_ext**-2)

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

			self.Args['DepFact'] = np.asfortranarray((2*np.pi)**2/Nx*\
			  np.cos(0.5*np.pi*kr_g/kr_g.max(0).max(0))**2*filt_bandpass*filt_antialias)
			self.Args['DepProj']   = (Rgrid,1./dx,1./dr,kx0)
			self.Args['FBCurrIn']  = (kx_env,InCurr)
			self.Args['FBIn']      = (kx_env,In)
			self.Args['FBout']     = (kx_env,Out)
			self.Args['FBoutFull'] = (kx_env,OutFull)
		else:
			self.Args['DepFact'] = np.asfortranarray( (2*np.pi)**2/Nx*\
			  np.cos(0.5*dr*kr_g)**2*np.cos(0.5*dx*kx_g)[:,:,None]**2)
			self.Args['DepProj']   = (Rgrid,1./dx,1./dr)
			self.Args['FBCurrIn']  = (kx,InCurr)
			self.Args['FBIn']      = (kx,In)
			self.Args['FBout']     = (kx,Out)
			self.Args['FBoutFull'] = (kx,OutFull)

		self.Data={}
		self.Data['EB'] = np.zeros((Nx,Nr,Mtot,6),dtype='complex',order='F')
		self.Data['J'] = np.zeros((Nx,Nr,Mtot,3),dtype='complex',order='F')
		self.Data['scl_spc'] = np.zeros((Nx,Nr,Mtot),dtype='complex',order='F')

		self.Data['EG_fb'] = np.zeros((Nx,Nkr,Mtot,6),dtype='complex',order='F')
		self.Data['vec_fb'] = np.zeros((Nx,Nkr,Mtot,3),dtype='complex',order='F')
		self.Data['scl_fb'] = np.zeros((Nx,Nkr,Mtot),dtype='complex',order='F')

		self.Data['J_fb']  = np.zeros_like(self.Data['vec_fb'])
		self.Data['B_fb']  = np.zeros((Nx,Nkr,Mtot_ext,3),dtype='complex',order='F')

		if 'SpaceCharge' in self.Configs['Features'] or 'StaticKick' in self.Configs['Features']:
			print 'Space charge is added'
			self.Data['Rho']        = np.zeros_like(self.Data['scl_spc'])
			self.Data['Rho_fb']     = np.zeros_like(self.Data['scl_fb'])
			self.Data['gradRho_fb_nxt'] = np.zeros_like(self.Data['vec_fb'])
			self.Data['gradRho_fb_prv'] = np.zeros_like(self.Data['vec_fb'])

		if 'StillAsBackground' in self.Configs['Features']:
			print '"Still" species will be treated as background'
			self.Data['BckGrndRho'] = np.zeros_like(self.Data['scl_spc'])

		if 'CoPropagative' in self.Configs:
			self.PSATD_coeffs(self.Configs['CoPropagative'])
		else:
			self.PSATD_coeffs()
		if 'NoPoissonCorrection' in self.Configs['Features']:	print 'Poisson correction will not be performed'

	def PSATD_coeffs(self,beta=1):
		kx_g,w = self.Args['kx_g'],self.Args['w']
		kx_g = beta*kx_g[:,:,None]
		dt = self.Configs['TimeStep']
		if 'SpaceCharge' in self.Configs['Features']:
			self.Data['PSATD_E'] = np.zeros(w.shape+(5,),dtype='double',order='F')
		elif 'KxShift' in self.Configs:
			self.Data['PSATD_E'] = np.zeros(w.shape+(3,),dtype='complex',order='F')
		else:
			self.Data['PSATD_E'] = np.zeros(w.shape+(3,),dtype='double',order='F')

		self.Data['PSATD_G'] = np.zeros_like(self.Data['PSATD_E'])

		self.Data['PSATD_E'][:,:,:,0] = np.cos(dt*w)
		self.Data['PSATD_E'][:,:,:,1] = np.sin(dt*w)/w
		self.Data['PSATD_G'][:,:,:,0] = -w*np.sin(dt*w)
		self.Data['PSATD_G'][:,:,:,1] = np.cos(dt*w)

		if 'KxShift' in self.Configs:
			self.Data['PSATD_E'][:,:,:,2] = np.exp(-0.5j*kx_g*dt)*1j*kx_g/(w**2-kx_g**2)*\
			  (1.0-np.exp(1j*kx_g*dt)*(w/(1j*kx_g)*np.sin(dt*w)+np.cos(dt*w)))
			self.Data['PSATD_G'][:,:,:,2]  = np.exp(-0.5j*kx_g*dt)*w**2/(w**2-kx_g**2)*\
			  (1.0+np.exp(1j*kx_g*dt)*(1j*kx_g/w*np.sin(dt*w)-np.cos(dt*w)))
		else:
			self.Data['PSATD_E'][:,:,:,2] = -np.sin(dt*w)/w
			self.Data['PSATD_G'][:,:,:,2] = 1-np.cos(dt*w)

		if 'SpaceCharge' in self.Configs['Features']:
			self.Data['PSATD_E'][:,:,:,3] = (dt*w*np.cos(dt*w)-np.sin(dt*w))/w**3/dt
			self.Data['PSATD_E'][:,:,:,4] = (np.sin(dt*w)-dt*w)/w**3/dt
			self.Data['PSATD_G'][:,:,:,3] = (1-np.cos(dt*w)-dt*w*np.sin(dt*w))/w**2/dt
			self.Data['PSATD_G'][:,:,:,4] = (np.cos(dt*w)-1)/w**2/dt

	def maxwell_solver(self):
		if 'SpaceCharge' in self.Configs['Features']:
			self.Data['EG_fb'] = chimera.maxwell_push_with_spchrg(self.Data['EG_fb'],self.Data['J_fb'],self.Data['gradRho_fb_prv'],\
			  self.Data['gradRho_fb_nxt'],self.Data['PSATD_E'],self.Data['PSATD_G'])
		else:
			self.Data['EG_fb'] = chimera.maxwell_push_wo_spchrg(self.Data['EG_fb'],self.Data['J_fb'],self.Data['PSATD_E'],self.Data['PSATD_G'])

	def poiss_corr(self):
		if 'NoPoissonCorrection' in self.Configs['Features']: return
		self.Data['vec_fb'][:] = self.Data['J_fb']
		self.FBGradDiv()
		if 'SpaceCharge' in self.Configs['Features']:
			self.Data['J_fb'] = chimera.poiss_corr_with_spchrg(self.Data['J_fb'],self.Data['vec_fb'],\
			  self.Data['gradRho_fb_prv'], self.Data['gradRho_fb_nxt'],self.Args['PoissFact'],1./self.Configs['TimeStep'])
		else:
			self.Data['J_fb'] = chimera.poiss_corr_wo_spchrg(self.Data['J_fb'],self.Data['vec_fb'],self.Args['PoissFact'])

	def maxwell_solver_stat(self,px0):
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
		self.Data['EG_fb'] = chimera.maxwell_init_push(self.Data['EG_fb'],self.Data['J_fb'],self.Data['gradRho_fb_nxt'],CPSATD1,CPSATD2)

	def poiss_corr_stat(self,px0):
		if 'NoPoissonCorrection' in self.Configs['Features']: return
		beta0 = px0/np.sqrt(1+px0**2)
		kx_g = beta0*self.Args['kx_g'][:,:,None,None]
		self.Data['vec_fb'][:] = self.Data['J_fb']
		self.FBGradDiv()
		self.Data['J_fb'] = chimera.poiss_corr_with_spchrg(self.Data['J_fb'],self.Data['vec_fb'],\
		  np.zeros_like(self.Data['gradRho_fb_nxt']),-1.j*kx_g*self.Data['gradRho_fb_nxt'],self.Args['PoissFact'],1.)

	def field_drift(self,px0):
		beta0 = px0/np.sqrt(1+px0**2)
		kx_g = beta0*self.Args['kx_g'][:,:,None,None]
		self.Data['EG_fb'] *= np.exp(-0.5j*kx_g*self.Configs['TimeStep'])

	def fb_curr_in(self):
		self.Data['J_fb'] = chimera.fb_vec_in(self.Data['J_fb'],self.Data['J'],self.Args['leftX'],*self.Args['FBCurrIn'])
		self.Data['J_fb'] = chimera.omp_mult_vec(self.Data['J_fb'],self.Args['DepFact'])

	def fb_dens_in(self):
		self.Data['Rho_fb'] = chimera.fb_scl_in(self.Data['Rho_fb'],self.Data['Rho'],self.Args['leftX'],*self.Args['FBCurrIn'])
		self.Data['Rho_fb'] = chimera.omp_mult_scl(self.Data['Rho_fb'],self.Args['DepFact'])
		self.FBGradDens()

	def fb_scl_spc_in(self):
		self.Data['scl_fb'] = chimera.fb_scl_in(self.Data['scl_fb'],self.Data['scl_spc'],self.Args['leftX'],*self.Args['FBIn'])

	def fb_fld_out(self):
		self.Data['EB'] = chimera.fb_eb_out(self.Data['EB'],self.Data['EG_fb'][:,:,:,:3],self.Data['B_fb'],\
		  self.Args['leftX'],*self.Args['FBout'])
		if 'KxShift' in self.Configs:
			self.Data['EB'] = chimera.eb_corr_axis_env(self.Data['EB'])
		else:
			self.Data['EB'] = chimera.eb_corr_axis(self.Data['EB'])

	def FBGrad(self):
		if 'KxShift' in self.Configs:
			self.Data['vec_fb'] = chimera.fb_grad_env(self.Data['vec_fb'],self.Data['scl_fb'],*self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_grad(self.Data['vec_fb'],self.Data['scl_fb'],*self.Args['FBDiff'])

	def FBDiv(self):
		if 'KxShift' in self.Configs:
			self.Data['scl_fb'] = chimera.fb_div_env(self.Data['scl_fb'],self.Data['vec_fb'],*self.Args['FBDiff'])
		else:
			self.Data['scl_fb'] = chimera.fb_div(self.Data['scl_fb'],self.Data['vec_fb'],*self.Args['FBDiff'])

	def FBGradDiv(self):
		if 'KxShift' in self.Configs:
			self.Data['vec_fb'] = chimera.fb_graddiv_env(self.Data['vec_fb'],*self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_graddiv_ext(self.Data['vec_fb'],*self.Args['FBDiff'])

	def FBDivGrad(self):
		self.Data['scl_fb'] = chimera.fb_divgrad(self.Data['scl_fb'],*self.Args['FBDiff'])

	def FBGradDens(self):
		if 'KxShift' in self.Configs:
			self.Data['gradRho_fb_nxt'] = chimera.fb_grad_env(self.Data['gradRho_fb_nxt'],self.Data['Rho_fb'],*self.Args['FBDiff'])
		else:
			self.Data['gradRho_fb_nxt'] = chimera.fb_grad_ext(self.Data['gradRho_fb_nxt'],self.Data['Rho_fb'],*self.Args['FBDiff'])

	def G2B_FBRot(self):
		if 'KxShift' in self.Configs:
			self.Data['B_fb'] = chimera.fb_rot_env(self.Data['B_fb'], self.Data['EG_fb'][:,:,:,3:],*self.Args['FBDiff'])
		else:
			self.Data['B_fb'] = chimera.fb_rot_plus(self.Data['B_fb'], self.Data['EG_fb'][:,:,:,3:],*self.Args['FBDiff'])
		self.Data['B_fb'] = chimera.omp_mult_vec(self.Data['B_fb'], self.Args['PoissFact_ext'])

	def B2G_FBRot(self):
		if 'KxShift' in self.Configs:
			self.Data['EG_fb'][:,:,:,3:] = chimera.fb_rot_env(self.Data['EG_fb'][:,:,:,3:],self.Data['B_fb'],*self.Args['FBDiff'])
		else:
			self.Data['EG_fb'][:,:,:,3:] = chimera.fb_rot_minus(self.Data['EG_fb'][:,:,:,3:],self.Data['B_fb'],*self.Args['FBDiff'])

	def add_gauss_beam(self,S):
		k0 = 2*np.pi*S['k0']
		a0 = 2*np.pi*S['a0']
		X_focus = S['x0']-S['x_foc']
		kx_g, kr_g,w = self.Args['kx_g'],self.Args['kr_g'],self.Args['w']
		w = w[:,:,:,None]

		if 'KxShift' in self.Configs:
			e_s0 = a0*0.5*(np.pi)**0.5*S['Lx']*S['LR']**2*self.Args['dkx']/self.Args['lengthR']**2
			self.Data['vec_fb'][:,:,self.Args['Nko'],2] = e_s0/j1(self.Args['lengthR']*kr_g[:,:,self.Args['Nko']])**2*\
			  np.exp(-1j*kx_g*S['x0'])*\
			  np.exp(-0.25*(kx_g-k0)**2*S['Lx']**2 -0.25*kr_g[:,:,self.Args['Nko']]**2*S['LR']**2 )
			DT = -1.j*w
		else:
			Xgrid,Rgrid = self.Args['Xgrid'],self.Args['Rgrid']	# cos laser phase
			self.Data['scl_spc'][:,:,0] = a0*np.cos(k0*(Xgrid[:,None]-S['x0']))*\
			  np.exp(-(Xgrid[:,None]-S['x0'])**2/S['Lx']**2-Rgrid[None,:]**2/S['LR']**2) *\
			  (abs(Xgrid[:,None]-S['x0'])< 3.0*S['Lx'])*(abs(Rgrid[None,:])< 3.0*S['LR'])
			self.Data['scl_spc'][:,0,0] = 0.0
			self.fb_scl_spc_in()
			self.Data['vec_fb'][:,:,:,2] = self.Data['scl_fb']/np.float(self.Args['Nx'])
			DT = -1.j*w*np.sign(kx_g[:,:,None,None] + (kx_g[:,:,None,None]==0))

		EE = self.Data['vec_fb'].copy()
		self.FBGradDiv()
		self.Data['vec_fb'] = chimera.omp_mult_vec(self.Data['vec_fb'], self.Args['PoissFact'])
		EE += self.Data['vec_fb']
		GG  = DT*EE

		self.Data['vec_fb'][:] = np.cos(w*X_focus)*EE + np.sin(w*X_focus)/w*GG
		GG = -w*np.sin(w*X_focus)*EE + np.cos(w*X_focus)*GG
		EE[:] = self.Data['vec_fb']
		EE *= np.exp(1.j*kx_g[:,:,None,None]*X_focus)
		GG *= np.exp(1.j*kx_g[:,:,None,None]*X_focus)

		self.Data['EG_fb'][:,:,:,:3] += EE
		self.Data['EG_fb'][:,:,:,3:] += GG
		self.Data['vec_fb'][:] = 0.0
		self.Data['scl_fb'][:] = 0.0

	def absorb_field(self,Lf,config='left'):
		Nfilt = int(Lf/self.Args['dx'])
		flt_gr = np.arange(Nfilt)
		filt_shape = (flt_gr>=0.75*Nfilt)*0.25*(1-np.cos(np.pi*(flt_gr-0.75*Nfilt)/(0.75*Nfilt)))**2
		filt = np.ones(self.Args['Nx'])
		if config=='left':
			filt[:Nfilt] = filt_shape
		elif config=='right':
			filt[-Nfilt:] = filt_shape[::-1]
		elif config=='both':
			filt[:Nfilt] = filt_shape
			filt[-Nfilt:] = filt_shape[::-1]
		self.Data['EG_fb'][:,:,:,:3] = chimera.fb_filtr(self.Data['EG_fb'][:,:,:,:3],self.Args['leftX'],\
		  self.Args['kx'],filt)
		self.Data['B_fb'][:] = chimera.fb_filtr(self.Data['B_fb'],self.Args['leftX'],\
		  self.Args['kx'],filt)
		self.B2G_FBRot()

	def FBRot(self):
		if 'KxShift' in self.Configs:
			self.Data['vec_fb'] = chimera.fb_rot_env(np.empty_like(self.Data['vec_fb']),self.Data['vec_fb'],*self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_rot(np.empty_like(self.Data['vec_fb']),self.Data['vec_fb'],*self.Args['FBDiff'])
