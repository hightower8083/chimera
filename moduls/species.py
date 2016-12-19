import numpy as np
from inspect import getargspec
import chimera.moduls.fimera as chimera
from scipy.constants import m_e,c,e,epsilon_0

class Specie:
	def __init__(self,PartSpcs):
		self.Configs = PartSpcs
		if 'Features' not in self.Configs: self.Configs['Features'] = ()
		leftX, rightX,lengthR, dx, dr = self.Configs['Grid']

		if 'Xchunked' in self.Configs:
			nthrds = self.Configs['Xchunked'][0]
			Nx = int(np.round(0.5/dx*(rightX - leftX)/nthrds)*2*nthrds)
		else:
			Nx = int(np.round(0.5/dx*(rightX - leftX))*2)

		Nr = int(np.round(lengthR/dr))
		Rgrid = dr*(np.arange(Nr)-0.5)
		Xgrid  = rightX - dx*np.arange(Nx)[::-1]
		leftX = Xgrid[0]

		if 'FixedCell' in self.Configs:
			self.Num_p = np.prod(self.Configs['FixedCell'])
			packX, packR, packO = np.mgrid[\
			  1:self.Configs['FixedCell'][0]:self.Configs['FixedCell'][0]*1j,\
			  1:self.Configs['FixedCell'][1]:self.Configs['FixedCell'][1]*1j,\
			  1:self.Configs['FixedCell'][2]:self.Configs['FixedCell'][2]*1j]
			packX = np.asfortranarray( (packX.ravel()-0.5)/self.Configs['FixedCell'][0])
			packR = np.asfortranarray( (packR.ravel()-0.5)/self.Configs['FixedCell'][1])
			packO = np.asfortranarray( np.exp(2.j*np.pi*(packO.ravel()-1)/self.Configs['FixedCell'][2]) )
			self.Pax = (packX,packR,packO)
		elif 'RandCell' in self.Configs:
			self.Num_p = self.Configs['RandCell']
		else:
			self.Num_p = 0

		if 'Density' in self.Configs:
			self.wght0 = self.Configs['Charge']*self.Configs['Density']*dr*dx*2*np.pi/self.Num_p
		else:
			self.wght0 = 0.0

		self.push_fact = 2*np.pi*self.Configs['Charge']/self.Configs['Mass']
		self.weight2pC = 4*np.pi**2*m_e*c**2*epsilon_0*1e6/e

		self.Args = {'Nx':Nx,'Nr':Nr,'Xgrid':Xgrid,'Rgrid':Rgrid,'leftX':leftX,'rightX':rightX,\
		  'lowerR':(Rgrid*(Rgrid>=0)).min(),'upperR':Rgrid.max(),'dx':dx,'dr':dr,'NpSlice':Nx*self.Num_p}

		if 'MomentaMeans' not in self.Configs:
			self.Configs['MomentaMeans'] = (0.0,0.0,0.0)
		if 'MomentaSpreads' not in self.Configs:
			self.Configs['MomentaSpreads'] = (0.0,0.0,0.0)

		if 'Devices' in self.Configs:
			self.Devices = self.Configs['Devices']
		else:
			self.Devices = ()

		self.Data = {}
		self.Data['EB'] = np.zeros((6,0,),order='F')
		self.Data['coords'] = np.zeros((3,0),order='F')
		self.Data['weights'] = np.zeros((0,),order='F')
		self.Data['coords_halfstep'] = np.zeros_like(self.Data['coords'])
		self.Data['momenta'] = np.zeros_like(self.Data['coords'])

	def gen_parts(self,Domain = None,Xsteps=None,ProfileFunc=None):
		Xgrid = self.Args['Xgrid']
		Rgrid = self.Args['Rgrid']

		if Domain!=None:
			parts_left, parts_right,parts_rad0,parts_rad1 = Domain
			if self.Args['leftX']>parts_right or self.Args['rightX']<parts_left or parts_rad1<self.Args['lowerR'] \
			  or parts_rad0>self.Args['upperR']: return np.zeros((8,0),order='F')

		if Domain!=None:
			ixb,ixe = (Xgrid<parts_left).sum()-1,(Xgrid<parts_right).sum()+1
			Xgrid = Xgrid[ixb:ixe]
			if parts_rad0<=Rgrid.min():
				irb=0
			else:
				irb=(Rgrid<parts_rad0).sum()-1
			if parts_rad1>=Rgrid.max():
				ire = Rgrid.shape[0]
			else:
				ire = (Rgrid<parts_rad1).sum()+1
			Rgrid = Rgrid[irb:ire]
		elif Xsteps!=None:
			Xgrid = Xgrid[-Xsteps:]

		coords = np.zeros((4,Xgrid.shape[0]*Rgrid.shape[0]*self.Num_p),order='F')
		if 'FixedCell' in self.Configs:
			RandPackO = np.random.rand(Xgrid.shape[0],Rgrid.shape[0]).astype('d',order='F')
			coords,Num_loc = chimera.genparts(coords,Xgrid,Rgrid,RandPackO,*self.Pax)
		elif 'RandCell' in self.Configs:
			coords,Num_loc = self.gen_randcell(coords,Xgrid,Rgrid)
		coords = coords[:,:Num_loc]
		Num_loc = coords.shape[1]

		if ProfileFunc == None:
			coords[-1] *= self.wght0
		elif len(getargspec(ProfileFunc).args)==1:
			coords[-1] *= self.wght0*ProfileFunc(coords[0])
		else:
			coords[-1] *= self.wght0*ProfileFunc(*coords[0:3])

		if 'FlatSpectrum' in self.Configs['Features']:
			rand_mom = 2*np.random.rand(3,Num_loc) - 1.0
		else:
			rand_mom = np.random.randn(3,Num_loc)

		px = self.Configs['MomentaMeans'][0] + self.Configs['MomentaSpreads'][0]*rand_mom[0]
		py = self.Configs['MomentaMeans'][1] + self.Configs['MomentaSpreads'][1]*rand_mom[1]
		pz = self.Configs['MomentaMeans'][2] + self.Configs['MomentaSpreads'][2]*rand_mom[2]
		momenta = np.vstack((px,py,pz)).astype('d',order='F')

		weights = coords[-1].astype('d',order='F')
		coords = coords[0:3].astype('d',order='F')
		return [coords,momenta,weights]

	def add_particles(self,coords,momenta,weights):
		Num2actl = self.Data['coords'].shape[-1]
		Num2add = coords.shape[-1]
		self.Data['coords'].resize((3,Num2actl+Num2add), refcheck=False)
		self.Data['coords'][:,-Num2add:] = coords
		self.Data['coords_halfstep'].resize((3,Num2actl+Num2add), refcheck=False)
		self.Data['coords_halfstep'][:,-Num2add:] = coords
		self.Data['momenta'].resize((3,Num2actl+Num2add), refcheck=False)
		self.Data['momenta'][:,-Num2add:] = momenta
		self.Data['weights'].resize((Num2actl+Num2add,), refcheck=False)
		self.Data['weights'][-Num2add:] = weights
 
	def make_field(self):
		if 'Still' in self.Configs['Features']: return
		if self.Data['EB'].shape[-1]!=self.Data['coords'].shape[-1]:
			self.Data['EB'].resize((6,self.Data['coords'].shape[1]), refcheck=False)
#			self.Data['EB'] = np.zeros((6,self.Data['coords'].shape[1]),order='F')
		self.Data['EB'][:] = 0.0

	def make_device(self,i_step=0):
		if 'Still' in self.Configs['Features']: return
		for device in self.Devices:
			pump_fld = device[0]
			self.Data['EB'] = pump_fld(self.Data['coords'],self.Data['EB'],i_step*self.Configs['TimeStep'],*device[1:])

	def push_velocs(self,dt=None):
		if dt==None:dt=self.Configs['TimeStep']
		if self.Data['coords'].shape[-1]==0 or ('Still' in self.Configs['Features']): return
		self.Data['momenta'] = chimera.push_velocs(self.Data['momenta'],self.Data['EB'],self.push_fact*dt)

	def push_coords(self,dt=None):
		dt=self.Configs['TimeStep']
		if self.Data['coords'].shape[1]==0 or ('Still' in self.Configs['Features']): return
		self.Data['coords'],self.Data['coords_halfstep'] = \
		  chimera.push_coords(self.Data['coords'], self.Data['momenta'], self.Data['coords_halfstep'], dt)

	def denoise(self,WaveNums2Kill):
		for k_supp in WaveNums2Kill:
			particles_mirror = self.Data['coords'].copy(order='F')
			particles_mirror[0] = particles_mirror[0] + 0.5/k_supp
			self.Data['coords'] = np.concatenate((self.Data['coords'],particles_mirror),axis=1)

			particles_mirror = self.Data['coords_halfstep'].copy(order='F')
			particles_mirror[0] = particles_mirror[0] + 0.5/k_supp
			self.Data['coords_halfstep'] = np.concatenate((self.Data['coords_halfstep'],particles_mirror),axis=1)

			self.Data['momenta'] = np.concatenate((self.Data['momenta'],self.Data['momenta']),axis=1)
			self.Data['weights'] = np.concatenate((self.Data['weights'],self.Data['weights']),axis=0)
			self.Data['weights'] *= 0.5

	def chunk_coords(self,position=None):
		if 'Xchunked' in self.Configs:
			if self.Data['coords'].shape[-1] == 0: return
			if position=='cntr':
				chnk_ind,self.chunks,outleft,outright  = chimera.chunk_coords(self.Data['coords_halfstep'],\
				  self.Args['Xgrid'],self.Configs['Xchunked'][0])
			else:
				chnk_ind,self.chunks,outleft,outright  = chimera.chunk_coords(self.Data['coords'],\
				  self.Args['Xgrid'],self.Configs['Xchunked'][0])
			if outright == 0:
				chnk_ind = chnk_ind.argsort()[outleft:]
			else:
				chnk_ind = chnk_ind.argsort()[outleft:-outright]
			if outleft!=0 or outright !=0: print('particles out', outleft,outright)

			self.Data['coords'] = chimera.align_data_vec(self.Data['coords'],chnk_ind)
			self.Data['coords_halfstep'] = chimera.align_data_vec(self.Data['coords_halfstep'],chnk_ind)
			self.Data['momenta'] = chimera.align_data_vec(self.Data['momenta'],chnk_ind)
			self.Data['weights'] = chimera.align_data_scl(self.Data['weights'],chnk_ind)

	def damp_particles(self,wind):
		if self.Data['coords'].shape[-1] == 0: return
		SimDom = np.asfortranarray([self.Args['leftX']+wind['AbsorbLayer'],self.Args['rightX'],\
		  0.0, self.Args['upperR']**2])

		if 'Xchunked' in self.Configs and 'NoSorting' not in wind['Features']:
			index2stay,self.chunks,go_out  = chimera.chunk_coords_boundaries(self.Data['coords'],SimDom,\
			  self.Args['Xgrid'],self.Configs['Xchunked'][0])
			index2stay = index2stay.argsort()[go_out:]
			num2stay = index2stay.shape[0]
		else:
			index2stay,num2stay = chimera.sortpartsout(self.Data['coords'],SimDom)
			index2stay = index2stay[:num2stay]

		self.Data['coords'] = chimera.align_data_vec(self.Data['coords'],index2stay)
		self.Data['coords_halfstep'] = chimera.align_data_vec(self.Data['coords_halfstep'],index2stay)
		self.Data['momenta'] = chimera.align_data_vec(self.Data['momenta'],index2stay)
		self.Data['weights'] = chimera.align_data_scl(self.Data['weights'],index2stay)

		self.Data['coords'].resize((3,num2stay), refcheck=False)
		self.Data['coords_halfstep'].resize((3,num2stay), refcheck=False)
		self.Data['momenta'].resize((3,num2stay), refcheck=False)
		self.Data['weights'].resize((num2stay,), refcheck=False)

	def get_dens_on_grid(self,Nko=0):
		VGrid = 2*np.pi*self.Args['dx']*self.Args['dr']*self.Args['Rgrid']
		VGrid = (VGrid+(self.Args['Rgrid']==0))**-1*(self.Args['Rgrid']>0.0)
		dens = np.zeros((self.Args['Nx'],self.Args['Nr'],Nko+1),dtype='complex',order='F')
		dens = chimera.dep_dens(self.Data['coords'],self.Data['weights'],dens,self.Args['leftX'],self.Args['Rgrid'],\
		  1./self.Args['dx'],1/self.Args['dr'])*VGrid[None,:,None]
		return dens

	def gen_randcell(self,coords,Xgrid,Rgrid):
		ip=0
		for ix in np.arange(Xgrid.shape[0]-1):
			for ir in np.arange(Rgrid.shape[0]-1):
				rand_cell = 2*(np.random.rand(3,self.Num_p)-0.5)
				xx_cell = Xgrid[ix] + self.Args['dx']*(np.arange(self.Num_p)+0.5*rand_cell[0])/self.Num_p
				rr_cell = Rgrid[ir] + self.Args['dr']*(np.arange(self.Num_p)+0.5*rand_cell[1])/self.Num_p
				oo_cell = 2*np.pi*(np.arange(self.Num_p)+0.5*rand_cell[2])/self.Num_p
				np.random.shuffle(xx_cell)
				np.random.shuffle(rr_cell)
				np.random.shuffle(oo_cell)
				coords[0,ip:ip+self.Num_p] = xx_cell
				coords[1,ip:ip+self.Num_p] = rr_cell*np.cos(oo_cell)
				coords[2,ip:ip+self.Num_p] = rr_cell*np.sin(oo_cell)
				coords[3,ip:ip+self.Num_p] = rr_cell
				ip += self.Num_p
		return coords,ip

	def beam_focus(self,x_foc):
		gg = (1+self.Data['momenta']**2).sum(0)**0.5
		pzmean = (self.Data['momenta'][0]*self.Data['weights']).sum()/self.Data['weights'].sum()
		self.Data['coords'][0] = self.Data['coords'][0] - \
		  x_foc*(self.Data['momenta'][0]-pzmean)/pzmean/gg**2
		self.Data['coords'][1:3] = self.Data['coords'][1:3] - self.Data['momenta'][1:3]/pzmean*x_foc
		self.Data['coords_halfstep'][:] = self.Data['coords']
