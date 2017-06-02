import numpy as np
import chimera.moduls.fimera as chimera
from scipy.constants import m_e,c,elementary_charge,epsilon_0,hbar
from scipy.constants import alpha as alpha_fs
from scipy.interpolate import griddata
from time import localtime

class SR:
	def __init__(self,sr_in):

		self.Configs = sr_in
		if 'Features' not in self.Configs:
			self.Configs['Features'] = ()
		if 'Mode' not in self.Configs:
			self.Configs['Mode']='incoherent'
		if 'Component' not in self.Configs:
			self.Configs['Component']=2

		self.chim_norm = 4e-6*np.pi**2*m_e*c**2*epsilon_0/elementary_charge**2
		self.J_in_um = 2e6*np.pi*hbar*c

		(omega_min,omega_max),(theta_min,theta_max), (phi_min,phi_max),\
		  (Nom,Nth,Nph) = self.Configs['Grid']

		self.Args = {}
		self.Data = {}

		if 'WavelengthGrid' in self.Configs['Features']:
			self.Args['wavelengths'] = np.r_[1./omega_max:1./omega_min:Nom*1j]
			self.Args['omega'] = 1./self.Args['wavelengths']
		else:
			self.Args['omega'] = np.r_[omega_min:omega_max:Nom*1j]
		self.Args['theta'] = np.r_[theta_min:theta_max:Nth*1j]
		self.Args['phi'] = phi_min + (phi_max-phi_min)/Nph*np.arange(Nph)

		if Nom>1:
			self.Args['dw'] = self.Args['omega'][1:] - self.Args['omega'][:-1]
			self.Args['dw'] = np.abs(self.Args['dw'])
		else:
			self.Args['dw'] = np.array([1.,])

		if Nth>1:
			self.Args['dth'] = self.Args['theta'][1]- self.Args['theta'][0]
		else:
			self.Args['dth'] = 1.

		if Nph>1:
			self.Args['dph'] = self.Args['phi'][1]- self.Args['phi'][0]
		else:
			self.Args['dph'] = 1.

		self.Args['dV'] = self.Args['dw']*self.Args['dth']*self.Args['dph']

		self.Args['DepFact'] = [self.Configs['TimeStep'], self.Args['omega'], \
		  np.sin(self.Args['theta']),np.cos(self.Args['theta']), \
		  np.sin(self.Args['phi']),np.cos(self.Args['phi']) ]

		if self.Configs['Mode'] == 'incoherent':
			self.Data['Rad_incoh'] = np.zeros((Nom,Nth,Nph),order='F')
		elif self.Configs['Mode'] == 'coherent':
			self.Data['Rad_coh'] = np.zeros((Nom,Nth,Nph),order='F')
		elif self.Configs['Mode'] == 'all':
			self.Data['Rad_coh'] = np.zeros((Nom,Nth,Nph),order='F')
			self.Data['Rad_incoh'] = np.zeros((Nom,Nth,Nph),order='F')

	def init_track(self,Steps,beam):
		Nparts = beam.Data['coords'].shape[-1]
		self.Data['coords'] = np.zeros((3,Steps,Nparts),order='F')
		self.Data['momenta_prv'] = np.zeros((3,Steps,Nparts),order='F')
		self.Data['momenta_nxt'] = np.zeros((3,Steps,Nparts),order='F')
		self.Data['weights'] = beam.Data['weights'].copy()
		self.Args['step'] = 0
		beam.Data['momenta_prv'] = beam.Data['momenta'].copy()


	def add_track(self,beam):
		self.Data['coords'][:,self.Args['step']] = beam.Data['coords']
		self.Data['momenta_prv'][:,self.Args['step']] = beam.Data['momenta_prv']
		self.Data['momenta_nxt'][:,self.Args['step']] = beam.Data['momenta']
		self.Args['step'] += 1
		beam.Data['momenta_prv'][:] = beam.Data['momenta'][:]

	def damp_track(self,out_folder='./'):
		tt = localtime()
		tm_sgn = (tt.tm_hour,tt.tm_min,tt.tm_sec,tt.tm_mday,tt.tm_mon,tt.tm_year)

		for comp in ['coords','momenta_prv','momenta_nxt']:
			fname = out_folder \
			  + 'track_{:}h{:}m{:}s_{:}-{:}-{:}_'+comp+'.npy'.format(*tm_sgn)
			np.save(fname,self.Data[comp])
		self.Args['step'] = 0 ## TBD

	def calculate_spectrum(self,comp='all'):
		if self.Configs['Mode'] == 'incoherent':
			if comp != 'all':
				comps = {'x':0, 'y':1, 'z':2}
				self.Data['Rad_incoh'] = chimera.sr_calc_incoh_comp(\
				  self.Data['Rad_incoh'], \
				  self.Data['coords'],self.Data['momenta_prv'],\
				  self.Data['momenta_nxt'],\
				  self.Data['weights'], comps[comp], *self.Args['DepFact'])
			else:
				self.Data['Rad_incoh'] = chimera.sr_calc_incoh_tot(\
				  self.Data['Rad_incoh'], \
				  self.Data['coords'],self.Data['momenta_prv'],\
				  self.Data['momenta_nxt'],\
				  self.Data['weights'], *self.Args['DepFact'])			

	def get_full_spectrum(self, spect_filter=None, chim_units=True, \
	  phot_num=False):
		if self.Configs['Mode'] == 'incoherent':
			val = alpha_fs/(4*np.pi**2)*self.Data['Rad_incoh']
			if spect_filter is not None:
				val *= spect_filter
			if chim_units:
				val *= self.chim_norm
			if phot_num:
				ax = self.Args['omega']
				val /= ax[:,None,None]
			return val

	def get_energy_spectrum(self, spect_filter = None, chim_units=True, \
	  phot_num=False):
		if self.Configs['Mode'] == 'incoherent':
			val = self.get_full_spectrum( \
			  spect_filter=spect_filter, chim_units=chim_units,phot_num=phot_num)
			val = 0.5*self.Args['dth']*self.Args['dph']*( (val[1:] + val[:-1]) \
			  *np.sin(self.Args['theta'][None,:,None]) ).sum(-1).sum(-1)
			return val

	def get_spot(self,spect_filter=None, chim_units=True, k0=None, \
	  phot_num=False):
		if self.Configs['Mode'] == 'incoherent':
			val = self.get_full_spectrum(\
			  spect_filter=spect_filter, chim_units=chim_units,phot_num=phot_num)
			if k0 is None:
				if val.shape[0]>1:
					val = 0.5*(val[1:] + val[:-1])					
				val = self.J_in_um*(val*self.Args['dw'][:,None,None]).sum(0)
			else:
				ax = self.Args['omega']
				indx = (ax<k0).sum()
				if np.abs(self.Args['omega'][indx+1]-k0) \
				  < np.abs(self.Args['omega'][indx]-k0):
					indx += 1
				val = self.J_in_um*val[indx]
			return val

	def get_spot_cartesian(self, th_part=1.0, bins=(200,200), \
	  spect_filter=None, chim_units=True, k0=None, phot_num=False):

		val = self.get_spot(spect_filter=spect_filter, \
		  chim_units=chim_units, k0=k0, phot_num=phot_num)

		th,ph = self.Args['theta'], self.Args['phi']
		ph,th = np.meshgrid(ph,th)
		coord = ((np.sin(th)*np.cos(ph)).flatten(),\
		  (np.sin(th)*np.sin(ph)).flatten())
		th_max = th_part*th.max()
		new_coord = np.mgrid[-th_max:th_max:bins[0]*1j,-th_max:th_max:bins[1]*1j]
		val = griddata(coord,val.flatten(),
		    (new_coord[0].flatten(), new_coord[1].flatten()),
		    fill_value=0., method='linear'
		  ).reshape(new_coord[0].shape)
		ext = np.array([-th_max,th_max,-th_max,th_max])
		return val, ext
		
	def get_energy(self,spect_filter=None, chim_units=True):
		if self.Configs['Mode'] == 'incoherent':
			val = self.get_energy_spectrum( \
			  spect_filter=spect_filter, chim_units=chim_units)
			val *= self.J_in_um
			val = (val*self.Args['dw']).sum()
		return val

	def get_spectral_axis(self):
		if self.Configs['Mode'] == 'incoherent':
			if 'WavelengthGrid' in self.Configs['Features']:
				ax = 0.5*(self.Args['wavelengths'][1:] \
				  + self.Args['wavelengths'][:-1])
			else:
				ax = 0.5*(self.Args['omega'][1:] + self.Args['omega'][:-1])
			return ax
