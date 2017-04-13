import numpy as np
import chimera.moduls.fimera as chimera
from scipy.constants import m_e,c, elementary_charge, epsilon_0
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

		if Nom>1: # TO SET UP WITH NON-UNIFORM OMEGA
			self.Args['dw'] = self.Args['omega'][1]- self.Args['omega'][0]
		else:
			self.Args['dw'] = 1.

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

	def add_track(self,beam):
		self.Data['coords'][:,self.Args['step']] = beam.Data['coords']
		self.Data['momenta_prv'][:,self.Args['step']] = beam.Data['momenta_prv']
		self.Data['momenta_nxt'][:,self.Args['step']] = beam.Data['momenta']
		self.Args['step'] += 1

	def damp_track(self,out_folder='./'):
		tt = localtime()
		tm_sgn = (tt.tm_hour,tt.tm_min,tt.tm_sec,tt.tm_mday,tt.tm_mon,tt.tm_year)

		for comp in ['coords','momenta_prv','momenta_nxt']:
			fname = out_folder \
			  + 'track_{:}h{:}m{:}s_{:}-{:}-{:}_'+comp+'.npy'.format(*tm_sgn)
			np.save(fname,self.Data[comp])
		self.Args['step'] = 0 ## TBD

	def get_spectrum_ser(self):
		if self.Configs['Mode'] == 'incoherent':
			self.Data['Rad_incoh'] = chimera.sr_calc_incoh_track_serial(self.Data['Rad_incoh'], \
			  self.Data['coords'],self.Data['momenta_prv'],self.Data['momenta_nxt'],\
			  self.Data['weights'],*self.Args['DepFact'])

	def get_spectrum(self):
		if self.Configs['Mode'] == 'incoherent':
			self.Data['Rad_incoh'] = chimera.sr_calc_incoh_track(self.Data['Rad_incoh'], \
			  self.Data['coords'],self.Data['momenta_prv'],self.Data['momenta_nxt'],\
			  self.Data['weights'],*self.Args['DepFact'])

	def project_current(self,beam,step):
		if self.Configs['Mode'] == 'incoherent':
			self.Data['Rad_incoh'] = chimera.sr_calc_incoh(self.Data['Rad_incoh'], \
			  beam.Data['coords'],beam.Data['momenta_prv'],beam.Data['momenta'],\
			  beam.Data['weights'],self.Configs['Component'],step,*self.Args['DepFact'])

	def project_current_ser(self,beam,step):
		if self.Configs['Mode'] == 'incoherent':
			self.Data['Rad_incoh'] = chimera.sr_calc_incoh_ser(self.Data['Rad_incoh'], \
			  beam.Data['coords'],beam.Data['momenta_prv'],beam.Data['momenta'],\
			  beam.Data['weights'],self.Configs['Component'],step,*self.Args['DepFact'])

	def get_energy(self):
		if self.Configs['Mode'] == 'incoherent':
			nrg = np.abs(self.Data['Rad_incoh'])**2*self.Args['theta'][None,:,None]*self.Args['dV']
			nrg = nrg.sum()
		return nrg
