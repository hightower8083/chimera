# This file is a part of CHIMERA software
# CHIMERA is a simulation code for FEL and laser plasma simulations
# Copyright (C)  2016 Igor A. Andriyash <igor.andriyash@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import chimera.moduls.fimera as chimera
from scipy.constants import m_e, c, e, epsilon_0, hbar
from scipy.constants import alpha as alpha_fs
from scipy.interpolate import griddata
from time import localtime

class SR:
	def __init__(self,sr_in):

		self.Args = sr_in
		if 'Features' not in self.Args:
			self.Args['Features'] = ()
		if 'Mode' not in self.Args:
			self.Args['Mode']='far'

		self.chim_norm = 4e-6*np.pi**2*m_e*c**2*epsilon_0/e**2
		self.J_in_um = 2e6*np.pi*hbar*c

		omega_min, omega_max = self.Args['Grid'][0]
		Nom, Nx, Ny = self.Args['Grid'][-1]

		self.Data = {}
		self.Data['Rad'] = np.zeros((Nom,Nx,Ny), order='F')

		if 'WavelengthGrid' in self.Args['Features']:
			self.Args['wavelengths'] = np.r_[1./omega_max:1./omega_min:Nom*1j]
			self.Args['omega'] = 1./self.Args['wavelengths']
		else:
			self.Args['omega'] = np.r_[omega_min:omega_max:Nom*1j]

		if Nom>1:
			self.Args['dw'] = self.Args['omega'][1:] - self.Args['omega'][:-1]
			self.Args['dw'] = np.abs(self.Args['dw'])
		else:
			self.Args['dw'] = np.array([1.,])

		if self.Args['Mode'] == 'far':
			Nth, Nph = Nx, Ny
			theta_min, theta_max  = self.Args['Grid'][1]
			phi_min, phi_max = self.Args['Grid'][2]

			self.Args['theta'] = np.r_[theta_min:theta_max:Nth*1j]
			self.Args['phi'] = phi_min + (phi_max-phi_min)/Nph*np.arange(Nph)

			if Nth>1:
				self.Args['dth'] = self.Args['theta'][1]- self.Args['theta'][0]
			else:
				self.Args['dth'] = 1.

			if Nph>1:
				self.Args['dph'] = self.Args['phi'][1]- self.Args['phi'][0]
			else:
				self.Args['dph'] = 1.

			self.Args['dV'] = self.Args['dw']*self.Args['dth']*self.Args['dph']

			self.Args['DepFact'] = [self.Args['TimeStep'], self.Args['omega'], \
			  np.sin(self.Args['theta']),np.cos(self.Args['theta']), \
			  np.sin(self.Args['phi']),np.cos(self.Args['phi']) ]

		elif self.Args['Mode'] == 'near':
			x_min, x_max  = self.Args['Grid'][1]
			y_min, y_max = self.Args['Grid'][2]
			self.Args['Z'] = self.Args['Grid'][3]

			self.Args['X'] = np.r_[x_min:x_max:Nx*1j]
			self.Args['Y'] = np.r_[y_min:y_max:Ny*1j]

			if Nx>1:
				self.Args['dx'] = self.Args['X'][1]- self.Args['X'][0]
			else:
				self.Args['dx'] = 1.

			if Ny>1:
				self.Args['dy'] = self.Args['Y'][1]- self.Args['Y'][0]
			else:
				self.Args['dy'] = 1.

			self.Args['dV'] = self.Args['dw']*self.Args['dx']*self.Args['dy']

			self.Args['DepFact'] = [self.Args['TimeStep'], self.Args['omega'], \
			  self.Args['X'], self.Args['Y'], self.Args['Z']]

		elif self.Args['Mode'] == 'near-circ':
			Nr, Nph = Nx, Ny

			r_min, r_max  = self.Args['Grid'][1]
			phi_min, phi_max = self.Args['Grid'][2]
			self.Args['Z'] = self.Args['Grid'][3]

			self.Args['R'] = np.r_[r_min:r_max:Nr*1j]
			self.Args['phi'] = phi_min + (phi_max-phi_min)/Nph*np.arange(Nph)

			if Nr>1:
				self.Args['dr'] = self.Args['R'][1]- self.Args['R'][0]
			else:
				self.Args['dr'] = 1.

			if Nph>1:
				self.Args['dph'] = self.Args['phi'][1]- self.Args['phi'][0]
			else:
				self.Args['dph'] = 1.

			self.Args['dV'] = self.Args['dw']*self.Args['dr']*self.Args['dph']

			self.Args['DepFact'] = [self.Args['TimeStep'], self.Args['omega'], \
			                        self.Args['R'], np.sin(self.Args['phi']), \
			                        np.cos(self.Args['phi']), self.Args['Z'] ]

	def init_track(self,Steps,beam):
		Nparts = beam.Data['coords'].shape[-1]
		self.Args['step'] = 0
		self.Data['coords'] = np.zeros((3,Steps,Nparts),order='F')
		self.Data['weights'] = beam.Data['weights'].copy()

		if self.Args['Mode'] == 'far':
			self.Data['momenta_nxt'] = np.zeros((3,Steps,Nparts),order='F')
			self.Data['momenta_prv'] = np.zeros((3,Steps,Nparts),order='F')
			beam.Data['momenta_prv'] = beam.Data['momenta'].copy()
		elif self.Args['Mode'] == 'near' or self.Args['Mode'] == 'near-circ':
			self.Data['momenta'] = np.zeros((3,Steps,Nparts),order='F')

	def add_track(self,beam):
		step = self.Args['step']
		if self.Args['Mode'] == 'far':
			self.Data['coords'][:,step] = beam.Data['coords']
			self.Data['momenta_prv'][:,step] = beam.Data['momenta_prv']
			self.Data['momenta_nxt'][:,step] = beam.Data['momenta']
			beam.Data['momenta_prv'][:] = beam.Data['momenta'][:]
		elif self.Args['Mode'] == 'near' or self.Args['Mode'] == 'near-circ':
			self.Data['coords'][:,step] = beam.Data['coords_halfstep']
			self.Data['momenta'][:,step] = beam.Data['momenta']

		self.Args['step'] += 1

	def damp_track(self,out_folder='./'):
		tt = localtime()
		tm_sgn = (tt.tm_hour,tt.tm_min,tt.tm_sec,tt.tm_mday,tt.tm_mon,tt.tm_year)
		if self.Args['Mode'] == 'far':
			comps = ['coords','momenta_prv','momenta_nxt']
		elif self.Args['Mode'] == 'near' or self.Args['Mode'] == 'near-circ':
			comps = ['coords','momenta']

		for comp in comps:
			fname = out_folder \
			  + 'track_{:}h{:}m{:}s_{:}-{:}-{:}_'+comp+'.npy'.format(*tm_sgn)
			np.save(fname,self.Data[comp])
		self.Args['step'] = 0 ## TBD

	def calculate_spectrum(self,comp='all'):
		if self.Args['Mode'] == 'far':
			self.calculate_spectrum_far(comp=comp)
		elif self.Args['Mode'] == 'near':
			self.calculate_spectrum_near(comp=comp)
		elif self.Args['Mode'] == 'near-circ':
			self.calculate_spectrum_near_circ(comp=comp)

	def calculate_spectrum_far(self,comp='all'):
		if comp != 'all':
			comps = {'x':1, 'y':2, 'z':3}
			self.Data['Rad'] = chimera.sr_calc_far_comp(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta_prv'],\
			  self.Data['momenta_nxt'],\
			  self.Data['weights'], comps[comp], *self.Args['DepFact'])
		else:
			self.Data['Rad'] = chimera.sr_calc_far_tot(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta_prv'],\
			  self.Data['momenta_nxt'],\
			  self.Data['weights'], *self.Args['DepFact'])

	def calculate_spectrum_near(self,comp='all'):
		if comp != 'all':
			comps = {'x':1, 'y':2, 'z':3}
			self.Data['Rad'] = chimera.sr_calc_near_comp(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta'],\
			  self.Data['weights'], comps[comp], *self.Args['DepFact'])
		else:
			self.Data['Rad'] = chimera.sr_calc_near_tot(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta'],\
			  self.Data['weights'], *self.Args['DepFact'])

	def calculate_spectrum_near_circ(self,comp='all'):
		if comp != 'all':
			comps = {'x':1, 'y':2, 'z':3}
			self.Data['Rad'] = chimera.sr_calc_nearcirc_comp(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta'],\
			  self.Data['weights'], comps[comp], *self.Args['DepFact'])
		else:
			self.Data['Rad'] = chimera.sr_calc_nearcirc_tot(\
			  self.Data['Rad'], \
			  self.Data['coords'],self.Data['momenta'],\
			  self.Data['weights'], *self.Args['DepFact'])

	def get_full_spectrum(self, spect_filter=None, chim_units=True, \
	  phot_num=False, lambda0_um = None):
		if self.Args['Mode'] == 'far':
			val = alpha_fs/(4*np.pi**2)*self.Data['Rad']
		elif self.Args['Mode'] == 'near' or self.Args['Mode'] == 'near-circ':
			val = alpha_fs*np.pi/4*self.Data['Rad']

		if spect_filter is not None:
			val *= spect_filter

		if phot_num:
			ax = self.Args['omega']
			val /= ax[:,None,None]
			if chim_units:
				if lambda0_um is None:
					print("Specify normalization wavelength in "+\
					  "microns (lambda0_um) ")
					return np.zeros_like(val)
				val *= self.chim_norm*lambda0_um
		else:
			val *= self.J_in_um
			if chim_units:
				val *= self.chim_norm
			else:
				if lambda0_um is None:
					print("Specify normalization wavelength in "+\
					  "microns (lambda0_um) ")
					return  np.zeros_like(val)
				val /= lambda0_um

		return val

	def get_energy_spectrum(self, spect_filter = None, chim_units=True, \
	  phot_num=False, lambda0_um = None):

		val = self.get_full_spectrum(spect_filter=spect_filter, \
		  chim_units=chim_units, phot_num=phot_num, lambda0_um=lambda0_um)

		if self.Args['Mode'] == 'far':
			val = 0.5*self.Args['dth']*self.Args['dph']*( (val[1:] + val[:-1]) \
			  *np.sin(self.Args['theta'][None,:,None]) ).sum(-1).sum(-1)
		elif self.Args['Mode'] == 'near':
			val = 0.5*self.Args['dx']*self.Args['dy'] \
			      *(val[1:] + val[:-1]).sum(-1).sum(-1)
		elif self.Args['Mode'] == 'near-circ':
			val = 0.5*self.Args['dr']*self.Args['dph']*( (val[1:] + val[:-1]) \
			  *self.Args['R'][None,:,None] ).sum(-1).sum(-1)

		return val

	def get_energy(self, spect_filter=None, chim_units=True, \
	  phot_num=False, lambda0_um = None):

		val = self.get_energy_spectrum(spect_filter=spect_filter, \
		  chim_units=chim_units, phot_num=phot_num, lambda0_um=lambda0_um)

		val = (val*self.Args['dw']).sum()
		return val

	def get_spot(self, k0=None, spect_filter=None, chim_units=True, \
	  phot_num=False, lambda0_um = None):

		val = self.get_full_spectrum(spect_filter=spect_filter, \
		  chim_units=chim_units, phot_num=phot_num, lambda0_um=lambda0_um)

		if k0 is None:
			if val.shape[0]>1:
				val = 0.5*(val[1:] + val[:-1])
			val = (val*self.Args['dw'][:,None,None]).sum(0)
		else:
			ax = self.Args['omega']
			indx = (ax<k0).sum()
			if np.abs(self.Args['omega'][indx+1]-k0) \
			  < np.abs(self.Args['omega'][indx]-k0):
				indx += 1
			val = val[indx]
		return val

	def get_spot_cartesian(self, k0=None, th_part=1.0, bins=(200,200), \
	  spect_filter=None, chim_units=True, \
	  phot_num=False, lambda0_um = None):

		val = self.get_spot(spect_filter=spect_filter, \
		  chim_units=chim_units, k0=k0, phot_num=phot_num, lambda0_um=lambda0_um)

		if self.Args['Mode'] == 'far':
			th,ph = self.Args['theta'], self.Args['phi']
		elif self.Args['Mode'] == 'near-circ':
			th,ph = self.Args['R'], self.Args['phi']
		else:
			print("This function is for 'far' and 'near-circ' modes only")

		ph,th = np.meshgrid(ph,th)
		th_max = th_part*th.max()

		coord = ((th*np.cos(ph)).flatten(), (th*np.sin(ph)).flatten())

		new_coord = np.mgrid[-th_max:th_max:bins[0]*1j,-th_max:th_max:bins[1]*1j]
		val = griddata(coord,val.flatten(),
		    (new_coord[0].flatten(), new_coord[1].flatten()),
		    fill_value=0., method='linear'
		  ).reshape(new_coord[0].shape)
		ext = np.array([-th_max,th_max,-th_max,th_max])
		return val, ext

	def get_spectral_axis(self):
		if 'WavelengthGrid' in self.Args['Features']:
			ax = 0.5*(self.Args['wavelengths'][1:] \
			  + self.Args['wavelengths'][:-1])
		else:
			ax = 0.5*(self.Args['omega'][1:] + self.Args['omega'][:-1])
		return ax
