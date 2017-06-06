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

from __future__ import print_function,division
import numpy as np
from numpy.linalg import inv as inv
from scipy.special import jn_zeros,jn,j1
from scipy.constants import m_e, c, e, epsilon_0
import chimera.moduls.fimera as chimera
from chimera.moduls.messages import msg

class Solver:
	def __init__(self,solver_in,msg=msg):
		"""
		Solver constructor for CHIMERA

		Parameters
		----------
		solver_in: dictionary
		  defines all solver parameters and features

		See Also
		--------
		examples of solver configuratins are given in ./doc/

		"""

		# Export input as a main arguments container
		self.Args = solver_in

		# Set up the basic constants
		if 'Features' not in self.Args:
			self.Args['Features'] = ()
		if 'TimeActive' not in self.Args:
			self.Args['TimeInterval'] = (0.0,np.inf)

		self.Args['Nko'] = self.Args['MaxAzimuthMode']
		self.Args['leftX'],self.Args['rightX'], self.Args['lengthR'], \
		  self.Args['dx'],self.Args['dr'] = self.Args['Grid']
		self.Args['dt'] = self.Args['TimeStep']

		self.Args['dx_inv'] = 1.0/self.Args['dx']
		self.Args['dr_inv'] = 1.0/self.Args['dr']
		self.Args['dt_inv'] = 1.0/self.Args['dt']

		# Define x-grid with nodes number divisible by 2 and number of chunks
		if 'Xchunked' in self.Args:
			self.Args['nthrds'] = self.Args['Xchunked'][0]
			self.Args['Nxchunk'] = 2*int(np.round(\
			  0.5/self.Args['nthrds']/self.Args['dx'] \
			  * (self.Args['rightX'] - self.Args['leftX'])))
			self.Args['Nx'] = self.Args['Nxchunk']*self.Args['nthrds']
		else:
			self.Args['Nx'] = int(2*np.round(0.5/self.Args['dx'] \
			  *(self.Args['rightX'] - self.Args['leftX'])))

		self.Args['Xgrid']  = self.Args['rightX'] \
		   - self.Args['dx']*np.arange(self.Args['Nx'])[::-1]
		self.Args['leftX'] = self.Args['Xgrid'][0]

		#  Define spectral x-grid and its parameters
		if 'KxShift' in self.Args:
			self.Args['kx0'] = 2*np.pi*self.Args['KxShift']
		else:
			self.Args['kx0'] = 0.0

		self.Args['kx_env'] = 2 * np.pi \
		  * np.fft.fftfreq(self.Args['Nx'],self.Args['dx'])
		self.Args['kx']     = self.Args['kx0'] + self.Args['kx_env']
		self.Args['dkx'] = (self.Args['kx'][1]-self.Args['kx'][0]) / (2.*np.pi)

		# Define r-grid with dr/2 offset
		self.Args['Nkr'] = int(np.round(self.Args['lengthR']/self.Args['dr']))
		self.Args['Nr'] = self.Args['Nkr'] + 1
		self.Args['RgridFull'] = self.Args['dr'] \
		                       * (np.arange(self.Args['Nr']) - 0.5)
		self.Args['lengthR'] = self.Args['RgridFull'][-1] + self.Args['dr']

		# Create the spectral r-grid and all DHT operators
		self.Args['InFull'], self.Args['InPlusFull'], self.Args['InMinusFull'], \
		  self.Args['OutFull'], self.Args['OutPlusFull'], \
		  self.Args['OutMinusFull'], self.Args['w'], self.Args['kr'], \
		  self.Args['kr_g'], self.Args['kx_g'], self.Args['idxM'] \
		    = self.get_spectral_operators()

		self.Args['Mtot'] = self.Args['InFull'].shape[-1]
		self.Args['PoissFact'] = 1. / self.Args['w']**2

		self.Args['VGrid'] = 2 * np.pi * self.Args['dx'] * self.Args['dr'] \
		  * self.Args['RgridFull']
		self.Args['VGrid'] = (self.Args['RgridFull']>0.0) / self.Args['VGrid']

		self.Args['InCurr'] = self.Args['InFull'] \
		  * self.Args['VGrid'][None,1:,None]
		self.Args['InCurrPlus'] = self.Args['InFullPlus'] \
		  * self.Args['VGrid'][None,1:,None]
		self.Args['InCurrMinus'] = self.Args['InFullMinus'] \
		  * self.Args['VGrid'][None,1:,None]

		# Print the reports
		msg.print_('sol/bounds',\
		  [self.Args[kw] for kw in ('leftX','rightX','lengthR')] )
		msg.print_('sol/resols',\
		  [self.Args[kw] for kw in ('dx','dr','dt')] )
		msg.print_('sol/actv',self.Args['TimeInterval'] )
		msg.print_('sol/grdsz',\
		  [self.Args[kw] for kw in ('Nx','Nkr','Mtot')] )
		if 'KxShift' in self.Args: msg.print_('sol/shift')

		# Truncate r-grid and corresponding operators
		if 'Rcut' in self.Args:
			indRcut = (self.Args['RgridFull']<self.Args['Rcut']).sum()

			for kw in ('Rgrid','Out', 'OutPlus', 'OutMinus',):
				self.Args[kw] = self.Args[kw+'Full'][:indRcut-1]

			for kw in ('In', 'InPlus', 'InMinus',):
				self.Args[kw] = self.Args[kw+'Full'][:,:indRcut-1]

			for kw in ('InCurr', 'InCurrPlus', 'InCurrMinus',):
				self.Args[kw] = self.Args[kw][:,:indRcut-1]

			self.Args['Rcut'] = self.Args['Rgrid'].max()
			msg.print_('sol/rcut', (self.Args['Rcut'],indRcut) )
		else:
			for kw in ('Rgrid','In', 'InPlus', 'InMinus', \
			           'Out', 'OutPlus', 'OutMinus',):
				self.Args[kw] = self.Args[kw+'Full']

		self.Args['Nr'] = self.Args['Rgrid'].shape[0]
		self.Args['lowerR'] = (self.Args['Rgrid'] * (self.Args['Rgrid']>=0)).min()
		self.Args['upperR'] = self.Args['Rgrid'].max()

		# make factor for energy integration
		tmp_idx = (np.abs(self.Args['idxM']) + 1)[None,None,:]
		self.Args['EnergyFact'] = \
		  0.5e-2 * (m_e*c**2/e)**2 * (4*np.pi*epsilon_0) \
		  * self.Args['lengthR']**2 / self.Args['dkx'] \
		  * jn(tmp_idx, self.Args['kr_g'] * self.Args['lengthR'])**2

		# Reshape spectral operators and fix all arrays to Fortran order
		for kw in ('In', 'InPlus', 'InMinus', \
		           'Out', 'OutPlus', 'OutMinus', \
			        'InCurr', 'InCurrPlus', 'InCurrMinus'):
			self.Args[kw] = np.swapaxes(self.Args[kw],0,1)
			if (kw+'Full') in self.Args.keys():
				self.Args[kw+'Full'] = np.swapaxes(self.Args[kw+'Full'],0,1)

		for kw in self.Args.keys():
			if type(self.Args[kw]) ==  np.ndarray:
				self.Args[kw] = np.asfortranarray(self.Args[kw])

		# Prepare the arguments packs, generate low-pass filter
		self.Args['DepProj'] = [self.Args[kw] \
		                        for kw in ('Rgrid','dx_inv','dr_inv')]

		if 'KxShift' in self.Args:
			filt_x = self.antialias()
			self.Args['DepProj'] += (self.Args['kx0'],)
			kx_base = [self.Args['kx_env'],]
		else:
			kx_max = np.abs(self.Args['kx_g']).max()
			filt_x = np.cos(np.pi/2*self.Args['kx_g']/kx_max)[:,:,None]**2
			kx_base = [self.Args['kx'],]

		kr_max = self.Args['kr_g'].max(0).max(0)
		self.Args['DepFact'] = (2*np.pi)**2  / self.Args['Nx'] * filt_x \
		                       * np.cos(np.pi/2*self.Args['kr_g']/kr_max)**2
		self.Args['DepFact'] =  np.asfortranarray(self.Args['DepFact'])

		self.Args['FBCurrIn']  = kx_base \
		     + [self.Args[kw] for kw ('InCurr', 'InCurrPlus', 'InCurrMinus')]

		self.Args['FBIn']  = kx_base \
		     + [self.Args[kw] for kw ('In', 'InPlus', 'InMinus')]

		self.Args['FBInFull']  = kx_base \
		     + [self.Args[kw] for kw ('InFull', 'InPlusFull', 'InMinusFull')]

		self.Args['FBout']     = kx_base \
		     + [self.Args[kw] for kw ('Out', 'OutPlus', 'OutMinus')]

		self.Args['FBoutFull'] = kx_base \
		     + [self.Args[kw] for kw ('OutFull', 'OutPlusFull', 'OutMinusFull')]

		# Define necessary data stuctures
		vectors_spc = ('E', 'B', 'J')
		vectors_fb = ('E_fb', 'B_fb', 'J_fb', 'vec_fb', 'G_fb')
		scalars_spc = ('scl_spc', )
		scalars_fb = ('scl_fb', )

		# ** Add data structures for charge density
		if ('SpaceCharge' in self.Args['Features']) and \
		   ('StaticKick' in self.Args['Features']):
			msg.print_('sol/spcchrg')

			scalars_spc+= ('Rho',)
			scalars_fb += ('Rho_fb',)
			vectors_fb += ('gradRho_fb_nxt','gradRho_fb_prv')
		# ** Add data structure for background charge density
		if 'StillAsBackground' in self.Args['Features']:
			msg.print_('sol/stillspcs')
			scalars_spc+= ('BckGrndRho',)

		# Define shapes of data stuctures
		baseSP_shape = tuple([self.Args[kw] for kw in ('Nx','Nr','Mtot')])
		baseFB_shape = tuple([self.Args[kw] for kw in ('Nx','Nkr','Mtot')])
		base_confg = {'dtype':'complex', 'order':'F'}

		# Allocate all data structures
		self.Data={}

		for kw in vectors_spc:
			self.Data[kw] = np.zeros(baseSP_shape + (3,), **base_confg)
		for kw in vectors_fb:
			self.Data[kw] = np.zeros(baseFB_shape + (3,), **base_confg)
		for kw in scalars_spc:
			self.Data[kw] = np.zeros(baseSP_shape, **base_confg)
		for kw in scalars_fb:
			self.Data[kw] = np.zeros(baseFB_shape, **base_confg)

		# Make coefficients for Maxwell equations
		if 'CoPropagative' in self.Args:
			self.PSATD_coeffs(self.Args['CoPropagative'])
		else:
			self.PSATD_coeffs()

		# Alarm if NoPoissonCorrection is activated
		if 'NoPoissonCorrection' in self.Args['Features']:
			msg.print_('sol/pois')

		# saving the log
		self.log = '\n'.join(msg.log)

	def PSATD_coeffs(self,beta=1.):
		"""
		Make the PSATD coefficients for Maxwell equations

		Parameters
		----------
		beta: float
		  reference frame velocity

		"""

		# Copy grid axis and scale (in case of a moving frame)
		kx_g,w = self.Args['kx_g'], self.Args['w']
		kx_g = beta*kx_g[:,:,None]
		dt = self.Args['TimeStep']

		# Prepare the data structures for coefficients
		if 'SpaceCharge' in self.Args['Features']:
			baseFB_shape = w.shape + (5,)
			base_confg = {'dtype':'double', 'order':'F'}
		elif 'KxShift' in self.Args:
			baseFB_shape = w.shape + (3,)
			base_confg = {'dtype':'complex', 'order':'F'}
		else:
			baseFB_shape = w.shape + (3,)
			base_confg = {'dtype':'double', 'order':'F'}

		for kw in ('PSATD_E','PSATD_G'):
			self.Data[kw] = np.zeros(baseFB_shape,**base_confg)

		# Add generic solver-independents coefficients
		self.Data['PSATD_E'][:,:,:,0] = np.cos(dt*w)
		self.Data['PSATD_E'][:,:,:,1] = np.sin(dt*w) / w
		self.Data['PSATD_G'][:,:,:,0] = -w * np.sin(dt*w)
		self.Data['PSATD_G'][:,:,:,1] = np.cos(dt*w)

		# Add solver-dependents coefficients
		if 'KxShift' in self.Args:
			self.Data['PSATD_E'][:,:,:,2] = 1j * kx_g * np.exp(-0.5j*kx_g*dt)  \
			     / (w**2-kx_g**2) * (1. -  np.exp(1j*kx_g*dt)                  \
			                          *(np.cos(dt*w) - 1j*w/kx_g*np.sin(dt*w)) \
			                         )
			self.Data['PSATD_G'][:,:,:,2]  = np.exp(-0.5j*kx_g*dt) * w**2    \
			     / (w**2-kx_g**2) * (1. + np.exp(1j*kx_g*dt)                 \
			                          *(1j*kx_g/w*np.sin(dt*w)-np.cos(dt*w)) \
			                         )
		else:
			self.Data['PSATD_E'][:,:,:,2] = -np.sin(dt*w) / w
			self.Data['PSATD_G'][:,:,:,2] = 1 - np.cos(dt*w)

		# ** add coefficients for charge density gradients (for full PIC)
		if 'SpaceCharge' in self.Args['Features']:
			self.Data['PSATD_E'][:,:,:,3] = (dt*w*np.cos(dt*w)-np.sin(dt*w)) \
			                                / w**3 / dt
			self.Data['PSATD_E'][:,:,:,4] = (np.sin(dt*w)-dt*w)/w**3/dt
			self.Data['PSATD_G'][:,:,:,3] = (1-np.cos(dt*w)-dt*w*np.sin(dt*w)) \
			                                 / w**2 / dt
			self.Data['PSATD_G'][:,:,:,4] = (np.cos(dt*w)-1) / w**2 / dt

	def maxwell_solver(self):
		"""
		Update the fields over one time-step

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'SpaceCharge' in self.Args['Features']:
			self.Data['E_fb'], self.Data['G_fb'] = \
			  chimera.maxwell_push_with_spchrg(self.Data['E_fb'],           \
			                                   self.Data['G_fb'],           \
			                                   self.Data['J_fb'],           \
			                                   self.Data['gradRho_fb_prv'], \
			                                   self.Data['gradRho_fb_nxt'], \
			                                   self.Data['PSATD_E'],        \
			                                   self.Data['PSATD_G'])
		else:
			self.Data['E_fb'] self.Data['G_fb'] = \
			  chimera.maxwell_push_wo_spchrg(self.Data['E_fb'],    \
			                                 self.Data['G_fb'],    \
			                                 self.Data['J_fb'],    \
			                                 self.Data['PSATD_E'], \
			                                 self.Data['PSATD_G'])

	def poiss_corr(self,poiss_corr_num = 2):
		"""
		Add the current Poisson correction

		Parameters
		----------
		poiss_corr_num: integer
		  number of iterations to do

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'NoPoissonCorrection' in self.Args['Features']:
			return

		for corr in range(poiss_corr_num):
			self.Data['vec_fb'][:] = self.Data['J_fb']
			self.FBGradDiv()
			if 'SpaceCharge' in self.Args['Features']:
				self.Data['J_fb'] = chimera.poiss_corr(           \
				                     self.Data['J_fb'],           \
				                     self.Data['vec_fb'],         \
				                     self.Data['gradRho_fb_prv'], \
				                     self.Data['gradRho_fb_nxt'], \
				                     self.Args['dt_inv'],         \
				                     self.Args['PoissFact'])
			else:
				self.Data['vec_fb'] = chimera.omp_mult_vec(self.Data['vec_fb'], \
				                                           self.Args['PoissFact'])
				self.Data['J_fb'] = chimera.omp_add_vec(self.Data['J_fb'], \
				                                        self.Data['vec_fb'])

	def maxwell_solver_stat(self,px0):
		"""
		Update the fields over one time-step, considering static solution

		Parameters
		----------
		px0 : float
		  average momentum of the species

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if ('SpaceCharge' not in self.Args['Features']) and \
		   ('StaticKick' not in self.Args['Features']):
			return

		beta0 = px0/np.sqrt(1+px0**2)
		kx_g, w = self.Args['kx_g'], self.Args['w']
		CPSATD1 = np.zeros(w.shape+(2,),dtype='complex',order='F')
		CPSATD2 = np.zeros(w.shape+(2,),dtype='complex',order='F')

		kx_g = beta0*kx_g[:,:,None]
		CPSATD1[:,:,:,0] = 1.j * kx_g / (w**2-kx_g**2)
		CPSATD1[:,:,:,1] = -1. / (w**2-kx_g**2)
		CPSATD2[:,:,:,0] = w**2 / (w**2-kx_g**2)
		CPSATD2[:,:,:,1] = 1.j * kx_g / (w**2-kx_g**2)
		self.Data['EG_fb'], self.Data['EG_fb'] = chimera.maxwell_init_push(    \
		                                          self.Data['E_fb'],           \
		                                          self.Data['G_fb'],           \
		                                          self.Data['J_fb'],           \
		                                          self.Data['gradRho_fb_nxt'], \
		                                          CPSATD1,                     \
		                                          CPSATD2)

	def poiss_corr_stat(self,px0):
		"""
		Add the current Poisson correction

		Parameters
		----------
		px0 : float
		  average momentum of the species

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if ('NoPoissonCorrection' in self.Args['Features']):
			return

		DT = -1.j * px0 / np.sqrt(1+px0**2) * self.Args['kx']
		self.Data['vec_fb'][:] = self.Data['J_fb']
		self.FBGradDiv()

		self.Data['J_fb'] = chimera.poiss_corr_stat(self.Data['J_fb'],          \
		                                            self.Data['vec_fb'],        \
		                                            self.Data['gradRho_fb_nxt'],\
		                                            DT, self.Args['PoissFact'])

	def field_drift(self,px0):
		"""
		Propagate the fields in free space over a timestep

		Parameters
		----------
		px0 : float
		  average momentum of the species

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		self.Data['EG_fb'] = chimera.field_drift(self.Data['EG_fb'],      \
		                                         self.Args['kx'],         \
		                                         px0 / np.sqrt(1+px0**2), \
		                                         self.Args['TimeStep'])

	def fb_curr_in(self):
		"""
		Project the current into the Fourier-Bessel space

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		self.Data['J_fb'] = chimera.fb_vec_in(self.Data['J_fb'],
		                                      self.Data['J'],         \
		                                      self.Args['leftX'],     \
		                                      *self.Args['FBCurrIn'])
		self.Data['J_fb'] = chimera.omp_mult_vec(self.Data['J_fb'],   \
		                                         self.Args['DepFact'])

	def fb_dens_in(self):
		"""
		Project the charge density into the Fourier-Bessel space

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		self.Data['Rho_fb'] = chimera.fb_scl_in(self.Data['Rho_fb'],   \
		                                        self.Data['Rho'],      \
		                                        self.Args['leftX'],    \
		                                        *self.Args['FBCurrIn'])

		self.Data['Rho_fb'] = chimera.omp_mult_scl(self.Data['Rho_fb'],\
		                                           self.Args['DepFact'])

	def fb_scl_spc_in(self):
		"""
		Project a scalar into the Fourier-Bessel space

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		self.Data['scl_fb'] = chimera.fb_scl_in(self.Data['scl_fb'],  \
		                                        self.Data['scl_spc'], \
		                                        self.Args['leftX'],   \
		                                        *self.Args['FBIn'])

	def fb_fld_out(self):
		"""
		Project the fields out from the Fourier-Bessel space

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		self.Data['E'], self.Data['B'] = chimera.fb_eb_out(self.Data['E'],     \
		                                                   self.Data['B'],     \
		                                                   self.Data['E_fb'],  \
		                                                   self.Data['B_fb'],  \
		                                                   self.Data['B_fb'],  \
		                                                   self.Args['leftX'], \
		                                                   *self.Args['FBout'])
		if 'KxShift' in self.Args:
			self.Data['EB'] = chimera.eb_corr_axis_env(self.Data['EB'])
		else:
			self.Data['EB'] = chimera.eb_corr_axis(self.Data['EB'])

	def FBGrad(self):
		"""
		Calculate the gradient of 'scl_fb' and write it into 'vec_fb'

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'KxShift' in self.Args:
			self.Data['vec_fb'] = chimera.fb_grad_env(self.Data['vec_fb'], \
			                                          self.Data['scl_fb'], \
		                                             *self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_grad(self.Data['vec_fb'], \
			                                      self.Data['scl_fb'], \
			                                      *self.Args['FBDiff'])

	def FBDiv(self):
		"""
		Calculate the divergence of 'vec_fb' and write it into 'scl_fb'

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'KxShift' in self.Args:
			self.Data['scl_fb'] = chimera.fb_div_env(self.Data['scl_fb'], \
			                                         self.Data['vec_fb'], \
			                                         *self.Args['FBDiff'])
		else:
			self.Data['scl_fb'] = chimera.fb_div(self.Data['scl_fb'], \
			                                     self.Data['vec_fb'], \
			                                     *self.Args['FBDiff'])

	def FBGradDiv(self):
		"""
		Calculate the gardient of divergence of 'vec_fb' and
		write it into 'vec_fb'

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""
		if 'KxShift' in self.Args:
			self.Data['vec_fb'] = chimera.fb_graddiv_env(\
			  self.Data['vec_fb'],*self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_graddiv(\
			  self.Data['vec_fb'],*self.Args['FBDiff'])

	def FBGradDens(self):
		"""
		Calculate the gradient of 'Rho_fb' and write it into 'gradRho_fb_nxt'

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'KxShift' in self.Args:
			self.Data['gradRho_fb_nxt'] = chimera.fb_grad_env(          \
			                               self.Data['gradRho_fb_nxt'], \
			                               self.Data['Rho_fb'],         \
			                               *self.Args['FBDiff'])
		else:
			self.Data['gradRho_fb_nxt'] = chimera.fb_grad(\
			                               self.Data['gradRho_fb_nxt'], \
			                               self.Data['Rho_fb'],         \
			                               *self.Args['FBDiff'])

	def G2B_FBRot(self):
		"""
		Calculate the magnetic field from G=rot(B);

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'KxShift' in self.Args:
			self.Data['B_fb'] = chimera.fb_rot_env(            \
			                     self.Data['B_fb'],            \
			                     self.Data['EG_fb'][:,:,:,3:], \
			                     *self.Args['FBDiff'])
		else:
			self.Data['B_fb'] = chimera.fb_rot(                \
			                     self.Data['B_fb'],            \
			                     self.Data['EG_fb'][:,:,:,3:], \
			                     *self.Args['FBDiff'])

		self.Data['B_fb'] = chimera.omp_mult_vec(self.Data['B_fb'], \
		                                         self.Args['PoissFact'])

	def add_gauss_beam(self,S):
		"""
		Add Gaussian pulse to the solver

		Parameters
		----------
		S: dict
		  Parameters of the pulse

		Comments
		--------
		examples of laser beam configuratins are given in ./doc/

		"""

		k0 = 2*np.pi*S['k0']
		a0 = 2*np.pi*S['a0']
		X_focus = S['x0']-S['x_foc']
		kx_g, kr_g,w = self.Args['kx_g'],self.Args['kr_g'],self.Args['w']
		w = w[:,:,:,None]

		if 'KxShift' in self.Args:
			e_s0 = a0 * 0.5 * (np.pi)**0.5 * S['Lx'] * S['LR']**2 \
			       * self.Args['dkx'] / self.Args['lengthR']**2

			self.Data['vec_fb'][:,:,self.Args['Nko'],2] =                        \
			   e_s0 / j1(self.Args['lengthR']*kr_g[:,:,self.Args['Nko']])**2     \
			   * np.exp(-1j*kx_g*S['x0']) * np.exp(-0.25*(kx_g-k0)**2*S['Lx']**2 \
			   -0.25*kr_g[:,:,self.Args['Nko']]**2*S['LR']**2 )

			DT = -1.j * w
		else:
			Xgrid,Rgrid = self.Args['Xgrid'],self.Args['Rgrid']	# sin phase
			self.Data['scl_spc'][:,:,0] = a0 * np.sin(k0*(Xgrid[:,None]-S['x0']))\
			  * np.exp(-(Xgrid[:,None]-S['x0'])**2/S['Lx']**2 \
			  - Rgrid[None,:]**2/S['LR']**2) * (abs(Rgrid[None,:])< 3.5*S['LR']) \
			  * (abs(Xgrid[:,None]-S['x0'])< 3.5*S['Lx'])
			self.Data['scl_spc'][:,0,0] = 0.0
			self.fb_scl_spc_in()
			self.Data['vec_fb'][:,:,:,2] =  self.Data['scl_fb']/self.Args['Nx']

			DT = -1.j*w*np.sign(kx_g[:,:,None,None] + (kx_g[:,:,None,None]==0))

		EE = self.Data['vec_fb'].copy()
		EE = self.div_clean(EE)
		GG  = DT*EE

		self.Data['vec_fb'][:] = np.cos(w*X_focus)*EE + np.sin(w*X_focus)/w*GG
		GG = -w*np.sin(w*X_focus)*EE + np.cos(w*X_focus)*GG
		EE[:] = self.Data['vec_fb']
		EE *= np.exp(1.j * kx_g[:,:,None,None] * X_focus)
		GG *= np.exp(1.j * kx_g[:,:,None,None] * X_focus)

		self.Data['EG_fb'][:,:,:,:3] += EE
		self.Data['EG_fb'][:,:,:,3:] += GG

		self.Data['vec_fb'][:] = 0.0
		self.Data['scl_fb'][:] = 0.0

	def get_damp_profile(self,Lf):
		Nfilt = int(Lf)
		flt_gr = np.arange(Nfilt)
		filt_shape = (flt_gr>=0.75*Nfilt) \
		  * (0.5 - 0.5*np.cos(np.pi*(flt_gr-0.75*Nfilt)/(0.25*Nfilt)))**2
		return filt_shape

	def damp_field(self,config='left'):
		mode = {'left':0,'right':1,'both':2}
		self.Data['EG_fb'][:,:,:,:3] = chimera.fb_filtr(\
		  self.Data['EG_fb'][:,:,:,:3],self.Args['leftX'],self.Args['kx'],\
		  self.Args['damp_profile'], mode[config])
		self.Data['EG_fb'][:,:,:,3:] = chimera.fb_filtr(\
		  self.Data['EG_fb'][:,:,:,3:],self.Args['leftX'],self.Args['kx'],\
		  self.Args['damp_profile'], mode[config])

	def FBRot(self):
		if 'KxShift' in self.Args:
			self.Data['vec_fb'] = chimera.fb_rot_env(\
			  np.empty_like(self.Data['vec_fb']),self.Data['vec_fb'],\
			  *self.Args['FBDiff'])
		else:
			self.Data['vec_fb'] = chimera.fb_rot(\
			  np.empty_like(self.Data['vec_fb']),self.Data['vec_fb'],\
			  *self.Args['FBDiff'])

	def div_clean(self,vec):
		self.Data['vec_fb'][:] = vec
		self.FBGradDiv()
		self.Data['vec_fb'] = chimera.omp_mult_vec(self.Data['vec_fb'],\
		  self.Args['PoissFact'])
		vec = chimera.omp_add_vec(vec,self.Data['vec_fb'])
		return vec

	def B2G_FBRot(self):
		if 'KxShift' in self.Args:
			self.Data['EG_fb'][:,:,:,3:] = chimera.fb_rot_env(\
			  self.Data['EG_fb'][:,:,:,3:],self.Data['B_fb'],*self.Args['FBDiff'])
		else:
			self.Data['EG_fb'][:,:,:,3:] = chimera.fb_rot(\
			  self.Data['EG_fb'][:,:,:,3:],self.Data['B_fb'],*self.Args['FBDiff'])

	def antialias(self):
		kx_env, kx, kx0, Nx = \
		  [self.Args[key] for key in ['kx_env', 'kx', 'kx0','Nx']]
		cutafter = 0.85
		fu_bandpass = lambda x : (x<cutafter) + (x>=cutafter) \
		  * np.cos(np.pi/2*(x-cutafter)/(1-cutafter))**2
		filt_bandpass = fu_bandpass(np.abs(kx_env) \
		  / np.abs(kx_env.max()))[:,None,None]
		filt_antialias = np.ones_like(filt_bandpass)

		if 'NoAntiEcho' not in self.Args['Features']:
			msg.print_('sol/echo_intro')

			if 'AntiEchoStrength' in self.Args['Features']:
				ae = self.Args['Features']['AntiEchoStrength']
			else:
				ae = 2

			num_echoes = np.int(np.abs(kx).max()/np.abs(kx_env).max()) + 1
			cell_echos = np.abs(kx_env).max()/kx0*np.arange(num_echoes) - 1.
			full_band = np.array([kx.min()/kx0, kx.max()/kx0]) - 1.

			fu_antialias = lambda x, x0, ae0 :\
			  1 - np.exp(-(x-x0)**2/(ae0*(x0+1.)/Nx)**2)

			echo_ind = 0
			echo_order = 0
			for cellecho in cell_echos:
				echo_order+=1
				if cellecho>full_band[0] and cellecho<full_band[1]:
					ech = cellecho/abs(full_band).max()

					echo_str = '**** {0:d}-th grid echo at {1:.5g} '\
					.format(echo_order, ech)

					if type(ae)==list or type(ae)==tuple:
						ae_loc = ae[echo_ind]
					else:
						ae_loc = ae
					echo_ind +=1
					if abs(cellecho)/abs(full_band).max()<0.4:
						echo_str += '(close to resonance) '

					if ae_loc <= 0:
						echo_str += 'no correction'
						print( echo_str)
						continue

					filt_antialias *= fu_antialias(\
					  kx/kx0 - 1, cellecho, ae_loc)[:,None,None]

					echo_str += 'correcting with strength {0:g}'.\
					format(ae_loc)
					print( echo_str)
		filt_tot = filt_antialias*filt_bandpass
		return filt_tot

	def get_spectral_operators(self, ext=0):
		Nx, Nkr, Nko, lengthR, kx, RgridFull = [self.Args[kw] \
		  for kw in ('Nx','Nkr','Nko','lengthR','kx','RgridFull')]

		if 'KxShift' in self.Args:
			Mmin, Mmax, Mtot = -Nko, Nko, 2*Nko + 1
		else:
			Mmin, Mmax, Mtot = 0, Nko, Nko + 1

		# Allocate the data in Args
		self.Args['kr'] = np.zeros((Nkr, Mtot))
		for kw in ('kr_g', 'w'):
			self.Args[kw] = np.zeros((Nx, Nkr, Mtot))
		for kw in ('OutFull','OutPlusFull', 'OutMinusFull', \
		        'InFull','InPlusFull', 'InMinusFull'):
			self.Args[kw] = np.zeros((Nkr, Nkr, Mtot))

		# Calculate the DHT axis and operators
		self.Args['m'] = np.arange(Mmin,Mmax+1)

		for jm in self.Args['m']:
			self.Args['kr'][:,jm] = jn_zeros(jm,Nkr) / lengthR

			self.Args['kr_g'][:,:,jm], self.Args['kx_g'] = \
			       np.meshgrid(self.Args['kr'][:,jm], kx)

			self.Args['w'][:,:,jm] = np.sqrt(self.Args['kx_g']**2 \
			                    + self.Args['kr_g'][:,:,jm]**2)

			self.Args['OutFull'][:,:,jm] = jn(jm, \
			                 RgridFull[1:,None] * self.Args['kr'][:,jm][None,:])
			self.Args['OutPlusFull'][:,:,jm] = jn(jm+1, \
			                 RgridFull[1:,None]*self.Args['kr'][:,jm][None,:])
			self.Args['OutMinusFull'][:,:,jm] = jn(jm-1, \
			                 RgridFull[1:,None]*self.Args['kr'][:,jm][None,:])

			self.Args['InFull'][:,:,jm] = inv(self.Args['OutFull'][:,:,jm])
			self.Args['InPlusFull'][:,:,jm] = inv(self.Args['OutPlusFull'][:,:,jm])
			self.Args['InMinusFull'][:,:,jm]=inv(self.Args['OutMinusFull'][:,:,jm])

		# Reorder modes for Fortran indexing (shifted mode only)
		if ('KxShift' in self.Args) and (Nko>0):
			for kw in ('InFull', 'InPlusFull', 'InMinusFull', \
			           'OutFull', 'OutPlusFull', 'OutMinusFull', 'w', 'kr_g'):
				self.Args[kw] = np.concatenate(\
			     (self.Args[kw][:,:,Mmin:], self.Args[kw][:,:,:Mmax+1]), axis=-1)

			self.Args['kr'] = np.concatenate(\
			(self.Args['kr'][:,Mmin:], self.Args['kr'][:,:Mmax+1]), axis=-1)

	def divG_clean(self):
		self.Data['vec_fb'][:] = self.Data['EG_fb'][:,:,:,3:]
		self.FBGradDiv()
		self.Data['vec_fb'] = chimera.omp_mult_vec(self.Data['vec_fb'],\
		  self.Args['PoissFact'])
		self.Data['EG_fb'][:,:,:,3:] = chimera.omp_add_vec(\
		  self.Data['EG_fb'][:,:,:,3:],self.Data['vec_fb'])

	def test_calibration(self):
		self.Data['vec_fb'][:] = self.Data['EG_fb'][:,:,:,:3]
		self.FBDiv()
		t1 = (self.Data['scl_fb'] - self.Data['Rho_fb']).copy()
		self.Data['vec_fb'][:] = self.Data['EG_fb'][:,:,:,3:]
		self.FBDiv()
		t2 = self.Data['scl_fb'].copy()
		self.Data['vec_fb'][:] = self.Data['B_fb']
		self.FBDiv()
		t3 = self.Data['scl_fb'].copy()
		return t1,t2,t3

