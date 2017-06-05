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
from scipy.constants import m_e,c,e,epsilon_0
from inspect import getargspec
import chimera.moduls.fimera as chimera

class Specie:
	def __init__(self, specie_in):
		"""
		Particle specie constructor for CHIMERA

		Parameters
		----------
		solver_in: dictionary
		  defines all solver parameters and features

		See Also
		--------
		examples of species configuratins are given in ./doc/

		"""

		# Export input as a main arguments container
		self.Args = specie_in

		# Define default arguments is not specified
		if 'Features' not in self.Args: self.Args['Features'] = ()
		if 'Charge' not in self.Args: self.Args['Charge'] = -1.0
		if 'Mass' not in self.Args: self.Args['Mass'] = 1.0
		if 'Density' not in self.Args: self.Args['Density'] = 0.0
		if 'MomentaMeans' not in self.Args:
			self.Args['MomentaMeans'] = (0.0,0.0,0.0)
		if 'MomentaSpreads' not in self.Args:
			self.Args['MomentaSpreads'] = (0.0,0.0,0.0)
		if 'Grid' not in self.Args:
			self.Args['Grid'] = ( 0.0, 0.0, 0.0, 1.0, 1.0 )

		if 'Devices' in self.Args:
			self.Devices = self.Args['Devices']
		else:
			self.Devices = ()

		# Set up the basic constants
		self.Args['leftX'], self.Args['rightX'], self.Args['lengthR'], \
		  self.Args['dx'], self.Args['dr'] = self.Args['Grid']

		self.Args['push_fact'] = 2*np.pi * self.Args['Charge'] / self.Args['Mass']
		self.Args['weight2pC'] = 4 * np.pi**2 * m_e * c**2 * epsilon_0 * 1e6 / e

		# Define x-grid with nodes number divisible by 2 and number of chunks
		if 'Xchunked' in self.Args:
			self.Args['nthrds'] = int(self.Args['Xchunked'][0])
			self.Args['Nx'] = 2 * self.Args['nthrds'] * int(np.round( \
			  0.5/self.Args['dx'] / self.Args['nthrds'] \
			  * (self.Args['rightX']-self.Args['leftX'])  ))
		else:
			self.Args['Nx'] = 2*int(np.round( 0.5 / self.Args['dx'] * \
			  (self.Args['rightX'] - self.Args['leftX']) ))

		if self.Args['Nx']>0:
			self.Args['Xgrid']  = self.Args['rightX'] - self.Args['dx'] * \
			  np.arange(self.Args['Nx'])[::-1]
		else:
			self.Args['Xgrid']  = np.array([0.,])
		self.Args['leftX'] = self.Args['Xgrid'][0]

		# Define r-grid with dr/2 offset
		self.Args['Nr'] = int(np.round(self.Args['lengthR']/self.Args['dr']))
		if self.Args['Nr']>0:
			self.Args['Rgrid'] = self.Args['dr'] * (np.arange(self.Args['Nr'])-0.5)
		else:
			self.Args['Rgrid'] = np.array([0.,])

		self.Args['lowerR'] = (self.Args['Rgrid'] * (self.Args['Rgrid']>=0)).min()
		self.Args['upperR'] = self.Args['Rgrid'].max()

		# prepare the arguments for particle generators
		if 'FixedCell' in self.Args:
			self.Args['Num_p'] = np.prod(self.Args['FixedCell'])
			packX, packR, packO = np.mgrid[\
			  1:self.Args['FixedCell'][0]:self.Args['FixedCell'][0]*1j,\
			  1:self.Args['FixedCell'][1]:self.Args['FixedCell'][1]*1j,\
			  1:self.Args['FixedCell'][2]:self.Args['FixedCell'][2]*1j]
			packX = np.asfortranarray( \
			  (packX.ravel()-0.5)/self.Args['FixedCell'][0])
			packR = np.asfortranarray( \
			  (packR.ravel()-0.5)/self.Args['FixedCell'][1])
			packO = np.asfortranarray( \
			  np.exp(2.j*np.pi*(packO.ravel()-1)/self.Args['FixedCell'][2]) )
			self.Args['Pax'] = (packX,packR,packO)
		elif 'RandCell' in self.Args:
			self.Args['Num_p'] = self.Args['RandCell']
		else:
			self.Args['Num_p'] = 0

		if self.Args['Num_p'] !=0 :
			self.Args['wght0'] = self.Args['Charge'] * self.Args['Density'] \
			  * self.Args['dr'] * self.Args['dx']* 2 * np.pi / self.Args['Num_p']
		else:
			self.Args['wght0'] = 0.0
		self.Args['NpSlice'] = self.Args['Nx']*self.Args['Num_p']

		# Prepare the data stuctures
		self.Data = {}
		self.Data['EB'] = np.zeros((6,0,),order='F')
		self.Data['coords'] = np.zeros((3,0),order='F')
		self.Data['weights'] = np.zeros((0,),order='F')
		self.Data['coords_halfstep'] = np.zeros_like(self.Data['coords'])
		self.Data['momenta'] = np.zeros_like(self.Data['coords'])

		# Add particle initial coordinates container if needed
		if 'KeepInitPos' in self.Args['Features']:
			self.Data['coords_init'] = np.zeros((3,0),order='F')

	def gen_parts(self,Domain = None,Xsteps=None,ProfileFunc=None):
		"""
		Macro-particle generators

		Parameters
		----------
		Domain: list
		  defines the cylindrical domain to fill with particles
		Xsteps: int
		  defines width of a right-most layer (in cell number) to fill with particles
		ProfileFunc: python function
		  define the density profile along Z (if single-argument),
		  in Z and R (if two arguments) or in X, Y, Z  (if three arguments)

		Returns
		--------
		[coords,momenta,weights] : list of ndarrays with shapes (3,Np), (3,Np), (Np,)
		  particle coordinates, momenta and weights container

		"""

		Xgrid = self.Args['Xgrid']
		Rgrid = self.Args['Rgrid']

		if Domain is not None:
			p_left, p_right,p_low,p_up = Domain
			if self.Args['leftX']>p_right or self.Args['rightX']<p_left\
			  or p_up<self.Args['lowerR'] or p_low>self.Args['upperR']:
				return np.zeros((8,0),order='F')

			ixb,ixe = (Xgrid<p_left).sum() - 1, (Xgrid<p_right).sum() + 1
			Xgrid = Xgrid[ixb:ixe]
			if p_low<=Rgrid.min():
				irb=0
			else:
				irb=(Rgrid<p_low).sum() - 1
			if p_up>=Rgrid.max():
				ire = Rgrid.shape[0]
			else:
				ire = (Rgrid<p_up).sum() + 1
			Rgrid = Rgrid[irb:ire]
		elif Xsteps is not None:
			Xgrid = Xgrid[-2*Xsteps:-Xsteps]

		coords = np.zeros( \
		  (4, Xgrid.shape[0] * Rgrid.shape[0] * self.Args['Num_p']), order='F')
		if 'FixedCell' in self.Args:
			RandPackO = np.random.rand(\
			  Xgrid.shape[0], Rgrid.shape[0]).astype('d',order='F')
			coords, Num_loc = chimera.genparts(\
			  coords, Xgrid, Rgrid, RandPackO,*self.Args['Pax'])
		elif 'RandCell' in self.Args:
			coords,Num_loc = self.gen_randcell(coords,Xgrid,Rgrid)
		coords = coords[:,:Num_loc]
		Num_loc = coords.shape[1]

		if ProfileFunc == None:
			coords[-1] *= self.Args['wght0']
		elif len(getargspec(ProfileFunc).args)==1:
			coords[-1] *= self.Args['wght0']*ProfileFunc(coords[0])
		elif len(getargspec(ProfileFunc).args)==2:
			coords[-1] *= self.Args['wght0'] \
			  *ProfileFunc( coords[0], np.sqrt( (coords[1:3]**2).sum(0) ) )
		elif len(getargspec(ProfileFunc).args)==3:
			coords[-1] *= self.Args['wght0']*ProfileFunc(*coords[0:3])
		else:
			print("can't understand the ProfileFunc")

		if 'FlatSpectrum' in self.Args['Features']:
			rand_mom = 2*np.random.rand(3,Num_loc) - 1.0
		else:
			rand_mom = np.random.randn(3,Num_loc)

		px = self.Args['MomentaMeans'][0] \
		  + self.Args['MomentaSpreads'][0]*rand_mom[0]
		py = self.Args['MomentaMeans'][1] \
		  + self.Args['MomentaSpreads'][1]*rand_mom[1]
		pz = self.Args['MomentaMeans'][2] \
		  + self.Args['MomentaSpreads'][2]*rand_mom[2]
		momenta = np.vstack((px,py,pz)).astype('d',order='F')

		weights = coords[-1].astype('d',order='F')
		coords = coords[0:3].astype('d',order='F')
		return [coords,momenta,weights]

	def add_particles(self,coords,momenta,weights):
		"""
		Add particles produced by gen_parts generators to the specie container

		Parameters
		----------
		coords: ndarray with shape (3, Np)
		  particles coordinates
		momenta: ndarray with shape (3, Np)
		  particles momenta
		weights: ndarray with shape (Np,)
		  particles weights

		"""

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
		if 'KeepInitPos' in self.Args['Features']:
			self.Data['coords_init'].resize((3,Num2actl+Num2add), refcheck=False)
			self.Data['coords_init'][:,-Num2add:] = coords

	def make_field(self):
		"""
		Set up the container for the E and B fields (on particles)

		"""

		if 'Still' in self.Args['Features']: return
		if self.Data['EB'].shape[-1] != self.Data['coords'].shape[-1]:
			self.Data['EB'].resize((6,self.Data['coords'].shape[1]), \
			  refcheck=False)
		self.Data['EB'][:] = 0.0

	def make_device(self,i_step=0):
		"""
		Caclulate E and B fields (on particles) from the external devices

		Parameters
		----------
		i_step: int
		  iteration number for the case of a dynamic device

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		if 'Still' in self.Args['Features']: return
		for device in self.Devices:
			pump_fld = device[0]
			self.Data['EB'] = pump_fld(self.Data['coords'], self.Data['EB'],\
			  i_step*self.Args['TimeStep'], *device[1:])

	def push_velocs(self, dt_frac=1.):
		"""
		Advance particles velocities using "Boris pusher" method

		Parameters
		----------
		dt_frac: float
		  fraction of time step (needed for initial half step)

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		dt = dt_frac*self.Args['TimeStep']
		if self.Data['coords'].shape[-1]==0 \
		  or ('Still' in self.Args['Features']): return
		self.Data['momenta'] = chimera.push_velocs(\
		  self.Data['momenta'],self.Data['EB'],self.Args['push_fact']*dt)

	def push_coords(self,dt_frac=1.):
		"""
		Advance particles coordinates (leap-frog method is used)

		Parameters
		----------
		dt_frac: float
		  fraction of time step (needed for initial half step)

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		dt = dt_frac*self.Args['TimeStep']
		if self.Data['coords'].shape[1]==0 \
		  or ('Still' in self.Args['Features']): return
		self.Data['coords'],self.Data['coords_halfstep'] = chimera.push_coords(\
		  self.Data['coords'], self.Data['momenta'], \
		  self.Data['coords_halfstep'], dt)

	def denoise(self,WaveNums2Kill):
		"""
		Suppress the shot-noise for the given wavenumbers by adding
		the mirror particles in anti-phase (doubles particle numbers each time)

		Parameters
		----------
		WaveNums2Kill: list of floats
		  wavenumbers of which the shot-noise should be suppressed

		"""

		for k_supp in WaveNums2Kill:
			particles_mirror = self.Data['coords'].copy(order='F')
			particles_mirror[0] = particles_mirror[0] + 0.5/k_supp
			self.Data['coords'] = np.concatenate((self.Data['coords'],\
			  particles_mirror),axis=1)

			particles_mirror = self.Data['coords_halfstep'].copy(order='F')
			particles_mirror[0] = particles_mirror[0] + 0.5/k_supp
			self.Data['coords_halfstep'] = np.concatenate((\
			  self.Data['coords_halfstep'],particles_mirror),axis=1)

			self.Data['momenta'] = np.concatenate((\
			  self.Data['momenta'],self.Data['momenta']),axis=1)
			self.Data['weights'] = np.concatenate((\
			  self.Data['weights'],self.Data['weights']),axis=0)
			self.Data['weights'] *= 0.5

	def chunk_and_damp(self, wind=None, SimDom=None, position='stag'):
		"""
		Impose the absorbing boundary conditions for the particles and,
		if parallel computing is enabled, divide them into chunks aligned along Z-axis

		Parameters
		----------
		wind: dict (optional)
		  if routine is called by a moving window, its configurations to be provided
		SimDom: list (optional)
		  the boundaries of the domain may be specified explicitly
		position: string
		  to consider either particles coordinates at t=n*dt (staggered from momenta)
		  or t=(n-1/2)*dt (sychronized with momenta)
		"""

		if self.Data['coords'].shape[-1] == 0: return
		comp = {'cntr':'coords_halfstep','stag':'coords'}
		if wind is None:
			wind = {'Features':(),'AbsorbLayer':0.0}

		if SimDom is None:
			SimDom = np.asfortranarray([\
			  self.Args['leftX'] + wind['AbsorbLayer']*self.Args['dx'],\
			  self.Args['rightX'], 0.0, self.Args['upperR']**2])

		if 'Xchunked' in self.Args and 'NoSorting' not in wind['Features']:
			index2stay,self.chunks,go_out  = chimera.chunk_coords_boundaries(\
			  self.Data[comp[position]],SimDom,self.Args['Xgrid'],\
			  self.Args['Xchunked'][0])
			index2stay = index2stay.argsort()[go_out:]
			num2stay = index2stay.shape[0]
		else:
			index2stay,num2stay = chimera.sortpartsout(self.Data['coords'],SimDom)
			index2stay = index2stay[:num2stay]

		comps = ['coords','coords_halfstep','momenta']
		if 'KeepInitPos' in self.Args['Features']:
			comps += ['coords_init']
		for comp in comps:
			self.Data[comp] = chimera.align_data_vec(\
			  self.Data[comp],index2stay)
			self.Data[comp].resize((3,num2stay), refcheck=False)

		self.Data['weights'] = chimera.align_data_scl(\
		  self.Data['weights'],index2stay)
		self.Data['weights'].resize((num2stay,), refcheck=False)

	def get_dens_on_grid(self,Nko=0):
		"""
		Calculate the species density on the grid (used by diagnostics)

		Parameters
		----------
		Nko: int
		  number of azimuthal modes to be used

		Returns
		----------
		dens: ndarray with dtype=np.complex
		  species density on the grid (angular Fourier decomposition)

		Comments
		--------
		wrapper for the OMP-vectorized Fortran subroutines

		"""

		VGrid = 2*np.pi*self.Args['dx']*self.Args['dr']*self.Args['Rgrid']
		VGrid = (VGrid+(self.Args['Rgrid']==0))**-1*(self.Args['Rgrid']>0.0)
		dens = np.zeros((self.Args['Nx'],self.Args['Nr'],Nko+1),\
		  dtype='complex',order='F')
		dens = chimera.dep_dens(self.Data['coords'],self.Data['weights'],\
		  dens,self.Args['leftX'],self.Args['Rgrid'],1./self.Args['dx'],\
		  1/self.Args['dr'])*VGrid[None,:,None]
		return dens

	def gen_randcell(self,coords,Xgrid,Rgrid):
		"""
		Fill the grid cells with particles using np.random.randn algo

		Parameters
		----------
		coords: ndarray with shape (3,Np)
		  pre-created container for the particles
		Xgrid: ndarray with shape (Nx,)
		  x-axis of the simulations grid
		Rgrid: ndarray with shape (Nr,)
		  r-axis of the simulations grid

		Returns
		----------
		coords: ndarray with shape (3,Np)
		  container with the particles coordinates
		ip: int
		  number of created particles

		"""

		ip=0
		for ix in np.arange(Xgrid.shape[0]-1):
			for ir in np.arange(Rgrid.shape[0]-1):
				rand_cell = 2*(np.random.rand(3,self.Args['Num_p'])-0.5)
				xx_cell = Xgrid[ix] + self.Args['dx'] \
				  * (np.arange(self.Args['Num_p'])+0.5*rand_cell[0]) \
				  / self.Args['Num_p']
				rr_cell = Rgrid[ir] + self.Args['dr']\
				  * (np.arange(self.Args['Num_p'])+0.5*rand_cell[1]) \
				  / self.Args['Num_p']
				oo_cell = 2*np.pi*(np.arange(self.Args['Num_p'])+0.5*rand_cell[2])\
				  / self.Args['Num_p']
				np.random.shuffle(xx_cell)
				np.random.shuffle(rr_cell)
				np.random.shuffle(oo_cell)
				coords[0,ip:ip+self.Args['Num_p']] = xx_cell
				coords[1,ip:ip+self.Args['Num_p']] = rr_cell*np.cos(oo_cell)
				coords[2,ip:ip+self.Args['Num_p']] = rr_cell*np.sin(oo_cell)
				coords[3,ip:ip+self.Args['Num_p']] = rr_cell
				ip += self.Args['Num_p']
		return coords,ip

	def beam_focus(self,x_foc):
		"""
		Add the drift of the species (for the relativistic beams)

		"""

		gg = (1+self.Data['momenta']**2).sum(0)**0.5
		pzmean = (self.Data['momenta'][0]*self.Data['weights']).sum()\
		  / self.Data['weights'].sum()
		self.Data['coords'][0] = self.Data['coords'][0] - \
		  x_foc*(self.Data['momenta'][0]-pzmean)/pzmean/gg**2
		self.Data['coords'][1:3] = self.Data['coords'][1:3] \
		  - self.Data['momenta'][1:3]/pzmean*x_foc
		self.Data['coords_halfstep'][:] = self.Data['coords']
