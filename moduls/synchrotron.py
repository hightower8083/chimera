########## SYNCHROTRON RADIATION CALCULATOR BLOCK #############
out_num = 1             # index of an output (just to mark) ###
species_index = 0       # particles species to do #############
start_ind = 0           # trajectory node to start from  ######
end_ind = None          # trajectory node to end at    ########
d_time = 1              # step over trajectory nodes ##########
d_part = 1              # particles step (reduce number) ######
Nouts = 1               # number of spectrum outputs ##########
out_step = 1            # step for total energy calculation ###
###############################################################


import numpy as np
from sys import path
path.append('./moduls/')
path.append('./src/')
from lib_names import *

# script args: regime ('coh','incoh','both'), component ('x','y','z','all'), number of particles to take

class Calculator:
	def __init__(self,comm,data,args):
		self.data = data
		self.comm = comm
		self.regime, self.comp, self.num_part2take = args[-3:]
		self.num_part2take = int(self.num_part2take)		
		self.out_num = '_'+str(self.data.out_num)
		self.calc = SychRadCalcs[self.regime]
		self.species_index = str(data.species_index)
		self.treat_track()
		self.define_indices()
		self.SPECT = 0.
		self.SPECT_c = 0.
		self.SPECT_i = 0.
		self.Nrg = 0.

	def treat_track(self):
		self.trak = np.load(self.data.out_folder+'trak_'+self.species_index+'.npy',\
		  mmap_mode='r')[self.data.start_ind:self.data.end_ind:self.data.d_time,\
		  :,::self.data.d_part]
		self.num_part2take_full = self.trak.shape[-1]
		self.divisor = np.nan_to_num(self.trak[:,-1,:].sum(axis=-1)).max()
		self.num_part = np.ceil(float(self.num_part2take_full)/self.comm.size)

	def define_indices(self):
		indx = self.comm.rank*self.num_part
		if self.comm.rank == self.comm.size-1:
			indx_last = self.num_part2take_full
		else:
			indx_last = (self.comm.rank+1)*self.num_part

		self.indexes = [indx]
		if (indx_last-indx)>self.num_part2take:
			while True:
				indx = indx + self.num_part2take
				if indx < indx_last:
					self.indexes.append(indx)
				else:
					break
		self.indexes.append(indx_last)

	def main_calc(self):
		for i in np.arange(len(self.indexes)-1):
			vals = [self.trak[:,:,self.indexes[i]:self.indexes[i+1]], comp_ind[self.comp],\
			  self.data.d_time*self.data.dt*self.data.track_step, 0., 2*np.pi, self.data.Nazim,\
	        0.,self.data.theta_max, self.data.Nelev, self.data.omega_min, self.data.omega_max,\
			  self.data.Nfreq, self.data.Nouts, self.data.out_step]
			if self.regime=='both':
				SPECT_c_bite, SPECT_i_bite,Nrg_bite = self.calc(*vals)
				self.SPECT_c = self.SPECT_c + SPECT_c_bite
				self.SPECT_i = self.SPECT_i + SPECT_i_bite
			else:
				SPECT_bite, Nrg_bite = self.calc(*vals)
				self.SPECT = self.SPECT + SPECT_bite
			self.Nrg = self.Nrg + Nrg_bite
		self.comm.barrier()

		if self.regime=='both':
			self.SPECT_c = self.comm.reduce(self.SPECT_c)
			self.SPECT_i = self.comm.reduce(self.SPECT_i)
		else:
			self.SPECT = self.comm.reduce(self.SPECT)
		self.Nrg = self.comm.reduce(self.Nrg)

	def writeout(self,timing):
		if self.comm.rank!=0: return
		if self.regime=='coh':
			spect0 = np.abs(self.SPECT)**2/self.divisor
			np.save(self.data.out_folder+'data_spect3D_coh_'+self.comp+self.out_num+'.npy', spect0)
			np.save(self.data.out_folder+'energy_stepwise_' +self.comp+self.out_num+'.npy',\
			  self.Nrg/self.divisor)
			print('done with %i particles in %f minutes'  ) % (self.num_part2take_full, timing/60.)

		elif self.regime=='both':
			spect_c = np.abs(self.SPECT_c)**2/self.divisor
			spect_i = self.SPECT_i/self.divisor
			np.save(self.data.out_folder+'data_spect3D_coh_'  +self.comp+self.out_num+'.npy', spect_c)
			np.save(self.data.out_folder+'data_spect3D_incoh_'+self.comp+self.out_num+'.npy', spect_i)
			np.save(self.data.out_folder+'energy_stepwise_'   +self.comp+self.out_num+'.npy', \
			  self.Nrg/self.divisor)
			print('done with %i particles in %f minutes'  ) % (self.num_part2take_full, timing/60.)

		elif self.regime=='incoh':
			spect0 = self.SPECT/self.divisor
			np.save(self.data.out_folder+'data_spect3D_incoh_'+self.comp+self.out_num+'.npy', spect0)
			np.save(self.data.out_folder+'energy_stepwise_'   +self.comp+self.out_num+'.npy',\
			  self.Nrg/self.divisor)
			print('done with %i particles in %f minutes'  ) % (self.num_part2take_full, timing/60.)

