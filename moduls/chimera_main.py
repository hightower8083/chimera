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
import chimera.moduls.fimera as chimera
import time, os, h5py

class ChimeraRun():
	def __init__(self,SimComps):
		if 'Particles' in SimComps:
			self.Particles = SimComps['Particles']
		else:
			self.Particles = ()
		if 'Solvers' in SimComps:
			self.Solvers = SimComps['Solvers']
		else:
			self.Solvers = ()
		if 'MovingFrames' in SimComps:
			self.MovingFrames = SimComps['MovingFrames']
		else:
			self.MovingFrames = ({},)
		self.init_Moving_Frames()
		self.init_damp_profile()
		self.make_halfstep()

	def init_Moving_Frames(self):
		for wind in  self.MovingFrames:
			if 'Features'   not in wind: wind['Features']   = ()
			if 'TimeActive' not in wind: wind['TimeActive'] = (0.0,np.inf)
			if 'Velocity'   not in wind: wind['Velocity']   = 1.0
			if 'TimeStep'   not in wind: wind['TimeStep']   = 1.0
			if 'Steps'   not in wind: wind['Steps']   = 1
			if 'Staged' in wind['Features']:
				wind['shiftX'] = 0.5*wind['Velocity']*wind['TimeStep']\
				  * wind['Steps']
			else:
				wind['shiftX'] = wind['Velocity']*wind['TimeStep']*wind['Steps']

	def init_damp_profile(self):
		for solver in self.Solvers:
			for wind in self.MovingFrames:
				if 'AbsorbLayer' in wind:
					solver.Args['damp_profile'] = solver.get_damp_profile(\
					  wind['AbsorbLayer'])

	def make_halfstep(self):
		for species in self.Particles:
			if len(self.Solvers)>0:
				args_tmp = self.Solvers[0].Args
			else:
				args_tmp = self.Particles[0].Args
			SimDom = np.asfortranarray([args_tmp['leftX'],args_tmp['rightX'], \
			  0,args_tmp['upperR']**2])
			if SimDom[0]!=SimDom[1]:
				species.chunk_and_damp(SimDom=SimDom)

		self.project_current()
		self.project_density()
		for solver in self.Solvers:
			for species in self.Particles:
				solver.maxwell_solver_stat(species.Args['MomentaMeans'][0])
		self.project_fields()
		for species in self.Particles: species.make_device()
		for species in self.Particles:
			species.push_velocs(dt_frac=0.5)

	def make_step(self,istep):
		self.frame_act(istep)
		for species in self.Particles: species.push_coords()
		self.chunk_particles(istep)
		self.project_current()
		self.frame_act(istep,'stage2')
		self.project_density()
		self.update_fields()
		self.project_fields()
		for species in self.Particles: species.make_device(istep)
		for species in self.Particles: species.push_velocs()

	def project_current(self):
		for solver in self.Solvers:
			if 'NoCurrent' in solver.Args['Features']:
				solver.Data['J'][:] = 0.0
				solver.Data['J_fb'][:] = 0.0
				continue
			self.dep_curr(solver)
			solver.fb_curr_in()

	def project_density(self):
		for solver in self.Solvers:
			if 'SpaceCharge' not in solver.Args['Features'] and \
			   'StaticKick'  not in solver.Args['Features']:
				continue

			if 'StaticKick' not in solver.Args['Features']:
				solver.Data['gradRho_fb_prv'][:] = solver.Data['gradRho_fb_nxt']

			self.dep_dens(solver)
			solver.fb_dens_in()
			solver.FBGradDens()

	def update_fields(self):
		for solver in self.Solvers:
			if 'StaticKick' in solver.Args['Features']:
				solver.Data['EG_fb'][:] = 0.0
				for species in self.Particles:
					PXmean = (species.Data['momenta'][0]\
					  *species.Data['weights']).sum()/species.Data['weights'].sum()
					solver.poiss_corr_stat(PXmean)
					solver.maxwell_solver_stat(PXmean)
					solver.field_drift(PXmean)
			else:
				solver.poiss_corr()
				solver.maxwell_solver()

	def project_fields(self,component='coords'):
		for species in self.Particles: species.make_field()
		for solver in self.Solvers:
			if 'NoPush' in solver.Args['Features']:
				solver.Data['EB'][:] = 0.0
				continue

			solver.G2B_FBRot()
			solver.fb_fld_out()
			for species in self.Particles:
				if species.Data[component].shape[1]==0 \
				  or 'Still' in species.Args['Features']: continue
				if 'KxShift' in solver.Args:
					species.Data['EB'] = chimera.proj_fld_env(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['EB'],species.Data['EB'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					species.Data['EB'] = chimera.proj_fld(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['EB'],species.Data['EB'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])

	def dep_curr(self,solver):
		solver.Data['J'][:] = 0.0
		for species in self.Particles:
			if species.Data['coords'].shape[1] == 0 \
			  or 'Still' in species.Args['Features']: continue
			if 'KxShift' in solver.Args:
				if 'Xchunked' in species.Args:
					solver.Data['J'] = chimera.dep_curr_env_chnk(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],species.chunks,\
					  species.Args['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr_env(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Args:
					solver.Data['J'] = chimera.dep_curr_chnk(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  species.chunks,species.Args['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_dens(self,solver,component='coords'):
		solver.Data['Rho'][:] = 0.0

		if 'StaticKick' in solver.Args['Features']:
			component='coords_halfstep'

		if 'StillAsBackground' in solver.Args['Features']:
			solver.Data['Rho'] += solver.Data['BckGrndRho']

		for species in self.Particles:
			if species.Data[component].shape[1] == 0 \
			  or 'Still' in species.Args['Features']: continue
			if 'KxShift' in solver.Args:
				if 'Xchunked' in species.Args:
					solver.Data['Rho'] = chimera.dep_dens_env_chnk(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],species.chunks,\
					  species.Args['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['Rho'] = chimera.dep_dens_env(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Args:
					solver.Data['Rho'] = chimera.dep_dens_chnk(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],species.chunks,\
					  species.Args['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['Rho'] = chimera.dep_dens(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])

	def dep_bg(self,solver):
		solver.Data['BckGrndRho'][:] = 0.0
		for species in self.Particles:
			if species.Data['coords'].shape[1]==0 \
			  or 'Still' not in species.Args['Features']: continue
			if 'KxShift' in solver.Args:
				if 'Xchunked' in species.Args:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env_chnk(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,\
					  species.Args['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'], solver.Args['leftX'],\
					  *solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Args:
					solver.Data['BckGrndRho'] = chimera.dep_dens_chnk(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,\
					  species.Args['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['BckGrndRho'] = chimera.dep_dens(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])

	def damp_fields(self,wind):
		if 'AbsorbLayer' not in wind: return
		if wind['AbsorbLayer']<=0: return
		for solver in self.Solvers:
			if 'StaticKick' not in solver.Args['Features']: 
				solver.damp_field()

	def damp_plasma(self,wind):
		if 'AbsorbLayer' not in wind: return
		for species in self.Particles: species.chunk_and_damp(wind=wind)

	def add_plasma(self,wind):
		if 'AddPlasma' not in wind: return
		if 'IonsOnTop' in wind['Features']:
			specie1, specie2 = self.Particles[:2]
			parts_to_add = specie1.gen_parts(\
			  Xsteps=int(wind['shiftX']/wind['TimeStep'])+1,\
			  ProfileFunc=wind['AddPlasma'])
			specie1.add_particles(*parts_to_add)
			parts_to_add[-1] *= -1
			specie2.add_particles(*parts_to_add)
		else:
			for species in self.Particles:
				species.add_particles(*species.gen_parts(\
				  Xsteps=int(wind['shiftX']/wind['TimeStep'])+1,\
				  ProfileFunc=wind['AddPlasma']))

	def postframe_corr(self,wind):
		if 'AbsorbLayer' or 'AddPlasma' in wind:
			for solver in self.Solvers:
				if 'StaticKick' in solver.Args['Features']: continue
				if 'SpaceCharge' in solver.Args['Features']:
					if 'StillAsBackground' in solver.Args['Features']:
						self.dep_bg(solver)
					self.dep_dens(solver)

	def move_frame(self,wind):
		for comp in (self.Solvers + self.Particles):
			comp.Args['Xgrid']  += wind['shiftX']
			comp.Args['leftX']  = comp.Args['Xgrid'][ 0]
			comp.Args['rightX'] = comp.Args['Xgrid'][-1]

	def frame_act(self,istep,act='stage1'):
		for wind in self.MovingFrames:
			if istep<wind['TimeActive'][0] or istep>wind['TimeActive'][1]:
				continue
			if np.mod( istep-wind['TimeActive'][0], wind['Steps'])!= 0: continue
			if act=='stage1':
				self.damp_fields(wind)
				self.move_frame(wind)
				self.add_plasma(wind)
				self.damp_plasma(wind)
				self.postframe_corr(wind)
			elif act=='stage2' and 'Staged' in wind['Features']:
				self.move_frame(wind)

	def chunk_particles(self,istep=0):
		for species in self.Particles:
			if  'Xchunked' not in species.Args \
			  or 'NoSorting' in species.Args['Features']: continue
			if np.mod(istep, species.Args['Xchunked'][1]+1)!= 0: continue
			if species.Data['coords'].shape[-1] == 0: continue
			if len(self.Solvers)>0:
				args_tmp = self.Solvers[0].Args
			else:
				args_tmp = self.Particles[0].Args
			SimDom = np.asfortranarray([args_tmp['leftX'],args_tmp['rightX'],\
			  0,args_tmp['Rgrid'].max()**2])
			species.chunk_and_damp(SimDom = SimDom, position='cntr',)



#########################################################################
########## TBD: CHANGE INTO openPMD DIAG ################################

	def drop_snap(self,fname='./snap_',verbose=False):
		fname += time.ctime().replace(' ','_').replace(':','-')+'.hdf5'
		myf = h5py.File(fname, mode='w')
		i=0
		for solver in self.Solvers:
			name = '/solver'+str(i)+'/'
			for key in solver.Data.keys():
				type_loc = type(solver.Data[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Data/'+key)
				myf[name+'Data/'+key] = solver.Data[key]

			for key in solver.Args.keys():
				type_loc = type(solver.Args[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Args/'+key)
				myf[name+'Args/'+key] = solver.Args[key]
			i+=1
		myf['/NumSolvers'] = i

		i=0
		for specie in self.Particles:
			name = '/specie'+str(i)+'/'
			for key in specie.Data.keys():
				type_loc = type(specie.Data[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Data/'+key)
				myf[name+'Data/'+key] = specie.Data[key]
			for key in specie.Args.keys():
				type_loc = type(specie.Args[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Args/'+key)
				myf[name+'Args/'+key] = specie.Args[key]
			i+=1
		myf['/NumSpecies'] = i
		myf['/NumThreads'] = int(os.environ['OMP_NUM_THREADS'])
		myf.close()
		return fname

	def read_snap(self, fname='./test.hdf5',verbose=False):
		myf = h5py.File(fname,mode='r')
		NumSpecies = myf['NumSpecies'].value
		NumSolvers = myf['NumSolvers'].value
		for i in range(NumSolvers):
			solver = self.Solvers[i]
			name = '/solver'+str(i)+'/'
			for key in solver.Data.keys():
				type_loc = type(solver.Data[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Data/'+key)
				solver.Data[key] = myf[name+'Data/'+key].value
			for key in solver.Args.keys():
				type_loc = type(solver.Args[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Args/'+key)
				solver.Args[key] = myf[name+'Args/'+key].value

		for i in range(NumSpecies):
			specie = self.Particles[i]
			name = '/specie'+str(i)+'/'
			for key in specie.Data.keys():
				type_loc = type(specie.Data[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Data/'+key)
				specie.Data[key] = myf[name+'Data/'+key].value
			for key in specie.Args.keys():
				type_loc = type(specie.Args[key])
				if type_loc==list or type_loc==tuple \
				or type_loc==set or type_loc==dict:
					continue
				if verbose: print(name+'Args/'+key)
				specie.Args[key] = myf[name+'Args/'+key].value
