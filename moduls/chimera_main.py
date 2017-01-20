import numpy as np
from sys import path
path.append('./moduls/')
import fimera as chimera

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
			self.init_Moving_Frames()
		else:
			self.MovingFrames = ()
		self.init_damp_profile()
		self.make_halfstep()


	def init_Moving_Frames(self):
		for wind in  self.MovingFrames:
			if 'Features'   not in wind: wind['Features']   = ()
			if 'TimeActive' not in wind: wind['TimeActive'] = (0.0,np.inf)
			if 'Velocity'   not in wind: wind['Velocity']   = 1.0
			if 'Staged' in wind['Features']:
				wind['shiftX'] = 0.5*wind['Velocity']*wind['TimeStep']\
				  * wind['Steps']
			else:
				wind['shiftX'] = wind['Velocity']*wind['TimeStep']*wind['Steps']

	def init_damp_profile(self):
		for solver in self.Solvers:
			solver.Args['damp_profile'] = np.ones(solver.Args['Nx'])
			for wind in self.MovingFrames:
				if 'AbsorbLayer' not in wind: continue
				if wind['AbsorbLayer']<=0: continue
				solver.Args['damp_profile'] *= solver.get_damp_profile(\
				  wind['AbsorbLayer'])

	def make_halfstep(self):
		for species in self.Particles:
			species.chunk_coords()
		self.project_current()
		self.project_density()
		for solver in self.Solvers:
			for species in self.Particles:
				solver.maxwell_solver_stat(species.Configs['MomentaMeans'][0])
			solver.G2B_FBRot()
		self.project_fields()
		for species in self.Particles: species.make_device()
		for species in self.Particles: 
			species.push_velocs(dt=0.5*species.Configs['TimeStep'])

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
			self.dep_curr(solver)
			solver.fb_curr_in()

	def project_density(self):
		for solver in self.Solvers:
			if 'SpaceCharge' not in solver.Configs['Features'] and \
			  'StaticKick' not in solver.Configs['Features']: continue
			if 'StaticKick' not in solver.Configs['Features']:
				solver.Data['gradRho_fb_prv'][:] = solver.Data['gradRho_fb_nxt']
			self.dep_dens(solver)
			solver.fb_dens_in()
			solver.FBGradDens()

	def update_fields(self):
		for solver in self.Solvers:
			if 'StaticKick' in solver.Configs['Features']:
				solver.Data['EG_fb'][:] = 0.0
				for species in self.Particles:
					PXmean = (species.Data['momenta'][0]\
					  *species.Data['weights']).sum()/species.Data['weights'].sum()
					solver.poiss_corr_stat(PXmean)
					solver.maxwell_solver_stat(PXmean)
					solver.field_drift(PXmean)
				solver.G2B_FBRot()
			else:
				solver.poiss_corr()
				solver.maxwell_solver()
				solver.G2B_FBRot()

	def project_fields(self,component='coords'):
		for species in self.Particles: species.make_field()
		for solver in self.Solvers:
			solver.fb_fld_out()
			for species in self.Particles:
				if species.Data[component].shape[1]==0 \
				  or 'Still' in species.Configs['Features']: continue
				if 'KxShift' in solver.Configs:
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
			  or 'Still' in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.Data['J'] = chimera.dep_curr_env_chnk(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],species.chunks,\
					  species.Configs['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr_env(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.Data['J'] = chimera.dep_curr_chnk(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  species.chunks,species.Configs['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr(\
					  species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],\
					  solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_dens(self,solver,component='coords'):
		solver.Data['Rho'][:] = 0.0
		if 'StaticKick' in solver.Configs['Features']: 
			component='coords_halfstep'
		if 'StillAsBackground' in solver.Configs['Features']: 
			solver.Data['Rho'] += solver.Data['BckGrndRho']
		for species in self.Particles:
			if species.Data[component].shape[1] == 0 \
			  or 'Still' in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.Data['Rho'] = chimera.dep_dens_env_chnk(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],species.chunks,\
					  species.Configs['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['Rho'] = chimera.dep_dens_env(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.Data['Rho'] = chimera.dep_dens_chnk(\
					  species.Data[component],species.Data['weights'],\
					  solver.Data['Rho'],species.chunks,\
					  species.Configs['Xchunked'][1],solver.Args['leftX'],\
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
			  or 'Still' not in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env_chnk(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,\
					  species.Configs['Xchunked'][1],solver.Args['leftX'],\
					  *solver.Args['DepProj'])
				else:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'], solver.Args['leftX'],\
					  *solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.Data['BckGrndRho'] = chimera.dep_dens_chnk(\
					  species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,\
					  species.Configs['Xchunked'][1],solver.Args['leftX'],\
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
			if 'StaticKick' in solver.Configs['Features']: continue
			solver.damp_field()

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

	def damp_plasma(self,wind):
		if 'AbsorbLayer' not in wind: return
		for species in self.Particles: species.damp_particles(wind)

	def postframe_corr(self,wind):
		if 'AbsorbLayer' or 'AddPlasma' in wind:
			for solver in self.Solvers:
				if 'StaticKick' in solver.Configs['Features']: continue
				if 'SpaceCharge' in solver.Configs['Features']:
					if 'StillAsBackground' in solver.Configs['Features']:	
						self.dep_bg(solver)
					self.dep_dens(solver)

	def move_frame(self,wind):
		for comp in self.Solvers+self.Particles:
			comp.Args['Xgrid']  += wind['shiftX']
			comp.Args['leftX']  = comp.Args['Xgrid'][ 0]
			comp.Args['rightX'] = comp.Args['Xgrid'][-1]

	def frame_act(self,istep,act='stage1'):
		for wind in self.MovingFrames:
			if istep<wind['TimeActive'][0] or istep>wind['TimeActive'][1]: 
				continue
			if act=='stage1': self.damp_fields(wind)
			if np.mod( istep-wind['TimeActive'][0], wind['Steps'])!= 0: continue
			if act=='stage1':
#				self.damp_fields(wind)
				self.move_frame(wind)
				self.add_plasma(wind)
				self.damp_plasma(wind)
				self.postframe_corr(wind)
			elif act=='stage2' and 'Staged' in wind['Features']:
				self.move_frame(wind)

	def chunk_particles(self,istep=0):
		for species in self.Particles:
			if  'Xchunked' not in species.Configs \
			  or 'NoSorting' in species.Configs['Features']: continue
			if np.mod(istep, species.Configs['Xchunked'][1]+1)!= 0: continue
			if species.Data['coords'].shape[-1] == 0: continue
		species.chunk_coords('cntr')
