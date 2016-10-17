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
		else:
			self.MovingFrames = ()

	def make_step(self,istep):
		for wind in self.MovingFrames: self.move_frame(wind,istep)
		for species in self.Particles: species.push_coords()
		self.chunk_particles(istep)
		self.project_current()
		self.project_density()
		self.update_fields()
		for wind in self.MovingFrames: self.move_frame_staged(wind,istep)
		for species in self.Particles: species.make_field(istep)
		self.project_fields()
		for species in self.Particles: species.push_velocs()

	def make_halfstep(self):
		for species in self.Particles:
			if species.Data['coords'].shape[1] != 0: species.chunk_coords()
		for species in self.Particles: species.make_field()
		self.project_current()
		self.project_density()
		for solver in self.Solvers:
			solver.maxwell_solver_stat(species.Configs['MomentaMeans'][0])
			solver.G2B_FBRot()
		self.project_fields()
		for species in self.Particles: species.push_velocs(dt=0.5*species.Configs['TimeStep'])

	def project_current(self):
		for solver in self.Solvers:
			self.dep_curr_on_grid(solver)
			solver.fb_curr_in()

	def project_density(self):
		for solver in self.Solvers:
			if 'SpaceCharge' not in solver.Configs['Features'] and \
			  'StaticKick' not in solver.Configs['Features']: continue
			if 'StaticKick' not in solver.Configs['Features']:
#				solver.Data['vec_fb'][:] = solver.Data['EG_fb'][:,:,:,:3]
#				solver.FBDiv()
#				solver.Data['Rho_fb'][:] = solver.Data['scl_fb']
#				solver.FBGradDens()
				solver.Data['gradRho_fb_prv'][:] = solver.Data['gradRho_fb_nxt']
			self.dep_dens_on_grid(solver)
			solver.fb_dens_in()

	def update_fields(self):
		for solver in self.Solvers:
			if 'StaticKick' in solver.Configs['Features']:
				solver.Data['EG_fb'][:] = 0.0
				for species in self.Particles:
					PXmean = (species.Data['momenta'][0]*species.Data['weights']).sum()/species.Data['weights'].sum()
					solver.poiss_corr_stat(PXmean)
					solver.maxwell_solver_stat(PXmean)
					solver.field_drift(PXmean)
				solver.G2B_FBRot()
			else:
				solver.poiss_corr()
				solver.maxwell_solver()
				solver.G2B_FBRot()

	def project_fields(self,component='coords'):
		for solver in self.Solvers:
			solver.fb_fld_out()
			for species in self.Particles:
				if species.Data[component].shape[1]==0 or 'Still' in species.Configs['Features']: continue
				if 'KxShift' in solver.Configs:
					species.Data['EB'] = chimera.proj_fld_env(species.Data[component],species.Data['weights'],solver.Data['EB'],\
					  species.Data['EB'],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					species.Data['EB'] = chimera.proj_fld(species.Data[component],species.Data['weights'],solver.Data['EB'],\
					  species.Data['EB'],solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_curr_on_grid(self,solver):
		solver.Data['J'][:] = 0.0
		for species in self.Particles:
			if species.Data['coords'].shape[1] == 0 or 'Still' in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.Data['J'] = chimera.dep_curr_env_chnk(species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],species.chunks,species.Configs['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr_env(species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.Data['J'] = chimera.dep_curr_chnk(species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],species.chunks,species.Configs['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['J'] = chimera.dep_curr(species.Data['coords_halfstep'],species.Data['momenta'],\
					  species.Data['weights'],solver.Data['J'],solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_dens_on_grid(self,solver,component='coords'):
		solver.Data['Rho'][:] = 0.0
		if 'StaticKick' in solver.Configs['Features']: component='coords_halfstep'
		for species in self.Particles:
			if species.Data[component].shape[1] == 0: continue
			if 'Still' in species.Configs['Features'] and 'StillAsBackground' in solver.Configs['Features']:
				solver.Data['Rho'] += solver.Data['BckGrndRho']
			else:
				if 'KxShift' in solver.Configs:
					if 'Xchunked' in species.Configs:
						solver.Data['Rho'] = chimera.dep_dens_env_chnk(species.Data[component],species.Data['weights'],\
						  solver.Data['Rho'],species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.Data['Rho'] = chimera.dep_dens_env(species.Data[component],species.Data['weights'],\
						  solver.Data['Rho'],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					if 'Xchunked' in species.Configs:
						solver.Data['Rho'] = chimera.dep_dens_chnk(species.Data[component],species.Data['weights'],\
						  solver.Data['Rho'],species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.Data['Rho'] = chimera.dep_dens(species.Data[component],species.Data['weights'],\
						  solver.Data['Rho'],solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_bg_on_grid(self,solver):
		solver.Data['BckGrndRho'][:] = 0.0
		for species in self.Particles:
			if species.Data['coords'].shape[1]==0 or 'Still' not in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env_chnk(species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['BckGrndRho'] = chimera.dep_dens_env(species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'], solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.Data['BckGrndRho'] = chimera.dep_dens_chnk(species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.Data['BckGrndRho'] = chimera.dep_dens(species.Data['coords'],species.Data['weights'],\
					  solver.Data['BckGrndRho'],solver.Args['leftX'],*solver.Args['DepProj'])

	def move_frame(self,wind,istep):
		if 'TimeActive' in wind:
			WindInterval = wind['TimeActive']
		else:
			WindInterval = (0.0,np.inf)
		if istep<WindInterval[0] or istep>WindInterval[1]: return
		if np.mod( istep-WindInterval[0], wind['Steps'])!= 0: return
		if 'Features' not in wind: wind['Features'] = ()
		if 'Velocity' not in wind:
			vw = 1.0
		else:
			vw = wind['Velocity']

		if 'Staged' in wind['Features']:
			shiftX = 0.5*vw*wind['TimeStep']*wind['Steps']
		else:
			shiftX = vw*wind['TimeStep']*wind['Steps']

		if 'AbsorbLayer' in wind and wind['AbsorbLayer']>0:
			for solver in self.Solvers:
				if 'StaticKick' in solver.Configs['Features']: continue
				solver.absorb_field(wind['AbsorbLayer'],'left')

		for comp in self.Solvers+self.Particles:
			comp.Args['Xgrid']  += shiftX
			comp.Args['leftX']  = comp.Args['Xgrid'][ 0]
			comp.Args['rightX'] = comp.Args['Xgrid'][-1]

		if 'AddPlasma' in wind:
			if 'IonsOnTop' in wind['Features']:
				specie1, specie2 = self.Particles[:2]
				parts_to_add = specie1.gen_parts(Xsteps=int(shiftX/wind['TimeStep'])+1, ProfileFunc=wind['AddPlasma'])
				specie1.add_particles(*parts_to_add)
				parts_to_add[-1] *= -1
				specie2.add_particles(*parts_to_add)
			else:
				for species in self.Particles:
					species.add_particles(*species.gen_parts(Xsteps=int(shiftX/wind['TimeStep'])+1, \
					  ProfileFunc=wind['AddPlasma']))
		if 'AbsorbLayer' in wind:
			for species in self.Particles:
				if species.Data['coords'].shape[-1]==0: continue
				SimDom = np.asfortranarray([species.Args['leftX']+wind['AbsorbLayer'],species.Args['rightX'],\
				  0.0, species.Args['upperR']**2])
				index2stay,num2stay = chimera.sortpartsout(species.Data['coords'],SimDom)
				index2stay = index2stay[:num2stay]

				species.Data['coords'] = chimera.align_data_vec(species.Data['coords'],index2stay)
				species.Data['coords_halfstep'] = chimera.align_data_vec(species.Data['coords_halfstep'],index2stay)
				species.Data['momenta'] = chimera.align_data_vec(species.Data['momenta'],index2stay)
				species.Data['weights'] = chimera.align_data_scl(species.Data['weights'],index2stay)

				species.Data['coords'] = species.Data['coords'][:,:num2stay]
				species.Data['coords_halfstep'] = species.Data['coords_halfstep'][:,:num2stay]
				species.Data['momenta'] = species.Data['momenta'][:,:num2stay]
				species.Data['weights'] = species.Data['weights'][:num2stay]
		if 'NoSorting' not in wind['Features']:
			for species in self.Particles:
				if species.Data['coords'].shape[-1] != 0: species.chunk_coords()
		if 'AbsorbLayer' or 'AddPlasma' in wind:
			for solver in self.Solvers:
				if 'StaticKick' in solver.Configs['Features']: continue
				if 'SpaceCharge' in solver.Configs['Features']:
					if 'StillAsBackground' in solver.Configs['Features']:	self.dep_bg_on_grid(solver)
					self.dep_dens_on_grid(solver)
					solver.fb_dens_in()

	def move_frame_staged(self,wind,istep):
		if 'Features' not in wind: wind['Features'] = ()
		if 'Staged' not in wind['Features']: return
		if 'TimeActive' in wind:
			WindInterval = wind['TimeActive']
		else:
			WindInterval = (0.0,np.inf)
		if istep<WindInterval[0] or istep>WindInterval[1]: return
		if np.mod( istep-WindInterval[0], wind['Steps'])!= 0: return
		if 'Velocity' not in wind:
			vw = 1.0
		else:
			vw = wind['Velocity']
		shiftX = 0.5*vw*wind['TimeStep']*wind['Steps']
		for comp in self.Solvers+self.Particles:
			comp.Args['Xgrid']  += shiftX
			comp.Args['leftX']  = comp.Args['Xgrid'][ 0]
			comp.Args['rightX'] = comp.Args['Xgrid'][-1]

	def chunk_particles(self,istep=0):
		for species in self.Particles:
			if  'Xchunked' not in species.Configs or 'NoSorting' in species.Configs['Features']: continue
			if np.mod(istep, species.Configs['Xchunked'][1]+1)!= 0: continue
			if 'Still' in species.Configs['Features']: continue
			if species.Data['coords'].shape[-1] == 0: continue
			species.chunk_coords('cntr')
