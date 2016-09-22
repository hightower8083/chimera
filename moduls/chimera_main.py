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
		for wind in self.MovingFrames: self.move_frame_staged(wind,istep)
		self.update_fields()
		for species in self.Particles: species.make_field(istep)
		self.project_fields()
		for species in self.Particles: species.push_velocs()

	def make_halfstep(self):
		for species in self.Particles: species.make_field()
		for species in self.Particles:
			if species.particles.shape[1] != 0: species.chunk_coords()
		self.get_init_fields()
		for solver in self.Solvers: solver.G2B_FBRot()
		self.project_fields()
		for species in self.Particles: species.push_velocs()

	def project_current(self):
		for solver in self.Solvers:
			if 'StaticKick' in solver.Configs['Features']: continue
			self.dep_curr_on_grid(solver)
			solver.fb_curr_in()

	def project_density(self):
		for solver in self.Solvers:
			if 'SpaceCharge' not in solver.Configs['Features']: continue
			if 'StaticKick' in solver.Configs['Features']: continue
			solver.vec_fb_aux0[:] = solver.vec_fb_aux1.copy()
			self.dep_dens_on_grid(solver)
			solver.fb_dens_in()

	def update_fields(self):
		#self.re_psatd()
		for solver in self.Solvers:
			if 'StaticKick' in solver.Configs['Features']:
				self.get_static_kick(solver)
				solver.G2B_FBRot()
			else:
				solver.poiss_corr()
				solver.maxwell_solver()
				solver.G2B_FBRot()

	def project_fields(self):
		for solver in self.Solvers:
			solver.fb_fld_out()
			for species in self.Particles:
				if species.particles.shape[1]==0 or 'Still' in species.Configs['Features']: continue
				if 'StaticKick' in solver.Configs['Features']:
					species.EB = chimera.proj_fld_cntr(species.particles_cntr,species.particles[-1],solver.EB,\
					  species.EB,solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					if 'KxShift' in solver.Configs:
						species.EB = chimera.proj_fld_env(species.particles,solver.EB,\
						  species.EB,solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						species.EB = chimera.proj_fld(species.particles,solver.EB,\
						  species.EB,solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_curr_on_grid(self,solver):
		solver.vec_spc[:] = 0.0
		for species in self.Particles:
			if species.particles.shape[1] == 0 or 'Still' in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.vec_spc = chimera.dep_curr_env_chnk(species.particles,\
					  species.particles_cntr,solver.vec_spc,species.chunks,species.Configs['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.vec_spc = chimera.dep_curr_env(species.particles,\
					  species.particles_cntr,solver.vec_spc,solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.vec_spc = chimera.dep_curr_chnk(species.particles,\
					  species.particles_cntr,solver.vec_spc,species.chunks,species.Configs['Xchunked'][1],\
					  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.vec_spc = chimera.dep_curr(species.particles,\
					  species.particles_cntr,solver.vec_spc,solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_dens_on_grid(self,solver):
		solver.scl_spc[:] = 0.0
		for species in self.Particles:
			if species.particles.shape[1] == 0: continue
			if 'Still' in species.Configs['Features'] and 'StillAsBackground' in solver.Configs['Features']:
				solver.scl_spc += solver.bg_spc
			else:
				if 'KxShift' in solver.Configs:
					if 'Xchunked' in species.Configs:
						solver.scl_spc = chimera.dep_dens_env_chnk(species.particles,solver.scl_spc,\
						  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.scl_spc = chimera.dep_dens_env(species.particles,solver.scl_spc,\
						  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					if 'Xchunked' in species.Configs:
						solver.scl_spc = chimera.dep_dens_chnk(species.particles,solver.scl_spc,\
						  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.scl_spc = chimera.dep_dens(species.particles,solver.scl_spc,\
						  solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_bg_on_grid(self,solver):
		solver.bg_spc[:] = 0.0
		for species in self.Particles:
			if species.particles.shape[1] == 0 or 'Still' not in species.Configs['Features']: continue
			if 'KxShift' in solver.Configs:
				if 'Xchunked' in species.Configs:
					solver.bg_spc = chimera.dep_dens_env_chnk(species.particles,solver.bg_spc,\
					  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.bg_spc = chimera.dep_dens_env(species.particles,solver.bg_spc,\
					  solver.Args['leftX'],*solver.Args['DepProj'])
			else:
				if 'Xchunked' in species.Configs:
					solver.bg_spc = chimera.dep_dens_chnk(species.particles,solver.bg_spc,\
					  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					solver.bg_spc = chimera.dep_dens(species.particles,solver.bg_spc,\
					  solver.Args['leftX'],*solver.Args['DepProj'])

	def dep_cntr_dens_curr_on_grid(self,solver,species):
		solver.scl_spc[:] = 0.0
		solver.scl_fb[:] = 0.0
		solver.vec_fb_aux1[:] = 0.0
		if species.particles.shape[1]==0: return
		if 'Xchunked' in species.Configs:
			solver.scl_spc = chimera.dep_dens_cntr_chnk(species.particles_cntr,species.particles[-1],solver.scl_spc,\
			  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
		else:
			solver.scl_spc = chimera.dep_dens_cntr(species.particles_cntr,species.particles[-1],solver.scl_spc,\
			  solver.Args['leftX'],*solver.Args['DepProj'])
		solver.fb_dens_in()

		solver.vec_fb[:] = 0.0
		solver.vec_spc[:] = 0.0
		if species.particles.shape[1] == 0 or 'Still' in species.Configs['Features']: return
		if 'Xchunked' in species.Configs:
			solver.vec_spc = chimera.dep_curr_chnk(species.particles,\
			  species.particles_cntr,solver.vec_spc,species.chunks,species.Configs['Xchunked'][1],\
			  solver.Args['leftX'],*solver.Args['DepProj'])
		else:
			solver.vec_spc = chimera.dep_curr(species.particles,\
			  species.particles_cntr,solver.vec_spc,solver.Args['leftX'],*solver.Args['DepProj'])
		solver.fb_curr_in()

	def get_static_kick(self,solver):
		solver.EG_fb[:] = 0.0
		for species in self.Particles:
			self.dep_cntr_dens_curr_on_grid(solver,species)
			PXmean = (species.particles[3]*species.particles[-1]).sum()/species.particles[-1].sum()
			solver.maxwell_solver_init(PXmean)
		solver.vec_fb[:] = 0.0
		solver.scl_fb[:] = 0.0
		solver.vec_fb_aux1[:] = 0.0
		solver.vec_fb_aux[:] = 0.0

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
		elif wind['Velocity'] == 'Follow':
			vw = (self.Particles[0].particles[3]/self.Particles[0].particles[-2]*\
			  self.Particles[0].particles[-1]).sum()/self.Particles[0].particles[-1].sum()
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
			####  EONS ON-TOP OF IONS
			specie1, specie2 = self.Particles
			rightX = specie1.Args['rightX']
			lowerR, upperR  = specie1.Args['lowerR'], specie1.Args['upperR']
			parts_to_add = specie1.gen_parts(Xsteps=int(shiftX/wind['TimeStep'])+1, ProfileFunc=wind['AddPlasma'])
			specie1.add_particles(parts_to_add)
			parts_to_add[-1] *= -1
			specie2.add_particles(parts_to_add)
			### GENERAL CASE
#			for species in self.Particles:
#				rightX = species.Args['rightX']
#				lowerR, upperR  = species.Args['lowerR'], species.Args['upperR']
#				species.add_particles(species.gen_parts((rightX-shiftX,rightX,lowerR,upperR),\
#				  wind['AddPlasma']))
		if 'AbsorbLayer' in wind:
			for species in self.Particles:
				if species.particles.shape[1]==0: continue
				SimDom = np.asfortranarray([species.Args['leftX']+wind['AbsorbLayer'],species.Args['rightX'],\
				  0.0, species.Args['upperR']**2])
				index2stay,num2stay = chimera.sortpartsout(np.asfortranarray(species.particles[:3]),SimDom)
				index2stay = index2stay[:num2stay]

				species.particles, species.particles_cntr = \
				  chimera.align_coords(species.particles, species.particles_cntr,index2stay)
				species.particles = species.particles[:,:index2stay.shape[0]]
				species.particles_cntr = species.particles_cntr[:,:index2stay.shape[0]]
		if 'NoSorting' not in wind['Features']:
			for species in self.Particles:
				if species.particles.shape[1] != 0: species.chunk_coords()
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
		elif wind['Velocity'] == 'Follow':
			vw = (self.Particles[0].particles[3]/self.Particles[0].particles[-2]*\
			  self.Particles[0].particles[-1]).sum()/self.Particles[0].particles[-1].sum()
		else:
			vw = wind['Velocity']
		shiftX = 0.5*vw*wind['TimeStep']*wind['Steps']
		for comp in self.Solvers+self.Particles:
			comp.Args['Xgrid']  += shiftX
			comp.Args['leftX']  = comp.Args['Xgrid'][ 0]
			comp.Args['rightX'] = comp.Args['Xgrid'][-1]

	def get_init_fields(self):
		for solver in self.Solvers:
			if 'SpaceCharge' not in solver.Configs['Features'] and 'StaticKick' not in solver.Configs['Features']: continue
			print 'adding inits'
			for species in self.Particles:
				solver.vec_spc[:] = 0.0
				solver.scl_spc[:] = 0.0
				if species.particles.shape[1]==0: continue

				if 'KxShift' in solver.Configs:
					if 'Xchunked' in species.Configs:
						solver.scl_spc = chimera.dep_dens_env_chnk(species.particles,solver.scl_spc,\
						  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.scl_spc = chimera.dep_dens_env(species.particles,solver.scl_spc,\
						  solver.Args['leftX'],*solver.Args['DepProj'])
				else:
					if 'Xchunked' in species.Configs:
						solver.scl_spc = chimera.dep_dens_chnk(species.particles,solver.scl_spc,\
						  species.chunks,species.Configs['Xchunked'][1],solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						solver.scl_spc = chimera.dep_dens(species.particles,solver.scl_spc,\
						  solver.Args['leftX'],*solver.Args['DepProj'])
				solver.fb_dens_in()
				if 'Still' not in species.Configs['Features']:
					if 'KxShift' in solver.Configs:
						if 'Xchunked' in species.Configs:
							solver.vec_spc = chimera.dep_curr_env_chnk(species.particles,\
							  species.particles_cntr,solver.vec_spc,species.chunks,species.Configs['Xchunked'][1],\
							  solver.Args['leftX'],*solver.Args['DepProj'])
						else:
							solver.vec_spc = chimera.dep_curr_env(species.particles,\
							  species.particles_cntr,solver.vec_spc,solver.Args['leftX'],*solver.Args['DepProj'])
					else:
						if 'Xchunked' in species.Configs:
							solver.vec_spc = chimera.dep_curr_chnk(species.particles,\
							  species.particles_cntr,solver.vec_spc,species.chunks,species.Configs['Xchunked'][1],\
							  solver.Args['leftX'],*solver.Args['DepProj'])
						else:
							solver.vec_spc = chimera.dep_curr(species.particles,\
							  species.particles_cntr,solver.vec_spc,solver.Args['leftX'],*solver.Args['DepProj'])
					solver.fb_curr_in()
				elif 'StillAsBackground' in solver.Configs['Features']:
					solver.bg_spc += solver.scl_spc
				solver.maxwell_solver_init(species.Configs['MomentaMeans'][0])
				solver.vec_spc[:] = 0.0
				solver.scl_spc[:] = 0.0
				solver.scl_fb[:] = 0.0
				solver.vec_fb[:] = 0.0
				if 'SpaceCharge' in solver.Configs['Features']:	solver.vec_fb_aux0[:] = 0.0

	def chunk_particles(self,istep=0):
		for species in self.Particles:
			if  'Xchunked' not in species.Configs or 'NoSorting' in species.Configs['Features']: continue
			if np.mod(istep, species.Configs['Xchunked'][1]+1)!= 0: continue
			if 'Still' in species.Configs['Features']: continue
			if species.particles.shape[1] == 0: continue
			species.chunk_coords('cntr')

	def re_psatd(self):
		for species in self.Particles:
			vb = (species.particles[3]/species.particles[-2]*species.particles[-1]).sum()/species.particles[-1].sum()
			for solver in self.Solvers:
				if 'KxShift' in solver.Configs:
					solver.PSATD_coeffs(vb)
