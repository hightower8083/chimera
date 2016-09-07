class StartFromSave:
	def __init__(self):
		self.Particles = []
		self.Solvers   = []
		self.Step2Start= 0
	def SaveRun(self,Chimera):
		for specie in Chimera.Particles:
			if 'Xchunked' in specie.Configs:
				self.Particles.append((specie.particles, specie.particles_cntr,specie.chunks, \
				  specie.EB, specie.Args['leftX'], specie.Args['rightX'],specie.Args['Xgrid']))
			else:
				self.Particles.append((specie.particles, specie.particles_cntr, \
				  specie.EB, specie.Args['leftX'], specie.Args['rightX'],specie.Args['Xgrid']))
		for solver in Chimera.Solvers:
			if 'SpaceCharge' in solver.Configs['Features']:
				if 'StillAsBackground' in solver.Configs['Features']:
					self.Solvers.append((\
					  solver.EB,solver.EG_fb,solver.vec_fb_aux,solver.vec_fb_aux1,solver.bg_spc,\
					  solver.Args['leftX'], solver.Args['rightX'],solver.Args['Xgrid']\
					  ))
				else:
					self.Solvers.append((\
					  solver.EB,solver.EG_fb,solver.vec_fb_aux,solver.vec_fb_aux1,\
					  solver.Args['leftX'], solver.Args['rightX'],solver.Args['Xgrid']\
					  ))
			else:
				self.Solvers.append((\
				  solver.EB,solver.EG_fb,solver,solver.vec_fb_aux, solver.Args['leftX'], \
				  solver.Args['rightX'],solver.Args['Xgrid']
				  ))

	def LoadRun(self,Chimera):
		for i in range(len(Chimera.Particles)):
			specie = Chimera.Particles[i]
			if 'Xchunked' in specie.Configs:
				specie.particles, specie.particles_cntr,specie.chunks,specie.EB, specie.Args['leftX'], \
				  specie.Args['rightX'],specie.Args['Xgrid'] = self.Particles[i]
			else:
				specie.particles, specie.particles_cntr,specie.EB, specie.Args['leftX'], \
				  specie.Args['rightX'],specie.Args['Xgrid'] = self.Particles[i]
		for i in range(len(Chimera.Solvers)):
			solver = Chimera.Solvers[i]
			if 'SpaceCharge' in solver.Configs['Features']:
				if 'StillAsBackground' in solver.Configs['Features']:
					solver.EB,solver.EG_fb,solver.vec_fb_aux, solver.vec_fb_aux1,solver.bg_spc,\
					  solver.Args['leftX'], solver.Args['rightX'],solver.Args['Xgrid'] = self.Solvers[i]
				else:
					solver.EB,solver.EG_fb,solver.vec_fb_aux, solver.vec_fb_aux1,\
					  solver.Args['leftX'], solver.Args['rightX'],solver.Args['Xgrid'] = self.Solvers[i]
			else:
				solver.EB,solver.EG_fb,solver.vec_fb_aux,solver.Args['leftX'], solver.Args['rightX'],\
				  solver.Args['Xgrid'] = self.Solvers[i]
