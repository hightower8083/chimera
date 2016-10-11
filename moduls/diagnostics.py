import numpy as np
import os
import fimera as chimera

PwrFctr = 0.5*0.511e6*1.6022e-19/2.818e-13*2.9979e10
Ntheta = 60

class Diagnostics:
	def __init__(self,Chimera, diags=(), out_folder=None):
		self.Chimera = Chimera
		self.diags = diags
		if out_folder!=None:
			self.out_folder = out_folder
			os.system('rm -rf ' + self.out_folder)
			os.system('mkdir '  + self.out_folder)

	def do_diags(self,i):
		for diag in self.diags:
			if np.mod(i,diag['Step'])!=0: continue
			self.istr = str(i)
			while len(self.istr)<7: self.istr='0'+self.istr
			if 'Features' not in diag: diag['Features'] = {}

			if diag['Type']=='Fields'   : self.fld_out(diag)
			if diag['Type']=='FieldsFB' : self.fldfb_out(diag)
			if diag['Type']=='Particles': self.phs_out(diag)
			if diag['Type']=='Density'  : self.dns_out(diag)
			if diag['Type']=='EnergyEM' : self.nrg_out(diag)
			if diag['Type']=='Power'    : self.pwr_out(diag)

	def fld_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			if 'Return' in diag['Features']:
				ToReturn.append(self.Chimera.Solvers[jj].Data['EB'][:,1:].copy())
			else:
				np.save(self.out_folder+'ee_'+str(jj)+'_'+self.istr+'.npy',self.Chimera.Solvers[jj].Data['EB'][:,1:])
		if 'Return' in diag['Features']: return ToReturn

	def fldfb_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			if 'Return' in diag['Features']:
				ToReturn.append(self.Chimera.Solvers[jj].Data['EG_fb'].copy())
			else:
				np.save(self.out_folder+'ee_'+str(jj)+'_'+self.istr+'.npy',self.Chimera.Solvers[jj].Data['EG_fb'])
		if 'Return' in diag['Features']: return ToReturn

	def dns_out(self,diag):
		if 'MaxMode' in diag['Features']:
			modes = diag['Features']['MaxMode']
		else:
			modes = 0
		if 'Return' in diag['Features']:ToReturn = []
		for jj in range(len(self.Chimera.Particles)):
			if 'Still' in self.Chimera.Particles[jj].Configs['Features']: continue
			if 'Return' in diag['Features']:
				ToReturn.append( self.Chimera.Particles[jj].get_dens_on_grid(modes)[:,1:].copy())
			else:
				np.save(self.out_folder+'dens_'+str(jj)+'_'+self.istr+'.npy', \
				  self.Chimera.Particles[jj].get_dens_on_grid(modes)[:,1:])
		if 'Return' in diag['Features']: return ToReturn

	def phs_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Particles)):
			if 'Still' in self.Chimera.Particles[jj].Configs['Features']: continue
#			indx = np.nonzero(self.Chimera.Particles[jj].particles[-2]>20.)[0] # LPA filter
			if 'Return' in diag['Features']:
				ToReturn.append(self.Chimera.Particles[jj].particles.copy())
			else:
				np.save(self.out_folder+'phs'+str(jj)+'_'+self.istr+'.npy',self.Chimera.Particles[jj].particles)
		if 'Return' in diag['Features']: return ToReturn

	def nrg_out(self,diag):
		if 'Return' in diag['Features']:	ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			sol = self.Chimera.Solvers[jj]
			dat = ((abs(sol.Data['EG_fb'][:,:,:,:3])**2).sum(-1)*sol.Args['EnergyFact']).sum(-1).sum(-1)
			dat = np.r_[dat[dat.shape[0]/2+1:], dat[:dat.shape[0]/2+1]]
			if 'Return' in diag['Features']:
				ToReturn.append(dat.copy())
			else:
				fout = open(self.out_folder+'nrg_'+str(jj)+'.txt', mode='a')
				np.savetxt(fout,dat[None,:])
				fout.close()
		if 'Return' in diag['Features']: return ToReturn

	def pwr_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			sol = self.Chimera.Solvers[jj]
			Rgrid,dr = sol.Args['RgridFull'],sol.Args['dr']
			dat = chimera.fb_vec_out(sol.Data['EG_fb'][:,:,:,:3],sol.Args['leftX'],*sol.Args['FBoutFull'])

			if 'Return' in diag['Features']:
				if 'Spot' in diag['Features']:
					ToReturn.append([PwrFctr*2*dr*((np.abs(dat)**2).sum(-1).sum(-1)*Rgrid[None,:]).sum(-1),\
					  chimera.intens_profo(dat,Ntheta)])
				else:
					ToReturn.append(PwrFctr*2*dr*((np.abs(dat)**2).sum(-1).sum(-1)*Rgrid[None,:]).sum(-1))
			else:
				fout = open(self.out_folder+'pwrX_'+str(jj)+'.txt',mode='a')
				np.savetxt(fout,PwrFctr*2*dr*((np.abs(dat)**2).sum(-1).sum(-1)*\
				  Rgrid[None,:]).sum(-1)[None])
				fout.close()
				if 'Spot' in diag['Features']:
					fout = open(self.out_folder+'pwrO_'+str(jj)+'.txt',mode='a')
					np.savetxt(fout,chimera.intens_profo(dat,Ntheta).ravel()[None])
					fout.close()
		if 'Return' in diag['Features']: return ToReturn
