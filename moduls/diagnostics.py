from __future__ import print_function,division
import numpy as np
import os
import chimera.moduls.fimera as chimera
#import fimera as chimera

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
				np.save(self.out_folder+'ee_'+str(jj)+'_'+self.istr \
				  + '.npy',self.Chimera.Solvers[jj].Data['EB'][:,1:])
		if 'Return' in diag['Features']: return ToReturn

	def fldfb_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			if 'Return' in diag['Features']:
				ToReturn.append(self.Chimera.Solvers[jj].Data['EG_fb'].copy())
			else:
				np.save(self.out_folder+'ee_'+str(jj)+'_'+self.istr \
				  + '.npy',self.Chimera.Solvers[jj].Data['EG_fb'])
		if 'Return' in diag['Features']: return ToReturn

	def dns_out(self,diag):
		if 'MaxMode' in diag['Features']:
			modes = diag['Features']['MaxMode']
		else:
			if len(self.Chimera.Solvers)>0:
				modes = self.Chimera.Solvers[0].Configs['MaxAzimuthMode']
			else:
				modes = 0
		if 'Return' in diag['Features']:ToReturn = []
		for jj in range(len(self.Chimera.Particles)):
			if 'Still' in self.Chimera.Particles[jj].Configs['Features']: continue
			dat = self.Chimera.Particles[jj].get_dens_on_grid(modes)[:,1:]
			if 'Return' in diag['Features']:
				ToReturn.append(dat.copy())
			else:
				np.save(self.out_folder+'dens_'+str(jj)+'_'+self.istr+'.npy',dat)
			dat = None
		if 'Return' in diag['Features']: return ToReturn

	def phs_out(self,diag):
		if 'Return' in diag['Features']: ToReturn = []
		for jj in range(len(self.Chimera.Particles)):
			if 'Still' in self.Chimera.Particles[jj].Configs['Features']: continue
			indx = np.nonzero( \
			  (1.+self.Chimera.Particles[jj].Data['momenta']**2).sum(0) \
			  > 20.**2)[0] # LPA filter
			PartsPack = np.concatenate((\
			  self.Chimera.Particles[jj].Data['coords_halfstep'][:,indx],\
			  self.Chimera.Particles[jj].Data['momenta'][:,indx],\
			  self.Chimera.Particles[jj].Data['weights'][indx][None,:]),axis=0)
			if 'Return' in diag['Features']:
				ToReturn.append(PartsPack.copy())
			else:
				np.save(self.out_folder+'phs'+str(jj)+'_'+self.istr \
				  + '.npy',PartsPack)
			PartsPack = None
			indx = None
		if 'Return' in diag['Features']: return ToReturn

	def nrg_out(self,diag):
		if 'Return' in diag['Features']:	ToReturn = []
		for jj in range(len(self.Chimera.Solvers)):
			sol = self.Chimera.Solvers[jj]
			dat = ((abs(sol.Data['EG_fb'][:,:,:,:3])**2).sum(-1) \
			  * sol.Args['EnergyFact']).sum(-1).sum(-1)
			dat = np.r_[dat[dat.shape[0]//2+1:], dat[:dat.shape[0]//2+1]]
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
			dat = chimera.fb_vec_out(sol.Data['EG_fb'][:,:,:,:3], \
			  sol.Args['leftX'],*sol.Args['FBoutFull'])

			if 'Return' in diag['Features']:
				if 'Spot' in diag['Features']:
					ToReturn.append([  PwrFctr*2*dr*\
					  ((np.abs(dat)**2).sum(-1).sum(-1)*Rgrid[None,:]).sum(-1),\
					  chimera.intens_profo(dat,Ntheta)  ])
				else:
					ToReturn.append(PwrFctr*2*dr*\
					  ((np.abs(dat)**2).sum(-1).sum(-1)*Rgrid[None,:]).sum(-1))
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

	def get_beam_envelops(self):
		ToReturn = []
		for jj in range(len(self.Chimera.Particles)):
			if 'Still' in self.Chimera.Particles[jj].Configs['Features']: 
				continue
			x,y,z = self.Chimera.Particles[jj].Data['coords_halfstep']
			px,py,pz = self.Chimera.Particles[jj].Data['momenta']
			w = self.Chimera.Particles[jj].Data['weights']

			xyz_0 = [\
			  (x*w).sum()/w.sum(),
			  (y*w).sum()/w.sum(),
			  (z*w).sum()/w.sum()]

			xyz_rms = [\
			  np.sqrt((x**2*w).sum()/w.sum()-xyz_0[0]**2),
			  np.sqrt((y**2*w).sum()/w.sum()-xyz_0[1]**2),
			  np.sqrt((z**2*w).sum()/w.sum()-xyz_0[2]**2)]

			xyz_emmit = [\
			  np.sqrt( (x**2*w).sum()/w.sum()*(px**2*w).sum()/w.sum() \
			          -(x*px*w).sum()**2/w.sum()**2 ),
			  np.sqrt( (y**2*w).sum()/w.sum()*(py**2*w).sum()/w.sum() \
			          -(y*py*w).sum()**2/w.sum()**2),
			  np.sqrt( (z**2*w).sum()/w.sum()*(pz**2*w).sum()/w.sum() \
			          -(z*pz*w).sum()**2/w.sum()**2)]

			PackEnvs = np.array([xyz_0,xyz_rms,xyz_emmit])
			ToReturn.append(PackEnvs)
			PackEnvs = None
		return ToReturn
