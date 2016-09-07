import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from chimera_tools import *

K0 = 1.71856	#
Lu = 2.			# cm
Periods = 102

g0 = 352.2505 	# gamma eons
Lx = 1.e-4		# cm
LR = 2.84e-4	# cm

StepsPerPeriod = 40 #

k_res = 2*g0**2/(1+K0**2/2)
gg = g0/(1.+K0**2/2)
vb = (1.-gg**-2)**0.5
Lx = Lx/Lu
LR = LR/Lu

out_folder = './../fel_test1/'
dt = 1./StepsPerPeriod
Steps2Do = int(Periods/dt)+1

fld_out_step = None
dns_out_step = None
fldfb_out_step = 500
phs_out_step = 500
nrg_out_step = int(1./dt)

densprof = lambda x,y,z: np.exp(-0.5*x**2/Lx**2)*np.exp(-0.5*(y**2+z**2)/LR**2)

BoxGrid     = (-50*Lx ,50*Lx ,450*LR, 0.5*Lx, 2.5*LR)
BoxGridBeam = (-4*Lx,4*Lx, 4*LR,0.1*Lx,0.05*LR)

solver_in = {'MaxAzimuthMode':0, 'Grid':BoxGrid, 'KxShift':  k_res,'TimeStep':dt, 'Rcut':250*LR, 'Features':('NoPoissonCorrection',)}

seed_in = {'a0':0.076,'k0':k_res,'x0':0.0,'x_foc':10.0,'Lx':0.0016,'LR':0.00892}

Devices = ({'DEV':['MagneticUndulator1D',K0, np.array([1.,0.,Periods,1.]),np.array([1.,-0.027,0.027])]},)

specie_in = {'FixedCell':(5,5,6), 'Charge':-1., 'Density':6e5, 'Mass':1., 'Grid':BoxGridBeam, 'TimeStep':dt,
'MomentaMeans':(g0,0.,0.), 'MomentaSpreads':(3.5225,0.35225,0.35225), 'NoiseRemv':(k_res,2*k_res),
'Features':(),'Devices':Devices}

MovingFrame = {'TimeStep':dt,'Steps':1,'Velocity':vb,'Features':('Staged',)}

solver = Solver(solver_in)
solver.add_gauss_beam(seed_in)
specie = Specie(specie_in)
specie.add_particles(specie.gen_parts(Domain=(-2e-4,2e-4,0.0,0.0005),ProfileFunc=densprof))
specie.correct_fel()
chimera_in = {'Solvers':(solver,),'Particles':(specie,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(chimera_in)
Chimera.make_halfstep()

if sys.argv[-1]=='sim':
	timestart = time.time()
	os.system('rm -rf '+out_folder)
	os.system('mkdir ' +out_folder)
	for i in xrange(Steps2Do):
		Chimera.make_step(i)
		if np.mod(i,Steps2Do/10)==0: print i

		if nrg_out_step != None and np.mod(i,nrg_out_step)==0:
			for jj in range(len(Chimera.Solvers)):
				sol = Chimera.Solvers[jj]
				fout = open(out_folder+'nrg_'+str(jj)+'.txt', mode='a')
				fout.write(  str((abs(sol.EG_fb[:,:,:,:3])**2*sol.Args['EnergyFact']).sum())+'\n')
				fout.close()
		if fldfb_out_step != None and np.mod(i,fldfb_out_step)==0:
			for jj in range(len(Chimera.Solvers)):
				sol = Chimera.Solvers[jj]
				istr = str(i)
				while len(istr)<7: istr='0'+istr
				np.save(out_folder+'eefb_'+istr+'.npy',sol.EG_fb)
		if phs_out_step != None and np.mod(i,phs_out_step)==0:
			for jj in range(len(Chimera.Particles)):
				if 'Still' in specie.Configs['Features']: continue
				istr = str(i)
				while len(istr)<7: istr='0'+istr
				np.save(out_folder+'phs_'+str(jj)+'_'+istr+'.npy', Chimera.Particles[jj].particles)

	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
