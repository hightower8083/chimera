import numpy as np
import sys, os,time
sys.path.append('./src/')
sys.path.append('./moduls/')
from species import Specie
from solvers import Solver
from chimera_main import ChimeraRun
from RunRestart import StartFromSave
import cPickle as pickle
runsave = StartFromSave()

out_folder = '../small_bubble1/'
dt = 0.04
Steps2Do = 601
BoxGrid = (-50.,0.0,20.,0.04,0.25)

densprof = lambda x: np.interp(x, [0.,40.,50.,300],[0.,1.,0.5,0.5])

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('SpaceCharge','StillAsBackground','Xchunked')}
seed_in = {'a0':3.0,'k0':1.0,'x0':-16.0,'x_foc':45.0,'Lx':4,'LR':4}

specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.005, 'Mass':1., 'Grid':BoxGrid, 'TimeStep':dt, 'Xchunked':(8,41), 
'Features':('NoSorting',)}
specie2_in = {'FixedCell':(2,2,4), 'Charge':1., 'Density':0.005, 'Mass':1886., 'Grid':BoxGrid, 'TimeStep':dt,'Xchunked':(8,41), 
'Features':('Still','NoSorting')}

MovingFrame = {'TimeStep':dt,'Steps':40,'AbsorbLayer':10.,'AddPlasma':densprof,'Features':()}

fld_out_step = 200
dns_out_step = 200
phs_out_step = None

if sys.argv[-1]=='sim':
	solver = Solver(solver_in)
	solver.add_gauss_beam(seed_in)
	specie1 = Specie(specie1_in)
	specie2 = Specie(specie2_in)
	chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':(MovingFrame,)}
	Chimera = ChimeraRun(chimera_in)
	Chimera.make_halfstep()

	timestart = time.time()
	os.system('rm -rf '+out_folder)
	os.system('mkdir ' +out_folder)
	for i in xrange(Steps2Do):
		Chimera.make_step(i)

		istr = str(i)
		while len(istr)<7: istr='0'+istr
		print istr

		if fld_out_step!= None and np.mod(i,fld_out_step)==0:
			for jj in range(len(Chimera.Solvers)):
				sol = Chimera.Solvers[jj]
				np.save(out_folder+'ee_'+str(jj)+'_'+istr+'.npy',solver.EB[:,1:])

		if dns_out_step!=None and np.mod(i,dns_out_step)==0:
			for jj in range(len(Chimera.Particles)):
				if 'Still' in Chimera.Particles[jj].Configs['Features']: continue
				np.save(out_folder+'edens_'+str(jj)+'_'+istr+'.npy', Chimera.Particles[jj].get_dens_on_grid(solver_in['MaxAzimuthMode'])[:,1:])

		if phs_out_step!=None and np.mod(i,phs_out_step)==0:
			for jj in range(len(Chimera.Particles)):
				if 'Still' in Chimera.Particles[jj].Configs['Features']: continue
				np.save(out_folder+'phs'+str(jj)+'_'+istr+'.npy',Chimera.Particles[jj].particles)

	runsave.SaveRun(Chimera)
	runsave.Step2Start = i+1
	file2save = open(out_folder+'savedrun_'+'.p','wb')
	pickle.dump(runsave, file2save, pickle.HIGHEST_PROTOCOL)
	print 'done in %f minutes' % ((time.time()-timestart)/60.,)
