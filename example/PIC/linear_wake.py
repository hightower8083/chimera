out_folder = '../linear_wake_01/'
dt = 0.075
Steps2Do = 5001
BoxGrid = (-90.,0.0,45.,0.075,0.4)

densprof = np.vectorize(lambda x,y,z: (x>=0.)*(x<=50.)*x/50.+ (x>50.))

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('Rsliced','SpaceCharge',)} 
seed_in = {'a0':0.01,'k0':1.0,'x0':-45.0,'x_foc':150.0,'Lx':12.0,'LR':12.0}

specie1_in = {'FixedCell':(3,3,4), 'Charge':-1., 'Density':0.001, 'Mass':1., 'Grid':BoxGrid, 
'TimeStep':dt, 'Features':('Rsliced',)}
specie2_in = {'FixedCell':(3,3,4), 'Charge':1., 'Density':0.001, 'Mass':1886., 'Grid':BoxGrid, 
'TimeStep':dt, 'Features':('Rsliced','Still')}

solver = Solver(comm,solver_in)
solver.add_gauss_beam(seed_in)
specie1 = Specie(comm,specie1_in)
specie2 = Specie(comm,specie2_in)

MovingFrame = {'TimeStep':dt,'Steps':20,'AbsorbLayer':7.5,'AddPlasma':densprof}
chimera_in = {'Solvers':(solver,),'Particles':(specie1,specie2,),'MovingFrames':(MovingFrame,)}
Chimera = ChimeraRun(comm,chimera_in)
