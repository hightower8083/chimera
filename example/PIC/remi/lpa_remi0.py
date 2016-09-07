out_folder = '../lpa_remi/'
dt = 0.04
Steps2Do = 20001
BoxGrid = (-120.,0.0,60.,0.04,0.3)

densprof = np.vectorize( lambda x,y,z: (x>=0.)*(x<=125.)*x/125. + (x>125.)*(x<=375.) + 
(x>375.)*(x<=500.)*(1-0.5*(x-375.)/125.) + 0.5*(x>500.))

solver_in = {'MaxAzimuthMode':1,'Grid':BoxGrid,'TimeStep':dt, 'Features':('Rsliced','SpaceCharge')}

seed_in = {'a0':4.0,'k0':1.0,'x0':-45.0,'x_foc':-45.0,'Lx':12.5,'LR':20.}

specie1_in = {'FixedCell':(2,2,4), 'Charge':-1., 'Density':0.0005734, 'Mass':1., 'Grid':BoxGrid, 
'TimeStep':dt, 'Features':('Rsliced',)}

specie2_in = {'FixedCell':(2,2,4), 'Charge': 1., 'Density':0.0005734, 'Mass':1886., 'Grid':BoxGrid, 
'TimeStep':dt, 'Features':('Rsliced','Still')}

MovingFrame = {'TimeStep':dt,'Steps':20,'AbsorbLayer':10.,'AddPlasma':densprof}
