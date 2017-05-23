import numpy as np
from scipy.constants import m_e,c,e

mc2_GeV = m_e*c**2/e*1e-9

def ocelot_to_chimera(p_arrays,beam,lam0,keep_orig=True,\
  monochrom=None, select_parts = None):
	"""
	Exports the list of ParticleArrays from OCELOT to CHIMERA

	Parameters
	----------
	p_arrays: list of oclt.ParticleArray objects
	  beam represented in the form of slices list
	beam : chimera.Species object
	  CHIMERA species obejct to populate with particles
	lam0: float
	  normalization length for CHIMERA in meters

	Returns
	-------
	beam : chimera.Species object
		CHIMERA species obejct populated with particles
	"""

	Np = np.int(np.sum([p_array.size() for p_array in p_arrays]))

	xx = np.hstack(([(p_array.s-p_array.tau())/lam0 for p_array in p_arrays]))
	yy = np.hstack(([p_array.y()/lam0 for p_array in p_arrays]))
	zz = np.hstack(([p_array.x()/lam0 for p_array in p_arrays]))

	gg = np.hstack(([(p_array.p()+1)*p_array.E/mc2_GeV for p_array in p_arrays]))
	oy = np.hstack(([p_array.py() for p_array in p_arrays]))
	oz = np.hstack(([p_array.px() for p_array in p_arrays]))

	qq = np.hstack(([p_array.q_array*1e12 for p_array in p_arrays]))

	px = np.sqrt( (gg**2-1.)/(1+oy**2+oz**2) )
	py = px*oy
	pz = px*oz
	qq = -qq/(beam.weight2pC*lam0*1e6)

	if monochrom is not None:
		emin,emax = monochrom
		indx = np.nonzero( (gg*mc2_GeV<emax)*(gg*mc2_GeV>emin)  )[0]
		xx = xx[indx]
		yy = yy[indx]
		zz = zz[indx]
		px = px[indx]
		py = py[indx]
		pz = pz[indx]
		qq = qq[indx]
		Np = indx.shape[0]

	if select_parts is not None:
		indx = np.arange(xx.shape[0])
		np.random.shuffle(indx)
		indx = indx[:select_parts]
		xx = xx[indx]
		yy = yy[indx]
		zz = zz[indx]
		px = px[indx]
		py = py[indx]
		pz = pz[indx]
		qq = qq[indx]*Np*1.0/select_parts
		Np = indx.shape[0]

	gg, oy, oz = None, None, None

	beam.Data['coords'] = np.zeros((3,Np))
	beam.Data['momenta'] = np.zeros((3,Np))
	beam.Data['weights'] = np.zeros(Np)

	beam.Data['coords'][0] = xx
	beam.Data['coords'][1] = yy
	beam.Data['coords'][2] = zz

	beam.Data['momenta'][0] = px
	beam.Data['momenta'][1] = py
	beam.Data['momenta'][2] = pz

	beam.Data['weights'][:] = qq

	beam.Data['coords'][0] -= (beam.Data['coords'][0]*beam.Data['weights']\
	  ).sum()/(beam.Data['weights']).sum()
	beam.Data['coords_halfstep'] = beam.Data['coords'].copy()
	if keep_orig is False:
		del p_arrays

	return  beam
