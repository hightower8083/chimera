import numpy as np
from numba import jit

@jit
def check_orbit_edges(up, un, w):
    for i in range(w.shape[0]):
        for p in range(w.shape[1]):
            if up[0,i,p] == 0.:
                up[:,i,p] = un[:,i,p]
            elif un[0,i,p] == 0.:
                un[:,i,p] = up[:,i,p]

def get_orbits(sr_calc_far, ts, pt, dStep, dNpart):
    rec_indx_start = 0
    Steps = len(ts.iterations) - rec_indx_start - 1
    Steps = Steps//dStep
    Nparts = int(np.ceil(pt.N_selected/dNpart))

    sr_calc_far.Data['coords'] = np.zeros((3, Steps, Nparts), order='F')
    sr_calc_far.Data['momenta_nxt'] = np.zeros((3, Steps, Nparts), order='F')
    sr_calc_far.Data['momenta_prv'] = np.zeros((3, Steps, Nparts), order='F')
    sr_calc_far.Data['weights'] = np.zeros((Steps, Nparts), order='F')

    for i in range(1, Steps):
        x0, y0, z0, ux0, uy0, uz0, w0 = ts.get_particle(['x','y','z','ux','uy','uz','w'], select=pt,
                                                        iteration=ts.iterations[rec_indx_start+i*dStep],
                                                        species='electrons' )
        if x0.size>1:
            x1, y1, z1 = ts.get_particle(['x','y','z',], select=pt,
                                         iteration=ts.iterations[rec_indx_start+(i+1)*dStep],
                                         species='electrons' )
            if x1.size>1:
                x0, y0, z0 = 0.5*(x0[::dNpart]+x1[::dNpart]), \
                             0.5*(y0[::dNpart]+y1[::dNpart]), \
                             0.5*(z0[::dNpart]+z1[::dNpart])
            else:
                x0, y0, z0  = x0[::dNpart], y0[::dNpart], z0[::dNpart]

            ux0, uy0, uz0, w0 = ux0[::dNpart], uy0[::dNpart], uz0[::dNpart], w0[::dNpart]
            sr_calc_far.Data['coords'][:, i, :] = np.array([z0, x0, y0])
            sr_calc_far.Data['momenta_prv'][:, i, :] =  np.array([uz0, ux0, uy0])
            sr_calc_far.Data['momenta_nxt'][:, i-1, :] = np.array([uz0, ux0, uy0])
            sr_calc_far.Data['weights'][i, :] = w0

    for key in ['coords', 'momenta_prv',  'momenta_nxt', 'weights']:
        sr_calc_far.Data[key] = np.nan_to_num(sr_calc_far.Data[key])

    check_orbit_edges(sr_calc_far.Data['momenta_prv'],
                      sr_calc_far.Data['momenta_nxt'],
                      sr_calc_far.Data['weights'])
    return sr_calc_far
