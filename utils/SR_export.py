import numpy as np
from scipy.ndimage import gaussian_filter
from tvtk.api import tvtk, write_data

def ExportToVTK(sr_calc_far, filter=None, filename='spectrum_total'):

    omega, theta, phi = sr_calc_far.Args['omega'], sr_calc_far.Args['theta'], sr_calc_far.Args['phi'], 
    
    phi = np.r_[phi, 2*np.pi]
    vals = np.concatenate((sr_calc_far.Data['Rad'], sr_calc_far.Data['Rad'][:,:,[0]]), 
                      axis=-1)

    if filter is not None:
        vals = gaussian_filter(vals, filter)
    
    Nom = omega.size
    Nth = theta.size
    Nph = phi.size
    
    omega = omega[:,None,None]
    theta = theta[None,:, None]
    phi = phi[None,None,:]
    
    x, y, z = (omega*np.sin(theta)*np.sin(phi)).ravel(), \
          (omega*np.sin(theta)*np.cos(phi)).ravel(), \
          (omega*np.cos(theta)*np.ones_like(phi)).ravel()
    
    spc_vtk = tvtk.StructuredGrid(dimensions=(Nph, Nth, Nom), 
                              points=np.vstack((x.ravel(),y.ravel(),z.ravel())).T)
    
    spc_vtk.point_data.scalars = vals.flatten()
    spc_vtk.point_data.scalars.name = 'spectrum'
    
    write_data(spc_vtk, filename)
    
