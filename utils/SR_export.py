import numpy as np
from scipy.ndimage import gaussian_filter
from tvtk.api import tvtk, write_data

def ExportToVTK( sr_calc_far, filter=None, filename='spectrum_total',
                 project=False, lambda0_um=1 ):

    omega, theta, phi = sr_calc_far.Args['omega'], sr_calc_far.Args['theta'], sr_calc_far.Args['phi']
    phi = np.r_[phi, 2*np.pi]

    if project is False:
        vals = sr_calc_far.Data['Rad']
        scalar_name = 'spectrum'
    else:
        vals = sr_calc_far.get_spot( chim_units=False, lambda0_um=lambda0_um )
        vals = vals[None, :, :]
        omega = omega[[-1]]
        filename += '_proj'
        scalar_name = 'spectrum_proj'

    vals = np.concatenate( (vals, vals[:, :, [0]]), axis=-1 )

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
    spc_vtk.point_data.scalars.name = scalar_name
    
    write_data(spc_vtk, filename)
    
