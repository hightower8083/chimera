from mayavi import mlab
import numpy as np
from scipy.interpolate import griddata
from input import *

b = np.load(out_folder+'data_spect3D_incoh_all_1.npy')[:,:,:,-1]
bins_w, bins_x, bins_y = 100, 70, 70

om, theta, phi = np.mgrid[0:omega_max:b.shape[0]*1j,0:theta_max:b.shape[1]*1j, 0:2*np.pi:b.shape[2]*1j]
x,y = np.cos(phi)*np.sin(theta), np.sin(phi)*np.sin(theta)
oo, xx,yy = np.mgrid[ 0:omega_max:bins_w*1j ,x.min():x.max():bins_x*1j, y.min():y.max():bins_y*1j]
bbb = griddata( np.vstack((om.flatten(), x.flatten(), y.flatten())).T, b.flatten(), (oo,xx,yy),fill_value=0.0,method='nearest')

fig = mlab.figure(size=(800,800),bgcolor=(1,1,1))
s = mlab.pipeline.scalar_field(bbb)
v_max = bbb.max()
mlab.pipeline.volume(s,vmin=0.4*v_max,vmax=0.8*v_max)

