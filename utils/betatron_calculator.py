from time import time
import numpy as np
import sys

from opmd_viewer import OpenPMDTimeSeries
from opmd_viewer import ParticleTracker

from chimera.moduls.SR_ext import SR
from betatron_calculator_methods import get_orbits



data_path = '/media/Data/RUNS/KTP/straight/density_scan/fbpic_n4e19_a1.92/'
name_base = 'betatron_spectrum_3D'

sys.path.append(data_path)
from lwfa_script import diag_period_often as diag_int_pic
from lwfa_script import dt as dt_pic
from lwfa_script import c

dt_pic *= c*1e6     # convert from fs to um

reference_iteration = 55000
dStep = 1
dNpart = 1
Nw, Nth, Nph = 400, 24, 24

T0 = time()

dt = dt_pic * diag_int_pic * dStep
sr_in_far = {'Grid':[(1., 6e4),(0., 0.12),(0, 2*np.pi),(Nw, Nth, Nph)],
             'TimeStep':dt, 'Features':(), }

sr_calc_far = SR(sr_in_far)

if __name__ == '__main__':

    ts = OpenPMDTimeSeries(data_path+'diags/hdf5/', check_all_files=False)
    pt = ParticleTracker(ts, iteration=reference_iteration, species='electrons',
                         select={'uz':[95, None]}, preserve_particle_index=True)

    sr_calc_far = get_orbits(sr_calc_far, ts, pt, dStep, dNpart)
    sr_calc_far.calculate_spectrum(comp='all')

    for key in ['coords', 'momenta_prv',  'momenta_nxt', 'weights']:
        sr_calc_far.Data[key] = None

    name = name_base + '_dT_' + str(dStep) + '_dNp_' + str(dNpart) + '.npy'
    np.save(data_path + name, sr_calc_far.Data)

    print('done in {:g} minutes'.format( (time()-T0)/60.) )
