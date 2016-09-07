import fimera as chimera

comp_ind = {'x':0,'y':1,'z':2,'all':3}

SychRadCalcs = {
'coh'   : chimera.spectrum_comp_coh,
'incoh' : chimera.spectrum_comp_incoh_wgh,
'both'  : chimera.spectrum_comp_both
}
