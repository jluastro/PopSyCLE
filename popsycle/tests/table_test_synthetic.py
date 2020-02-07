import warnings
warnings.filterwarnings('ignore')
from popsycle import synthetic
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import h5py


synthetic.perform_pop_syn(ebf_file = 'example.ebf',   
                          output_root = 'test',  
                          iso_dir = '/u/casey/scratch/work/microlens/popsycle_test/isochrones/',  
                          bin_edges_number = None, overwrite = True, seed=42);

synthetic.calc_events(hdf5_file = 'test.h5', 
                      output_root2 = 'test', 
                      radius_cut = 2, 
                      obs_time = 1000, 
                      n_obs = 11, 
                      theta_frac = 2, 
                      blend_rad = 0.65, 
                      overwrite = True, 
                      n_proc = 1, seed=42)

synthetic.refine_events(input_root = 'test', 
                        filter_name = 'i', 
                        red_law = 'Damineli16', 
                        overwrite = True, 
                        output_file = 'default')


hfr = h5py.File('trial_1.h5', 'r')
print((list(hfr.keys())))
dsetr = hfr['l0b0']
print(dsetr)
hfr.close()

hft = h5py.File('test.h5', 'r+')
print(list(hft.keys()))
dsett = hft['l0b0']
print(dsett)
dsett.resize((27, 176660))
print(dsett)
print('test==trial', dsetr == dsett)
hft.close()

tabr = Table.read('trial_1_refined_events_i_Damineli16.fits')
tabt = Table.read('test_refined_events_i_Damineli16.fits')
print('The col names of the trial 1 fits are:', tabr.colnames)
print('The col names of the test fits are:', tabt.colnames)
tabt.remove_columns(['teff_S', 'grav_S', 'lum_S', 'teff_L', 'grav_L', 'lum_L'])
print('Check if tables differ:', tabr == tabt)
