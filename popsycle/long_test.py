import warnings
warnings.filterwarnings('ignore')
from popsycle import synthetic, utils, binary_utils
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from spisea.imf import multiplicity
import h5py
from astropy import units
synthetic.calc_events(hdf5_file = 'example-multi-clean.h5',
                      hdf5_file_comp = 'example-multi-clean_companions.h5',
                      output_root2 = 'example-multi-clean-new-LONG', 
                      radius_cut = 2, 
                      obs_time = 1000, 
                      theta_frac = 2, 
                      blend_rad = 0.65, 
                      overwrite = False, 
                      n_proc = 1)