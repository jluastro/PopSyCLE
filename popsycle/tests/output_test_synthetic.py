import warnings
warnings.filterwarnings('ignore')
from popsycle import synthetic
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import h5py

def test_h5_output(ebf_file, reference_h5_file, extra_col= False):
    """"
    Parameters
    ----------
     ebf_file : str
        Name of the ebf file used to generate the reference h5 file
    
     reference_h5_file : str
        Name of the file to compare new output to (should be run with seed=42 on the ebf_file)

     extra_col : boolean, defaults to False
        Tells the code whether or not the new h5 file will have additional columns (ie does the new version of
        popsycle give more information than before

    """
     
    #create the new h5 file by running popsycle
    synthetic.perform_pop_syn(ebf_file = ebf_file,
                              output_root = 'test',
                              iso_dir = '/u/casey/scratch/work/microlens/popsycle_test/isochrones/',
                              bin_edges_number = None, overwrite = True, seed=42);

    #read in the data from the reference h5 file
    hfr = h5py.File(reference_h5_file, 'r')
    ref_dset = np.concatenate((hfr['l0b0'], hfr['l0b1'], hfr['l1b0'], hfr['l1b1']),
                              axis=1)
    hfr.close()
    
    #read in the data from the test h5 file created by popsycle
    hft = h5py.File('test.h5', 'r')
    test_dset = np.concatenate((hft['l0b0'], hft['l0b1'], hft['l1b0'], hft['l1b1']),
                              axis=1)
    hft.close()

    #see if we have the right number of columns
    if test_dset.shape[0] != ref_dset.shape[0]:
        if not extra_col:
            print("the h5 files are not the same size. Run again with extra_col=True if you have added columns")

    #test to see whether the files are the same
    matched_col=0 #initialize matched_col counter
    for i in range(0, ref_dset.shape[0]):
        test_col = test_dset[i,:]
        ref_col = ref_dset[i, :]
        if test_col.all() == ref_col.all():
            matched_col = matched_col+1

        #check to see if disagreements are because of nans
        else:
            bad_idxs = np.where(ref_col != test_col)
            ref_nan_idx = np.where(ref_col == np.nan)
            test_nan_idx = np.where(test_col == np.nan)
            if test_nan_idx.all() == ref_nan_idx.all() and bad_idxs.all() == ref_nan_idx.all():
                matched_col = matched_col+1
            else:
                matched_col= matched_col
                print('Test failed in column', i)

    if matched_col == ref_dset.shape[0]:
        print("The new test h5 file matched the reference file!")

    else:
        print("The new test h5 file does not match the reference file")

    return

    
            

    
         




 







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
                        photometric_system = 'ubv',
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
