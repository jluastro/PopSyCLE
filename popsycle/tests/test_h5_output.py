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
    if test_dset.shape[0] != ref_dset.shape[0] and not extra_col:
            assert test_dset.shape[0] == ref_dset.shape[0], "the h5 files are not the same size. Run again with extra_col=True if you have added columns")

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
                assert test_nan_idx.all() == ref_nan_idx.all(), "Files do not have nan values at the same indices"
                assert bad_idxs.all() == ref_nan_idx.all(), "Coulumns disagree at non-nan values"

    assert matched_col == ref_dset.shape[0], "The new test h5 file does not match the reference file!")

  

    return
