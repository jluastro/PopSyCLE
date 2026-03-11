import numpy as np
import warnings
from astropy.table import Table, Column
from ast import literal_eval
import h5py
import pandas as pd
import os
from popsycle import synthetic


def add_magnitudes(mags):
    """
    Adds a list of magnitudes
    
    Parameters
    ----------
    mags: array-like
        List or array of magnitudes
    
    Returns
    -------
    m_sum : float
        Sum of input magnitudes
    """
    mags = np.array(mags)
    fluxes = 10**(-0.4*mags)
    fluxes = np.nan_to_num(fluxes)
    
    if np.all(fluxes == 0):
        m_sum = np.nan
    else:
        m_sum = -2.5*np.log10(np.sum(fluxes, axis = 0))
    
    return m_sum

def subtract_magnitudes(m1, m2):
    """
    Subtracts two magnitudes
    Expects m1 to be brighter than m2.
    Will return nan if the opposite is the case
    or if m1 is nan.
    
    Parameters
    ----------
    m1 : float or array-like
        Brighter magnitude
    m2 : float or array-like
        Dimmer magnitude
    
    Returns
    -------
    m_diff : float or array-like
        Difference in the corresponding flux of 
        m1 - m2
    """
    
    m1 = np.array(m1)
    m2 = np.array(m2)
        
    f1 = 10**(-0.4*m1)
    f2 = 10**(-0.4*m2)
        
    f1 = np.nan_to_num(f1)
    f2 = np.nan_to_num(f2)
    
    m_diff = -2.5*np.log10(f1 - f2)
    
    diff_criteria =  np.abs(f1 - f2) < 1e-12
    if np.shape(diff_criteria) == ():
        if diff_criteria == True:
            m_diff = np.nan
    else:
        m_diff[diff_criteria] = np.nan
    
    if np.any(f1 - f2 < 0):
        warnings.warn("Warning: first magnitude dimmer than second magnitude. Result will be nan.")
    
    return m_diff

def primary_mag_from_system_mag(system_mag, companion_mag):
    return subtract_magnitudes(system_mag, companion_mag)

def event_table_companion_idxs_to_lists(events):
    """
    The event table from refine_binary_events() is saved with
    'companion_idx_list' as strs so it can be saved as a .fits file.
    This function switches that column to lists.
    (Note this changes in place, so input table is changed). 
    
    Parameters
    ----------
    events : Astropy Table
        Event table from refine_binary_events with 'companion_idx_list' as strs
    
    Returns
    -------
    events : Astropy Table
        Event table with 'companion_idx_list' as lists
    """
    
    lists = []
    for i in events['companion_idx_list'].astype('str'):
        lists.append(literal_eval(i))
        
    events['companion_idx_list'] = lists
    
    return events

def add_observable_peaks_column(t_prim, t_comp_rb, t_comp_rb_mp, t_lightcurves, delta_m_cutoff = 0.1, delta_m_cutoff_secondary = 0.02,
                               match_by_field = True):
    """
    Adds new column to event table with the number of observable peaks

    Parameters
    ----------
    t_prim : Astropy table
        Events table from refine_binary_events.

    t_comp_rb : Astropy table
        Companion table from refine_binary_events.

    t_comp_rb_mp : Astropy table
        Multi peak table from refine binary events 
        (each row corresponds to a peak in a lightcurve).

    t_lightcurves : Astropy table
        Lightcurve table from refine_binary_events
        (each row corresponds to a generated lightcurve).

    delta_m_cutoff : float, optional
        Bump magnitude cutoff for primary peak 
        (if primary peak < delta_m_cutoff, the number of observable peaks is 0).
        Default is 0.1.

    delta_m_cutoff_secondary : float, optional
        Bump magnitude cutoff for secondary peaks. 
        Default is 0.02.

    match_by_field : bool, optional
        If there are multiple fields, looks for a column called "field_id"
        and requires the companions and primaries to be in the same field.
        Default is True.

    Returns
    -------
    t_prim : Astropy table
        Events table with additional column with the number of observable peaks.
    """
    
    fixed_peaks = []
    multiplt_comp_rbe_comp_same_SL_pair_counter = 0
    idx_multiple_comp_same_SL_pair_counter = 0
    
    n_peaks_col = 'n_peaks'
    bin_delta_m_col = 'bin_delta_m'
    observable_n_peaks_col = 'observable_n_peaks'

    try:
        t_prim.add_column( Column(np.zeros(len(t_prim), dtype=float), name=observable_n_peaks_col) )
    except ValueError:
        t_prim.replace_column(observable_n_peaks_col, np.zeros(len(t_prim)))
    t_prim[observable_n_peaks_col][:] = np.nan
    # If 1 or 0 peaks, append those to observable_n_peaks
    # if the single peak is less than a delta_m threshold, 
    # it will be cut later by bin_delta_m
    t_prim[observable_n_peaks_col][t_prim[n_peaks_col] <= 1] = t_prim[t_prim[n_peaks_col] <= 1][n_peaks_col]

    t_prim_multi_peaks = t_prim[t_prim[n_peaks_col] > 1]

    # only one companion per position in multipeak table. Keeps track for speed
    total_pos = np.full(len(t_comp_rb_mp['companion_idx']), False)
    
    for i in range(len(t_prim_multi_peaks)):
        matched_multipeaks = (t_comp_rb_mp['obj_id_L'] == t_prim_multi_peaks['obj_id_L'][i]) & (t_comp_rb_mp['obj_id_S'] == t_prim_multi_peaks['obj_id_S'][i])
        #matched_multipeaks = np.where(np.logical_and(t_comp_rb_mp['obj_id_L'] == t_prim['obj_id_L'][i], t_comp_rb_mp['obj_id_S'] == t_prim['obj_id_S'][i]))[0]
        
        # grabs companion list for the used lightcurve (max binary_delta_m as determined in refine_binary_events)
        # there would be more than one companion to choose between if there are any triples
        matched_lightcurves = (t_lightcurves['obj_id_L'] == t_prim_multi_peaks['obj_id_L'][i]) & (t_lightcurves['obj_id_S'] == t_prim_multi_peaks['obj_id_S'][i])
        matched_multipeaks_companion_idx_L = t_lightcurves['companion_id_L'][matched_lightcurves & (t_lightcurves['used_lightcurve'] == True)].value.data[0]
        matched_multipeaks_companion_idx_S = t_lightcurves['companion_id_S'][matched_lightcurves & (t_lightcurves['used_lightcurve'] == True)].value.data[0]
        
        if match_by_field:
            matched_multipeaks_field = t_lightcurves['field_id'][matched_lightcurves & (t_lightcurves['used_lightcurve'] == True)].value[0]
            
        if t_lightcurves['class'][i] == 'BSPL':
            matched_multipeaks_companions_list = [matched_multipeaks_companion_idx_S]
        elif t_lightcurves['class'][i] == 'PSBL':
            matched_multipeaks_companions_list = [matched_multipeaks_companion_idx_L]
        elif t_lightcurves['class'][i] == 'BSBL':
            matched_multipeaks_companions_list = [matched_multipeaks_companion_idx_L, matched_multipeaks_companion_idx_S]

        #t_matched_multipeaks = t_comp_rb_mp[matched_multipeaks & (t_comp_rb_mp['companion_idx'] == matched_multipeaks_companions_list)]

        # brightest peak not included in mp table, so checks if it's bright enough and if so
        # adds an additional peak to observable_peaks
        if t_prim_multi_peaks[bin_delta_m_col][i] < delta_m_cutoff:
            #t_prim[observable_n_peaks_col][t_prim[n_peaks_col] > 1][i] = 0 
            fixed_peaks.append(0)

        else:
            if match_by_field:
                in_matched_field = t_comp_rb_mp['field_id'] == matched_multipeaks_field
                iterate_companion_idxs = t_comp_rb_mp['companion_idx'][in_matched_field]
            else:
                iterate_companion_idxs = t_comp_rb_mp['companion_idx']
            pos = np.full(len(iterate_companion_idxs), False)
            for j in range(len(iterate_companion_idxs)):
                # only one companion per position in multipeak table
                # skips over matched ones for speed
                if total_pos[j] == True:
                    continue
                current_pos = np.array_equiv(iterate_companion_idxs.data[j], matched_multipeaks_companions_list)
                pos[j] = current_pos
                total_pos[j] = current_pos
                # Assumes all instances are adjacent and kills loop if you've made it to them
                if current_pos == False and sum(pos) > 0:
                    break
    
            if match_by_field:
                t_matched_multipeaks = t_comp_rb_mp[in_matched_field][matched_multipeaks[in_matched_field] & pos]
            else:
                t_matched_multipeaks = t_comp_rb_mp[pos]
            
            observable_peaks = len(np.where(t_matched_multipeaks['delta_m'] > delta_m_cutoff_secondary)[0]) + 1
            fixed_peaks.append(observable_peaks)
            
    t_prim[observable_n_peaks_col][t_prim[n_peaks_col] > 1] = fixed_peaks

    return t_prim



def make_bhs_single(hdf5_file, hdf5_comp_file, bh_binary_frac = 0.1, phots = ['ubv_I', 'ubv_K', 'ubv_J', 'ubv_U', 'ubv_R', 'ubv_B', 'ubv_V', 'ubv_H'],
                    new_hdf5_file = None, new_hdf5_file_comp = None, symlink_aux_files = True):
    """
    This makes some fraction of BHs singles.
    Currently no binary star evolution, so all BHs end up in binaries.
    We drop the companions from the companion table and set the BH parameters
    to those of a single BH.
    These are saved to a new file.

    To be run after perform_pop_syn.

    Parameters
    ----------
    hdf5_file : str
        File name of hdf5 file to modify.

    hdf5_comp_file : str
        File name of hdf5 companion file to modify.

    bh_binary_frac : float
        Binary fraction to be kept of BHs.
        Default is 0.1.

    phots : list of str
        Photometry values to be made nan for BHs.
        Should include all photometry values in table.
        Default is all those in ubv system.

    new_hdf5_file : str or None
        New hdf5 file name.
        Default is None which saves it as 
        hdf5_file[:-3] + '_{}_bhb_frac.h5'.format(bh_binary_frac).
        
    new_hdf5_file_comp : str or None
        New hdf5 file name.
        Default is None which saves it as
        hdf5_comp_file[:-3] + '_{}_bhb_frac.h5'.format(bh_binary_frac).

    symlink_aux_files : bool
        Makes symbolic links to the following necessary auxiliary files with the new root:
        _perform_pop_syn.log
        _galaxia.log
        _galaxia_params.txt
        Default is True.
    """
    
    tmp_prim = h5py.File(hdf5_file, 'r')
    keys = tmp_prim.keys()  
    tmp_comp = h5py.File(hdf5_comp_file, 'r')
    keys_comp = tmp_comp.keys()

    if new_hdf5_file is None:
        new_hdf5_file = hdf5_file[:-3] + '_{}_bhb_frac.h5'.format(bh_binary_frac)
    if new_hdf5_file_comp is None:
        new_hdf5_file_comp = hdf5_comp_file[:-13] + '{}_bhb_frac_companions.h5'.format(bh_binary_frac)
    if symlink_aux_files:
        os.symlink(hdf5_file[:-3] + '_galaxia.log', new_hdf5_file[:-3] + '_galaxia.log')
        os.symlink(hdf5_file[:-3] + '_galaxia_params.txt', new_hdf5_file[:-3] + '_galaxia_params.txt')
        os.symlink(hdf5_file[:-3] + '_perform_pop_syn.log', new_hdf5_file[:-3] + '_perform_pop_syn.log')

    prim_copy = h5py.File(new_hdf5_file, 'w')
    prim_copy[list(keys)[-2]] = tmp_prim[list(keys)[-2]][:]
    prim_copy[list(keys)[-1]] = tmp_prim[list(keys)[-1]][:]
    prim_copy.close()
    
    comp_copy = h5py.File(new_hdf5_file_comp, 'w')
    comp_copy[list(keys)[-2]] = tmp_comp[list(keys)[-2]][:]
    comp_copy[list(keys)[-1]] = tmp_comp[list(keys)[-1]][:]
    comp_copy.close()

    del tmp_prim
    del tmp_comp

    for i in list(keys)[1:-2]:
        if i[0] == 'l':
            prim = pd.read_hdf(hdf5_file, i).set_index(['obj_id'])
            bbh_prim_crit = (prim['rem_id'] == 103) & (prim['isMultiple'] == 1)
            bh_prim = prim[bbh_prim_crit]
            #idxs of bhs that will be made single
            bh_prim_singlify_idxs = np.random.choice(bh_prim.index, size = int(len(bh_prim)*(1-bh_binary_frac)), replace = False)
            
            bh_prim.loc[bh_prim_singlify_idxs, ('isMultiple')] = 0
            bh_prim.loc[bh_prim_singlify_idxs, ('N_companions')] = 0
            bh_prim.loc[bh_prim_singlify_idxs, ('systemMass')] = bh_prim.loc[bh_prim_singlify_idxs, ('mass')]
            # set photometry to nan
            for phot in phots:
                bh_prim.loc[bh_prim_singlify_idxs, (phot)] = np.nan
    
            prim[bbh_prim_crit] = bh_prim
            
            comp = pd.read_hdf(hdf5_comp_file, i).set_index(['system_idx'])
            comp.drop(index=bh_prim_singlify_idxs, axis=0, inplace=True)

            # Verify number of companions in table same as accounted for in primary table
            assert(len(comp) == np.sum(prim['N_companions']))

            prim.reset_index(inplace=True)
            comp.reset_index(inplace=True)

            prim_hdf5 = h5py.File(new_hdf5_file, 'r+')
            compound_dtype = synthetic._generate_compound_dtype(prim.dtypes.to_dict())
            save_data = np.empty(len(prim), dtype=compound_dtype)
            for colname in prim.keys():
                save_data[colname] = prim[colname].to_numpy()
            dataset = prim_hdf5.create_dataset(i, shape=(0,),
                                        chunks=(1e4,),
                                        maxshape=(None,),
                                        dtype=compound_dtype)
            dataset.resize((len(prim),))
            prim_hdf5[i][:] = save_data
            prim_hdf5.close()

            del prim, save_data

            comp_hdf5 = h5py.File(new_hdf5_file_comp, 'r+')
            compound_dtype = synthetic._generate_compound_dtype(comp.dtypes.to_dict())
            save_data = np.empty(len(comp), dtype=compound_dtype)
            for colname in comp.keys():
                save_data[colname] = comp[colname].to_numpy()
            dataset = comp_hdf5.create_dataset(i, shape=(0,),
                                        chunks=(1e4,),
                                        maxshape=(None,),
                                        dtype=compound_dtype)
            dataset.resize((len(comp),))
            comp_hdf5[i][:] = save_data
            comp_hdf5.close()

            del comp, save_data
            
                #prim_hdf5.create_dataset(i)#, data=prim_np.astype("|V256"))

            #with h5py.File(new_hdf5_file_comp, 'r+') as comp_hdf5:
            #    comp_np = comp.reset_index().to_numpy()
            #    comp_hdf5.create_dataset(i)#, data=comp_np.astype("|V256"))
                
    return




#def primary_mag_from_system_mag_hdf5(prim_hdf5, comp_hdf5, mag_colnames):
#    grouped_companions = companions_table.group_by(['system_idx'])
#    for mag_col in mag_colnames:
#        companions_mag = grouped_companions[mag_col].groups.aggregate(add_magnitudes)
#        prim_table[mag_col + '_prim'] = subtract_magnitudes(prim_table[mag_col], companions_mag)
#    return prim_table

#def primary_mag_from_system_mag_Table(prim_table, companions_table, mag_colnames):
    
#    prim_df = prim_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'])
#    companion_df = companion_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'])
#    joined = companion_df.join(prim_df, lsuffix='_comp', rsuffix='_prim', how='outer')
#    
#    grouped_companions = companions_table.group_by(['system_idx'])
#    for mag_col in mag_colnames:
#        companions_mag = grouped_companions[mag_col].groups.aggregate(add_magnitudes)
#        prim_table[mag_col + '_prim'] = subtract_magnitudes(prim_table[mag_col], companions_mag)
#    return prim_table

    
