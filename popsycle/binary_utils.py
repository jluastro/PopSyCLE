import numpy as np
import warnings
from astropy.table import Table, Column
from ast import literal_eval


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
    
            t_matched_multipeaks = t_comp_rb_mp[in_matched_field][matched_multipeaks[in_matched_field] & pos]
            
            observable_peaks = len(np.where(t_matched_multipeaks['delta_m'] > delta_m_cutoff_secondary)[0]) + 1
            fixed_peaks.append(observable_peaks)
            
    t_prim[observable_n_peaks_col][t_prim[n_peaks_col] > 1] = fixed_peaks

    return t_prim

def cut_Mruns(t_prim, t_comp_rb, t_comp_rb_mp, min_mag, delta_m_cut, u0_cut, ubv_filter, S_LSN):
    """
    Make observational cuts on PopSyCLE runs with multiple systems

    Parameters
    ----------
    t_prim : Astropy table
        Events table from refine_binary_events.
        Must contain 'observable_n_peaks' column.

    t_comp_rb : Astropy table
        Companion table from refine_binary_events.

    t_comp_rb_mp : Astropy table
        Multi peak table from refine binary events 
        (each row corresponds to a peak in a lightcurve).

    min_mag : float
        Minimum baseline or source magnitude (specified by S_LSN).

    delta_m_cut : float or None.
        Minimum bump magnitude.

    u0_cut : float
        Maximum u0.

    ubv_filter : str
        Filter name used when cutting on min_mag and delta_m_cut.

    S_LSN : str
        'S' for source mag cut or 'LSN' for baseline magnitude cut.

    Returns
    -------
    t_both_mcut : Astropy table
        Table with specified observational cuts.
        
    t_both_mcut_one_peak : Astropy table
        Table with specified observational cuts and only single peaked events.
        
    t_multiples_mcut_multi_peak : Astropy table
        Table with specified observational cuts and only multipeaked events
        containing a multiple system.
    """
    #S_LSN is source or baseline mag cut
    if S_LSN == 'S':
        mag_cut = t_prim['ubv_{}_app_S'.format(ubv_filter)] <= min_mag
    elif S_LSN == 'LSN':
        mag_cut = t_prim['ubv_{}_app_LSN'.format(ubv_filter)] <= min_mag
    
    u0_cut = np.abs(t_prim['u0']) < u0_cut
    
    binary_filt = (t_prim['isMultiple_L'] == 1) | (t_prim['isMultiple_S'] == 1)
    single_filt = (t_prim['isMultiple_L'] == 0) & (t_prim['isMultiple_S'] == 0)
    assert(len(t_prim) == (sum(binary_filt) + sum(single_filt)))

    if delta_m_cut is not None:
        delta_m_cut = ((t_prim['bin_delta_m'] > 0.1) & binary_filt) | ((t_prim['delta_m_{}'.format(ubv_filter)] > 0.1) & single_filt)
        total_cut = mag_cut & u0_cut & delta_m_cut
    else:
        total_cut = mag_cut & u0_cut

    t_both_mcut = t_prim[total_cut]
    binary_filt_cut = (t_both_mcut['isMultiple_L'] == 1) | (t_both_mcut['isMultiple_S'] == 1)
    single_filt_cut = (t_both_mcut['isMultiple_L'] == 0) & (t_both_mcut['isMultiple_S'] == 0)
    assert(len(t_both_mcut) == (sum(binary_filt_cut) + sum(single_filt_cut)))

    t_mult_mcut_no_peaks = t_both_mcut[binary_filt_cut & (t_both_mcut['observable_n_peaks'] == 0)] 
    t_both_mcut_one_peak = t_both_mcut[(binary_filt_cut & (t_both_mcut['observable_n_peaks'] == 1)) | single_filt_cut]
    t_multiples_mcut_multi_peak = t_both_mcut[binary_filt_cut & (t_both_mcut['observable_n_peaks'] > 1)]
    assert(len(t_both_mcut) == (len(t_mult_mcut_no_peaks) + len(t_both_mcut_one_peak) + len(t_multiples_mcut_multi_peak)))

    return t_both_mcut, t_both_mcut_one_peak, t_multiples_mcut_multi_peak





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

    
