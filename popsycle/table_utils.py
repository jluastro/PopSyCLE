import pandas as pd
import numpy as np
from astropy.table import Table, vstack, Column
import os

def combine_re_popsycle_tables(base_dir, filter_str, red_law, modifier_str = '', with_multiples = False):
    """
    Combine together tables from different fields in PopSyCLE after refine events
    and adds a field id column.

    Inputs
    ------
    base_dir : str
        Base directory of field folders.
        (example : /g3/PopSyCLE_sims/roman_v2025/).

    filter_str : str
        Filter string (i.e. 'ubv_I', 'multi_filt')

    red_law : str
        Reddening law (i.e. Damineli16)

    modifier_str : str, Optional
        String after field name (i.e. '_0.1_bhb_frac').
        Default is ''.

    with_multiples : bool, Optional
        Also combines companions table if True.

    Returns
    -------
    events : Astropy table
        Combined PopSyCLE event table.

    companions : Astropy table, Optional
        Combined PopSyCLE companions table.
        Only returns if with_multiples == True.
    
    """
    fields = os.listdir(base_dir)
    
    event_tables = []
    
    if with_multiples:
        companion_tables = []
    for field in fields:
        try:
            base_dir_field = base_dir + '{}/{}'.format(field, field) + modifier_str
            events = Table.read(base_dir_field + '_refined_events_' + filter_str + '_' + red_law + '.fits')
            events.add_column( Column((np.repeat(field, len(events))), name='field_id') )
            event_tables.append(events)

            if with_multiples:
                companions = Table.read(base_dir_field + '_refined_events_'  + filter_str + '_' + red_law + '_companions.fits')
                companions.add_column( Column((np.repeat(field, len(companions))), name='field_id') )
                companion_tables.append(companions)
                
        except FileNotFoundError:
            print('Files not found for field folder: ', field)
            pass
        
    events = vstack(event_tables)
    if with_multiples:
        companions = vstack(companion_tables)
        return events, companions
    else:
        return events

def combine_rbe_popsycle_tables(base_dir, filter_str, red_law, modifier_str = ''):
    """
    Combine together tables from different fields in PopSyCLE after refine binary events
    and adds a field id column.

    Inputs
    ------
    base_dir : str
        Base directory of field folders.
        (example : /g3/PopSyCLE_sims/roman_v2025/).

    filter_str : str
        Filter string (i.e. 'ubv_I', 'multi_filt')

    red_law : str
        Reddening law (i.e. Damineli16)

    modifier_str : str, Optional
        String after field name (i.e. '_0.1_bhb_frac').
        Default is ''.

    Returns
    -------
    events : Astropy table
        Combined PopSyCLE event table.

    companions : Astropy table
        Combined PopSyCLE companions table.

    multi_peak : Astropy table
        Combined PopSyCLE table of the multi peaks.

    lcs : Astropy table
        Combined PopSyCLE table of lightcurves.
    
    """
    fields = os.listdir(base_dir)
    
    event_tables = []
    companion_tables = []
    multi_peak_tables = []
    lc_tables = []
    for field in fields:
        try:
            base_dir_field = base_dir + '{}/{}'.format(field, field) + modifier_str
            events = Table.read(base_dir_field + '_refined_events_' + filter_str + '_' + red_law + '_rb.fits')
            companions = Table.read(base_dir_field + '_refined_events_'  + filter_str + '_' + red_law + '_companions_rb.fits')
            multi_peak = Table.read(base_dir_field + '_refined_events_'  + filter_str + '_' + red_law + '_companions_rb_multi_peaks.fits')
            lcs = Table.read(base_dir_field + '_refined_events_'  + filter_str + '_' + red_law +  '_rb_lightcurves.fits')
        
            events.add_column( Column((np.repeat(field, len(events))), name='field_id') )
            companions.add_column( Column((np.repeat(field, len(companions))), name='field_id') )
            multi_peak.add_column( Column((np.repeat(field, len(multi_peak))), name='field_id') )
            lcs.add_column( Column((np.repeat(field, len(lcs))), name='field_id') )
    
            event_tables.append(events)
            companion_tables.append(companions)
            multi_peak_tables.append(multi_peak)
            lc_tables.append(lcs)
        except FileNotFoundError:
            print('Files not found for field folder: ', field)
            pass
        
    events = vstack(event_tables)
    companions = vstack(companion_tables)
    multi_peak = vstack(multi_peak_tables)
    lcs = vstack(lc_tables)

    return events, companions, multi_peak, lcs

def trim_popsycle_tables(idx, event_tab, comps_tab, lcurv_tab):
    """
    Trim down PopSyCLE event, companions, and lightcurve tables
    to contain only those with the specified indices=idx.
    
    The companions and lightcurve tables are matched against the
    trimmed event table's obj_id_L + obj_id_S.

    Inputs
    ------
    idx : array, ints
        indices returned from a where statement.
    event_tab : astropy.table
        PopSyCLE events table.
    comps_tab : astropy.table
        PopSyCLE companions table.
    lcurve_tab : astropy.table
        PopSyCLE lightcurves table.
    
    This code was generated by Gemini AI and edited by J. Lu.
    """
    import pandas as pd
    import numpy as np

    # --- Step 1: Convert to pandas DataFrames and apply initial trim ---
    # Convert the tables to pandas DataFrames for faster operations
    event_df = event_tab.to_pandas()
    comps_df = comps_tab.to_pandas()
    lcurv_df = lcurv_tab.to_pandas()

    # Apply the initial index trim
    event_df_t = event_df.iloc[idx].copy() # Use iloc for index-based selection

    # --- Step 2: Create a DataFrame of the keys to search for ---
    # This is the list of unique (obj_id_L, obj_id_S) pairs we are interested in
    search_keys = event_df_t[['obj_id_L', 'obj_id_S']].drop_duplicates().copy()

    # --- Step 3: Vectorized filtering using merge or isin() (Recommended) ---
    # The most efficient way is to create a combined key and use isin()
    # Create a combined key column for comparison (e.g., 'key_L_S')
    # Convert columns to string and concatenate them
    search_keys['key_L_S'] = search_keys['obj_id_L'].astype(str) + '_' + search_keys['obj_id_S'].astype(str)
    comps_df['key_L_S'] = comps_df['obj_id_L'].astype(int).astype(str) + '_' + comps_df['obj_id_S'].astype(int).astype(str)
    lcurv_df['key_L_S'] = lcurv_df['obj_id_L'].astype(int).astype(str) + '_' + lcurv_df['obj_id_S'].astype(int).astype(str)

    # Filter the tables using the isin() method on the combined key
    comps_df_t = comps_df[comps_df['key_L_S'].isin(search_keys['key_L_S'])]
    lcurv_df_t = lcurv_df[lcurv_df['key_L_S'].isin(search_keys['key_L_S'])]
    print(f'Trimmed to N_events = {len(event_df_t)}, N_comp = {len(comps_df_t)}, N_lcurve = {len(lcurv_df_t)}')

    # --- Step 4 (Optional): Clean up and convert back to Astropy Table if needed ---
    # Remove the temporary key column
    comps_df_t = comps_df_t.drop(columns=['key_L_S'])
    lcurv_df_t = lcurv_df_t.drop(columns=['key_L_S'])

    # If the result must be an Astropy Table:
    # from astropy.table import Table
    event_tab_t = Table.from_pandas(event_df_t)
    comps_tab_t = Table.from_pandas(comps_df_t)
    lcurv_tab_t = Table.from_pandas(lcurv_df_t)
    
    return event_tab_t, comps_tab_t, lcurv_tab_t


def cut_Mruns(t_prim, t_comp_rb, t_comp_rb_mp, min_mag, delta_m_cut, u0_cut, photometric_system, filter_name, S_LSN):
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

    photometric_system : str
        Photometric system when cutting on min_mag.

    filter_name : str
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
        mag_cut = t_prim['{}_{}_app_S'.format(photometric_system, filter_name)] <= min_mag
    elif S_LSN == 'LSN':
        mag_cut = t_prim['{}_{}_app_LSN'.format(photometric_system, filter_name)] <= min_mag
    
    u0_cut = np.abs(t_prim['u0']) < u0_cut
    
    binary_filt = (t_prim['isMultiple_L'] == 1) | (t_prim['isMultiple_S'] == 1)
    single_filt = (t_prim['isMultiple_L'] == 0) & (t_prim['isMultiple_S'] == 0)
    assert(len(t_prim) == (sum(binary_filt) + sum(single_filt)))

    if delta_m_cut is not None:
        delta_m_cut = ((t_prim['bin_delta_m'] > 0.1) & binary_filt) | ((t_prim['delta_m_{}'.format(filter_name)] > 0.1) & single_filt)
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
