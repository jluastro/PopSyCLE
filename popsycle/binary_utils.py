import numpy as np
import warnings
from astropy.table import Table
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


#def primary_mag_from_system_mag(system_mag, companion_mag):
#    return subtract_magnitudes(system_mag, companion_mag)

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

    
