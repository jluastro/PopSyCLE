import numpy as np
from popsycle.synthetic import add_magnitudes


def primary_mag_from_system_mag(system_mag, companion_mag):
    return subtract_magnitudes(system_mag, companion_mag)

def primary_mag_from_system_mag_hdf5(prim_hdf5, comp_hdf5, mag_colnames):
    grouped_companions = companions_table.group_by(['system_idx'])
    for mag_col in mag_colnames:
        companions_mag = grouped_companions[mag_col].groups.aggregate(add_magnitudes)
        prim_table[mag_col + '_prim'] = subtract_magnitudes(prim_table[mag_col], companions_mag)
    return prim_table

def primary_mag_from_system_mag_Table(prim_table, companions_table, mag_colnames):
    
    prim_df = prim_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'])
    companion_df = companion_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'])
    joined = companion_df.join(prim_df, lsuffix='_comp', rsuffix='_prim', how='outer')
    
    grouped_companions = companions_table.group_by(['system_idx'])
    for mag_col in mag_colnames:
        companions_mag = grouped_companions[mag_col].groups.aggregate(add_magnitudes)
        prim_table[mag_col + '_prim'] = subtract_magnitudes(prim_table[mag_col], companions_mag)
    return prim_table
    
    
    
    #for mag_col in mag_colnames:
    #    joined[mag_col + '_primary'] = subtract_magnitudes(prim_df[mag_col + '_prim']
    

def subtract_magnitudes(m1, m2):
    
    m1 = np.array(m1)
    m2 = np.array(m2)
    
    m_diff = -2.5*np.log10(10**(-0.4*m1) - 10**(-0.4*m2))
    
    return m_diff