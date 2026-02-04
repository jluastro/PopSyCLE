import pdb

import numpy as np
from bagle import model
from astropy.coordinates import SkyCoord
from astropy import units as unit
from astropy.table import Table
from astropy import table
from popsycle import synthetic, binary_utils, phot_utils
from multiprocessing import Pool, Value, Lock

# def get_used_lightcurve_companions(event_table, comp_table, lcurve_table):
# TODO: Finish this... not working at all yet.
#     """Trim the companions and lightcurve tables down to just those that
#     were involved in the 'used_lightcurve'. This is useful for ensuring you will
#     only have one lightcurve per event in the event table. Get rid of all
#     the irrelevant companions that only contribute blend flux.
#
#     Parameters
#     ----------
#     event_table : astropy.table.Table
#     comp_table : astropy.table.Table
#     lcurve_table : astropy.table.Table
#
#     Returns
#     -------
#
#     """
#     # Set event table index for easier cross-matching.
#     event_table_df = event_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'], drop=False)
#
#     # Convert other tables to pandas.
#     lcurv_table_df = lcurve_table.to_pandas()
#     comps_table_df = comp_table.to_pandas()
#
#     # Group lightcurves associated with the same event.
#     grouped_lcurv = lcurv_table_df.groupby(['obj_id_L', 'obj_id_S'])
#
#     # Filter each group to just the companions used in the lightcurve.
#     # Should be 1 lens companion for PSBL, 1 source companion for BSPL, 1 lens + 1 source for BSBL
#     lcurv_used = lcurv_table_df.loc[grouped_lcurv['used_lightcurve'].idxmax()]
#
#     # Clean up columns that should be ints not floats and set index.
#     lcurv_used['obj_id_L'] = lcurv_used['obj_id_L'].astype('int')
#     lcurv_used['obj_id_S'] = lcurv_used['obj_id_S'].astype('int')
#     lcurv_used.set_index(['obj_id_L', 'obj_id_S'], inplace=True)
#
#     comps_table_df['companion_idx'] = comps_table_df['companion_idx'].astype('float')
#
#     # Group the companions by the index for easier access.
#     grouped_comps = comps_table_df.groupby(['obj_id_L', 'obj_id_S'])
#
#     def f_is_used_companion(group):
#         lcurve_grp = lcurv_used.loc[[group['obj_id_L'], group['obj_id_S']]]
#         group['companion_idx'] == lcurve_grp['companion_id_L']
#
#     in_used = grouped_comps.
#
#     grouped_comps_used = grouped_comps.loc[]
#
#     for i in range(len(event_table_df)):
#         event_i = Table.from_pandas(event_table_df.iloc[[i]]) # return astropy table, not series
#         index_i = event_table_df.index[i]
#
#         try:
#             # Get the used lightcurve row.
#             lcurve_i = lcurv_used.loc[index_i]
#
#             # Get the individual companions associated with this used lightcurve.
#             comps_i_all = grouped_comps.get_group(index_i)
#             comps_i = comps_i_all.loc[(comps_i_all['companion_idx'] == lcurve_i['companion_id_L']) |
#                                       (comps_i_all['companion_idx'] == lcurve_i['companion_id_S'])]
#             comps_i = Table.from_pandas(comps_i)
#
#         except KeyError:
#             comps_i = None
#
#     return


def get_bagle_model_list(event_table, comp_table, lcurve_table,
                         photometric_system, filter_name, red_law, n_multi_proc=6):
    """
    Create BAGLE model instances for table of events.

    Parameters
    ----------
    event_table : astropy.table.Table
        Event table from synthetic.py refine_events.

    comp_table : astropy.table.Table
        Companions table that contains all lens or source companions for the events table.

    lcurve_table : astropy.table.Table
        Lightcurves generated for the binary events (from refine_binary_events). This is needed
        to figure out which of the companions (in the case of triples) is used in advance.

    photometric_system : str
        Name of the photometric system, i.e. 'ubv'.

    filter_name : str
        Name of filter associated with photometric system, i.e. 'I'.

    red_law : str
        Name of reddening law in filt_dict list above, i.e. 'Damineli16'.

    Returns
    -------
    A list of BAGLE model instances.

    """
    # Set event table index for easier cross-matching.
    event_table_df = event_table.to_pandas().set_index(['obj_id_L', 'obj_id_S'], drop=False)

    if lcurve_table is None:
        pass
    else:
        # Convert other tables to pandas.
        lcurv_table_df = lcurve_table.to_pandas()

        # Group lightcurves associated with the same event.
        grouped_lcurv = lcurv_table_df.groupby(['obj_id_L', 'obj_id_S'])

        # Filter each group to just the companions used in the lightcurve.
        # Should be 1 lens companion for PSBL, 1 source companion for BSPL, 1 lens + 1 source for BSBL
        lcurv_used = lcurv_table_df.loc[grouped_lcurv['used_lightcurve'].idxmax()]

        # Clean up columns that should be ints not floats and set index.
        lcurv_used['obj_id_L'] = lcurv_used['obj_id_L'].astype('int')
        lcurv_used['obj_id_S'] = lcurv_used['obj_id_S'].astype('int')
        lcurv_used.set_index(['obj_id_L', 'obj_id_S'], inplace=True)

    if comp_table is None:
        pass
    else:
        # Convert other tables to pandas
        comps_table_df = comp_table.to_pandas()
    
        comps_table_df['companion_idx'] = comps_table_df['companion_idx'].astype('float')

        # Group the companions by the index for easier access.
        grouped_comps = comps_table_df.groupby(['obj_id_L', 'obj_id_S'])

    # Split up the events and companions for use in multiprocessing
    inputs = np.empty(len(event_table), dtype=object)

    for i in range(len(event_table_df)):
        event_i = Table.from_pandas(event_table_df.iloc[[i]]) # return astropy table, not series
        index_i = event_table_df.index[i]

        if ((lcurve_table is None) or (comp_table is None)):
            inputs[i] = [event_i, None, photometric_system, filter_name, red_law]
        else:
            try:
                # Get the used lightcurve row.
                lcurve_i = lcurv_used.loc[index_i]

                # Get the individual companions associated with this used lightcurve.
                comps_i_all = grouped_comps.get_group(index_i)
                comps_i = comps_i_all.loc[(comps_i_all['companion_idx'] == lcurve_i['companion_id_L']) |
                                      (comps_i_all['companion_idx'] == lcurve_i['companion_id_S'])]
                comps_i = Table.from_pandas(comps_i)

            except KeyError:
                comps_i = None

            inputs[i] = [event_i, comps_i, photometric_system, filter_name, red_law]

    if n_multi_proc > 1:
        # Set up the multiprocessing
        pool = Pool(n_multi_proc)

        # Generate model instances for all the events.
        results = pool.starmap(get_bagle_model, inputs)
        pool.close()
        pool.join()

        all_models = results
    else:
        all_models = []
        for i in range(len(inputs)):
            all_models.append( get_bagle_model(*inputs[i]) )

    return all_models


def get_bagle_model(event, companions, photometric_system, filter_name, red_law):
    """
    Get a BAGLE model instance for a single event (and its associated companions).

    Parameters
    ----------
    event : astropy.table.Table
        A row from an astropy Table of events. Usually this is the produce of refine_events
        and refine_binary_events(). The length of the table should be 1.

    companions : astropy.table.Table
        A table of the companions associated with the above event. This can be None if
        no binaries were simulated. Note, the companions should only be those needed
        to generate the event. If an event has triples involved, the irrelevant (or less
        important) companions should be trimmed first.

    photometric_system : str

    filter_name : str

    red_law : str

    Returns
    -------

    """
    event = event[0]

    model_name, parameter_dict = get_bagle_model_name_and_params(event, companions,
                                                                 photometric_system, filter_name, red_law)

    mod_class = getattr(model, model_name)
    mod = mod_class(**parameter_dict)

    return mod

def coords_and_prop_motion(event):
    # Get the coordinates of this event.
    L_coords = SkyCoord(l=event['glon_L'] * unit.degree,
                        b=event['glat_L'] * unit.degree,
                        pm_l_cosb=event['mu_lcosb_L'] * unit.mas / unit.year,
                        pm_b=event['mu_b_L'] * unit.mas / unit.year, frame='galactic')
    S_coords = SkyCoord(l=event['glon_S'] * unit.degree,
                        b=event['glat_S'] * unit.degree,
                        pm_l_cosb=event['mu_lcosb_S'] * unit.mas / unit.year,
                        pm_b=event['mu_b_S'] * unit.mas / unit.year, frame='galactic')

    raL = L_coords.icrs.ra.value  # Lens R.A.
    decL = L_coords.icrs.dec.value  # Lens dec
    muL = np.array([L_coords.icrs.pm_ra_cosdec.value, L_coords.icrs.pm_dec.value])  # lens proper motion mas/year
    muS = np.array([S_coords.icrs.pm_ra_cosdec.value, S_coords.icrs.pm_dec.value])  # source proper motion mas/year

    return raL, decL, muL, muS

def get_pspl_lightcurve_parameters(events, filter_dict, event_id = None):
    """
    Find the parameters for PSPL_PhotAstrom_Par_Param1 from 
    event_table.

    Parameters
    ----------
    events : Astropy table
        Table containing the events calculated from refine_events.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    event_id : float or None, optional
        Index of event table of event. If len(events) > 1, this must be specified.
        Default is None.
        
    Returns
    -------
    parameter_dict : dict
        Dictionary of the PSPL_PhotAstrom_Par_Param1 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event

    model_name : str
        Name of model associated with event.
    """
    
    if len(events) == 1 or type(events) == table.row.Row:
        event = events
    else:
        if event_id == None:
            raise Exception('If you input more than one event, must specify event id')
        else:
            event = events[event_id]
            
    obj_id_L = event['obj_id_L']
    obj_id_S = event['obj_id_S']

    raL, decL, muL, muS = coords_and_prop_motion(event)
            
    model_name = 'PSPL_PhotAstrom_Par_Param1'

    mL = event['mass_L']  # Msun (Primary lens current mass)
    t0 = event['t0']  # mjd
    xS0 = np.array([0, 0])  # arbitrary offset (arcsec)
    beta = event['u0'] * event['theta_E']  # mas
    dL = event['rad_L'] * 10 ** 3  # Distance to lens
    dS = event['rad_S'] * 10 ** 3  # Distance to source
    mag_src = []
    b_sff = []
    for photometric_system in filter_dict:
        for filter_name in filter_dict[photometric_system]:
            mag_src += [event['%s_%s_app_S' % (photometric_system, filter_name)]]
            b_sff += [event['f_blend_%s' % filter_name]]

    parameter_dict = {'raL': raL, 'decL': decL, 'mL': mL,
                      't0': t0, 'beta': beta, 'dL': dL, 'dL_dS': dL / dS,
                      'xS0_E': xS0[0], 'xS0_N': xS0[1],
                      'muL_E': muL[0], 'muL_N': muL[1], 'muS_E': muS[0], 'muS_N': muS[1],
                      'b_sff': b_sff, 'mag_src': mag_src}

    return parameter_dict, obj_id_L, obj_id_S, model_name

def get_psbl_lightcurve_parameters(events, companions, comp_idx_L, filter_dict, event_id = None):
    """
    Find the parameters for PSBL_PhotAstrom_Par_EllOrbs_Param7 from 
    event_table and comp_table.

    Parameters
    ----------
    events : Astropy table
        Table containing the events calculated from refine_events.
    
    companions : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx_L : int
        Index into the comp_table of the companion for which the psbl is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id.
        Default is None.
        
    Returns
    -------
    parameter_dict : dict
        Dictionary of the PSBL_PhotAstrom_EllOrbs_Par_Param7 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event

    model_name : str
        Name of model associated with event.
        
    """
    if (type(comp_idx_L) != int) and (type(comp_idx_L) != np.int64):
        raise Exception('comp_idx_L must be an integer')
        
    if type(events) == table.row.Row:
        event = events
        obj_id_L = event['obj_id_L']
        obj_id_S = event['obj_id_S']
    else:
        if event_id == None:
            obj_id_L = companions['obj_id_L'][comp_idx_L]
            obj_id_S = companions['obj_id_S'][comp_idx_S]
            event_id = (np.where(np.logical_and((events['obj_id_L'] == obj_id_L), (events['obj_id_S'] == obj_id_S)))[0])[0]
            event = events[event_id]
        else:
            event = events[event_id]
            obj_id_L = event['obj_id_L']
            obj_id_S = event['obj_id_S']

    model_name = 'PSBL_PhotAstrom_Par_EllOrbs_Param7'
    
    raL, decL, muL, muS = coords_and_prop_motion(event)

    mLp = event['mass_L']  # msun (Primary lens current mass)
    mLs = companions['mass'][comp_idx_L]  # msun (Companion lens current mass)
    t0_p = event['t0']  # mjd
    xS0 = np.array([0, 0])  # arbitrary offset (arcsec)
    beta_p = event['u0'] * event['theta_E']  # 5.0
    dL = event['rad_L'] * 10 ** 3  # Distance to lens
    dS = event['rad_S'] * 10 ** 3  # Distance to source
#     mag_src = [event['%s_%s_app_S' % (photometric_system, filter_name)]]
#     b_sff = [event['f_blend_%s' % filter_name]]  # ASSUMES ALL BINARY LENSES ARE BLENDED
    omega = companions['omega'][comp_idx_L]
    big_omega = companions['Omega'][comp_idx_L]
    i = companions['i'][comp_idx_L]
    e = companions['e'][comp_idx_L]
    tp = companions['tp'][comp_idx_L]
    a = 10**(companions['log_a'][comp_idx_L])
#     dmag_Lp_Ls = [event['%s_%s_L' % (photometric_system, filter_name)] - companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_L]]
    
    #FIX HERE
    mag_src = []
    b_sff = []
    dmag_Lp_Ls = []
    for photometric_system in filter_dict:
        for filter_name in filter_dict[photometric_system]:
            filt = photometric_system+'_'+filter_name
            mag_src += [event['%s_%s_app_S' % (photometric_system, filter_name)]]
            b_sff += [event['f_blend_%s' % filter_name]]
            dmag_Lp_Ls += [event['%s_%s_L' % (photometric_system, filter_name)] - companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_L]]
            

    parameter_dict = {'raL': raL, 'decL': decL,
                      'mLp': mLp, 'mLs': mLs, 't0_p': t0_p,
                      'xS0_E': xS0[0], 'xS0_N': xS0[1], 'beta_p': beta_p,
                      'muL_E': muL[0], 'muL_N': muL[1], 'muS_E': muS[0], 'muS_N': muS[1],
                      'dL': dL, 'dS': dS, 
                      #'sep': sep, 'alpha': alpha, 
                      'mag_src': mag_src, 'b_sff': b_sff, 'omega_pri': omega, 'big_omega_sec': big_omega,
                     'i': i, 'e': e, 'tp': tp, 'a': a, 'dmag_Lp_Ls': dmag_Lp_Ls}

    return parameter_dict, obj_id_L, obj_id_S, model_name

def get_bspl_lightcurve_parameters(events, companions, comp_idx_S, filter_dict, red_law, event_id = None):
    """
    Find the parameters for BSPL_PhotAstrom_Par_EllOrbs_Param4 from 
    event_table and comp_table.

    Parameters
    ----------
    events : Astropy table
        Table containing the events calculated from refine_events.
    
    companions : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx_S : int
        Index into the comp_table of the companion for which the bspl is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    red_law : str
        Redenning law
    
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id
        
    Returns
    -------
    parameter_dict : dict
        Dictionary of the BSPL_PhotAstrom_Par_EllOrbs_Param4 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event

    model_name : str
        Name of model associated with event.
    """
    if (type(comp_idx_S) != int) and (type(comp_idx_S) != np.int64):
        raise Exception('comp_idx_S must be an integer')
        
    if type(events) == table.row.Row:
        event = events
        obj_id_L = event['obj_id_L']
        obj_id_S = event['obj_id_S']
    else:
        if event_id == None:
            obj_id_L = companions['obj_id_L'][comp_idx_L]
            obj_id_S = companions['obj_id_S'][comp_idx_S]
            event_id = (np.where(np.logical_and((events['obj_id_L'] == obj_id_L), (events['obj_id_S'] == obj_id_S)))[0])[0]
            event = events[event_id]
        else:
            event = events[event_id]
            obj_id_L = event['obj_id_L']
            obj_id_S = event['obj_id_S']
    
    raL, decL, muL, muS = coords_and_prop_motion(event)
    
    model_name = 'BSPL_PhotAstrom_Par_EllOrbs_Param4'

#     filt_dict = phot_utils.make_filt_dict()
#     f_i = filt_dict[photometric_system + '_' + filter_name][red_law]
#     abs_mag_sec = companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_S]

    mL = event['mass_L']  # msun (Lens current mass)
    t0_p = event['t0']  # mjd
    beta_p = event['u0'] * event['theta_E']  # 5.0
    dL = event['rad_L'] * 10 ** 3  # Distance to lens
    dL_dS = dL / (event['rad_S'] * 10 ** 3)  # Distance to lens/Distance to source
    xS0 = np.array([0, 0])  # arbitrary offset (arcsec)
#     mag_src_sec = synthetic.calc_app_mag(event['rad_S'], abs_mag_sec, event['exbv_S'], f_i)
#     mag_src_pri = binary_utils.subtract_magnitudes(
#         event['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec)
#     b_sff = event['f_blend_%s' % filter_name]  # ASSUMES THAT SOURCE BINARIES ARE BLENDED
    omega = companions['omega'][comp_idx_S]
    big_omega = companions['Omega'][comp_idx_S]
    i = companions['i'][comp_idx_S]
    e = companions['e'][comp_idx_S]
    log_a = companions['log_a'][comp_idx_S]
    tp = companions['tp'][comp_idx_S]
    mass_source_p = event['mass_S']
    mass_source_s = companions['mass'][comp_idx_S]
    
    #FIX HERE
    b_sff = []
    mag_src_sec = []
    mag_src_pri = []
    filt_dict = phot_utils.make_filt_dict()
    for photometric_system in filter_dict:
        for filter_name in filter_dict[photometric_system]:
            filt = photometric_system+'_'+filter_name
            f_i = filt_dict[photometric_system + '_' + filter_name][red_law]
            abs_mag_sec = companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_S]
            mag_src_sec_tmp = synthetic.calc_app_mag(event['rad_S'], abs_mag_sec, event['exbv_S'], f_i)
            
            mag_src_sec += [mag_src_sec_tmp]
            mag_src_pri += [binary_utils.subtract_magnitudes(
                event['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec_tmp)] 
            b_sff += [event['f_blend_%s' % filter_name]]
    
    parameter_dict = {'raL': raL, 'decL': decL, 'mL': mL,
                      't0': t0_p, 'beta': beta_p, 'dL': dL, 'dL_dS': dL_dS,
                      'xS0_E': xS0[0], 'xS0_N': xS0[1],
                      'muL_E': muL[0], 'muL_N': muL[1], 'muS_E': muS[0], 'muS_N': muS[1],
                      'mag_src_pri': mag_src_pri, 'mag_src_sec': mag_src_sec, 'b_sff': b_sff, 'omega_pri': omega, 'big_omega_sec': big_omega, 
                      'i': i, 'e': e, 'log_a': log_a, 'mass_source_p': mass_source_p, 'mass_source_s': mass_source_s, 'tp': tp}

    return parameter_dict, obj_id_L, obj_id_S, model_name

def get_bsbl_lightcurve_parameters(events, companions, comp_idx_L, comp_idx_S, filter_dict, red_law, event_id = None):
    """
    Find the parameters for BSBL_PhotAstrom_Par_EllOrbs_Param3 from 
    event_table and comp_table.

    Parameters
    ----------
    events : Astropy table
        Table containing the events calculated from refine_events.
    
    companions : Astropy table
        Table containing the companions calculated from refine_events.
    
    comp_idx_L : int
        Index into the comp_table of the lens companion for which the model is being calculated.
        
    comp_idx_S : int
        Index into the comp_table of the source companion for which the model is being calculated.
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.
    
    red_law : str
        Redenning law
        
    event_id : float or None, optional
        Corresponding event_id in event_table to companion id.
        Default is None.
        
    Returns
    -------
    parameter_dict : dict
        Dictionary of the BSBL_PhotAstrom_Par_EllOrbs_Param3 parameters
        
    obj_id_L : int
        Object id of the lens associated with event
        
    obj_id_S : int
        Object id of the source associated with event

    model_name : str
        Name of model associated with event.
    """
    if ((type(comp_idx_L) != int) and (type(comp_idx_L) != np.int64)) or ((type(comp_idx_S) != int) and (type(comp_idx_S) != np.int64)):
        raise Exception('comp_idx_L and comp_idx_S must be integers')
    
    if type(events) == table.row.Row:
        event = events
        obj_id_L = event['obj_id_L']
        obj_id_S = event['obj_id_S']
    else:
        if event_id == None:
            obj_id_L = companions['obj_id_L'][comp_idx_L]
            obj_id_S = companions['obj_id_S'][comp_idx_S]
            event_id = (np.where(np.logical_and((events['obj_id_L'] == obj_id_L), (events['obj_id_S'] == obj_id_S)))[0])[0]
            event = events[event_id]
        else:
            event = events[event_id]
            obj_id_L = event['obj_id_L']
            obj_id_S = event['obj_id_S']
    
    raL, decL, muL, muS = coords_and_prop_motion(event)

    model_name = 'BSBL_PhotAstrom_Par_EllOrbs_Param3'

#     f_i = synthetic.filt_dict[photometric_system + '_' + filter_name][red_law]
#     abs_mag_sec = companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_S]
    
    mLp = event['mass_L']  # msun (Lens current mass)
    mLs = companions['mass'][comp_idx_L]  # msun (Companion lens current mass)
    t0_p = event['t0']  # mjd
    beta_p = event['u0'] * event['theta_E']  # 5.0
    dL = event['rad_L'] * 10 ** 3  # Distance to lens
    dS = event['rad_S'] * 10 ** 3  # Distance to source
    xS0_E = 0.0  # arbitrary offset (arcsec)
    xS0_N = 0.0  # arbitrary offset (arcsec)
#     mag_src_sec = synthetic.calc_app_mag(event['rad_S'], abs_mag_sec, event['exbv_S'], f_i)
#     mag_src_pri = binary_utils.subtract_magnitudes(
#         event['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec)
#     b_sff = event['f_blend_%s' % filter_name]  # ASSUMES THAT SOURCE BINARIES ARE BLENDED
    omegaL = companions['omega'][comp_idx_L]
    big_omegaL = companions['Omega'][comp_idx_L]
    iL = companions['i'][comp_idx_L]
    eL = companions['e'][comp_idx_L]
    tpL = companions['tp'][comp_idx_L]
    aL = 10**(companions['log_a'][comp_idx_L])
    omegaS = companions['omega'][comp_idx_S]
    big_omegaS = companions['Omega'][comp_idx_S]
    iS = companions['i'][comp_idx_S]
    eS = companions['e'][comp_idx_S]
    tpS = companions['tp'][comp_idx_S]
    aS = 10**(companions['log_a'][comp_idx_S])
#     dmag_Lp_Ls = [event['%s_%s_L' % (photometric_system, filter_name)] - companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_L]]
    mass_source_p = event['mass_S']
    mass_source_s = companions['mass'][comp_idx_S]
    
    #FIX HERE
    b_sff = []
    mag_src_sec = []
    mag_src_pri = []
    dmag_Lp_Ls = []
    filt_dict = phot_utils.make_filt_dict()
    for photometric_system in filter_dict:
        for filter_name in filter_dict[photometric_system]:
            filt = photometric_system+'_'+filter_name
            f_i = filt_dict[photometric_system + '_' + filter_name][red_law]
            abs_mag_sec = companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_S]
            
            mag_src_sec_tmp = synthetic.calc_app_mag(event['rad_S'], abs_mag_sec, event['exbv_S'], f_i)
            mag_src_sec += [mag_src_sec_tmp]
            mag_src_pri += [binary_utils.subtract_magnitudes(
                event['%s_%s_app_S' % (photometric_system, filter_name)], mag_src_sec_tmp)]
            b_sff += [event['f_blend_%s' % filter_name]]
            dmag_Lp_Ls += [event['%s_%s_L' % (photometric_system, filter_name)] - companions['m_%s_%s' % (photometric_system, filter_name)][comp_idx_L]]
            

    parameter_dict = {'raL': raL, 'decL': decL, 'mLp': mLp, 'mLs': mLs,
                      't0_p': t0_p, 'xS0_E': xS0_E, 'xS0_N': xS0_N, 'beta_p': beta_p,
                      'muL_E': muL[0], 'muL_N': muL[1], 'muS_E': muS[0], 'muS_N': muS[1],
                      'dL': dL, 'dS': dS, 
                      'mag_src_pri': mag_src_pri, 'mag_src_sec': mag_src_sec, 'b_sff': b_sff,
                     'omegaL_pri': omegaL, 'big_omegaL_sec': big_omegaL, 'iL': iL, 'eL': eL, 'tpL': tpL, 'aL': aL,
                     'omegaS_pri': omegaS, 'big_omegaS_sec': big_omegaS, 'iS': iS, 'eS': eS, 'tpS': tpS, 'aS': aS,
                     'dmag_Lp_Ls': dmag_Lp_Ls, 'mass_source_p': mass_source_p, 'mass_source_s': mass_source_s}
    

    return parameter_dict, obj_id_L, obj_id_S, model_name

def get_bagle_model_name_and_params(event, companions, photometric_system, filter_name, red_law):
    """
    For a single event and its associated companions, get the BAGLE model name
    and parameters (in a dictionary).

    Parameters
    ----------
    event : astropy.table.Table
        A row from an astropy Table of events. Usually this is the produce of refine_events
        and refine_binary_events(). The length of the table should be 1.

    companions : astropy.table.Table
        A table of the companions associated with the above event. This can be None if
        no binaries were simulated. Note, the companions should only be those needed
        to generate the event. If an event has triples involved, the irrelevant (or less
        important) companions should be trimmed first.

    photometric_system : str

    filter_name : str

    red_law : str

    Returns
    -------
    model_name : str
        Model name.

    parameter_dict : dict
        Dictionary of parameters (both fit and fixed).
    """
    # Determine the class of model to use.
    multi_L = event['isMultiple_L']
    multi_S = event['isMultiple_S']

    if multi_L == 1 and multi_S == 1:
        event_type = 'BSBL'
    elif multi_L == 1 and multi_S == 0:
        event_type = 'PSBL'
    elif multi_L == 0 and multi_S == 1:
        event_type = 'BSPL'
    else:
        event_type = 'PSPL'

    if event_type == 'PSPL':
        parameter_dict, obj_id_L, obj_id_S, model_name = get_pspl_lightcurve_parameters(event, photometric_system, filter_name)

    if event_type == 'PSBL':
        # There should only be a single lens companion.
        comp_idxs_L = np.where(companions['prim_type'] == b"L")[0]

        if len(comp_idxs_L) > 1:
            raise RuntimeError('Found too many lens companions. Triples not supported.')
        else:
            comp_idx_L = comp_idxs_L[0]

        parameter_dict, obj_id_L, obj_id_S, model_name = get_psbl_lightcurve_parameters(event, companions, comp_idx_L, photometric_system, filter_name)

    if event_type == 'BSPL':
        comp_idxs_S = np.where(companions['prim_type'] == b"S")[0]
    
        if len(comp_idxs_S) > 1:
            raise RuntimeError('Found too many source companions. Triples not supported.')
        else:
            comp_idx_S = comp_idxs_S[0]

        parameter_dict, obj_id_L, obj_id_S, model_name = get_bspl_lightcurve_parameters(event, companions, comp_idx_S, photometric_system, filter_name, red_law)

    if event_type == 'BSBL':
        # There should only be a single source companion.
        comp_idxs_S = np.where(companions['prim_type'] == b"S")[0]
        comp_idxs_L = np.where(companions['prim_type'] == b"L")[0]

        if len(comp_idxs_S) > 1:
            raise RuntimeError('Found too many source companions. Triples not supported.')
        else:
            comp_idx_S = comp_idxs_S[0]

        if len(comp_idxs_L) > 1:
            raise RuntimeError('Found too many lens companions. Triples not supported.')
        else:
            comp_idx_L = comp_idxs_L[0]

        parameter_dict, obj_id_L, obj_id_S, model_name = get_bsbl_lightcurve_parameters(event, companions, comp_idx_L, comp_idx_S, photometric_system, filter_name, red_law)

    return model_name, parameter_dict
