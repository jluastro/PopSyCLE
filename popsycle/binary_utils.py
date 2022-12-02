import numpy as np
from astropy.table import Table


from microlens.jlu import model


def refine_bspl_events(events, companions, photometric_system, filter_name,
                         overwrite = False, output_file = 'default',
                         save_phot = False, phot_dir = None):
    """
    Takes the output Astropy table from refine_events (both primaries and companions) and from that
    calculates the binary light curves.

    Parameters
    ----------
    events : str
        fits file containing the events calculated from refine_events
    
    companions : str
        fits file containing the companions calculated from refine_events
    
    photometric_system : str
        The name of the photometric system in which the filter exists.
    
    filter_name : str
        The name of the filter in which to calculate all the
        microlensing events. The filter name convention is set
        in the global filt_dict parameter at the top of this module.

    Optional Parameters
    -------------------
    overwrite : bool
        If set to True, overwrites output files. If set to False, exists the
        function if output files are already on disk.
        Default is False.

    output_file : str
        The name of the final refined_events file.
        If set to 'default', the format will be:
            <input_root>_refined_events_<photometric_system>_<filt>_<red_law>.fits
            
    save_phot : bool
        If set to True, saves the photometry generated instead of just parameters.
        Default is False
    
    phot_dir : str
        Name of the directory photometry is saved if save_phot = True.
        This parameters is NOT optional if save_phot = True.
        Default is None.
    
    Output:
    ----------
    A file will be created named
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>_companions_rb.fits
    that contains all the same objects, only now with lots of extra
    columns of data. (rb stands for refine binaries).
    
    A file will be created named
    <input_root>_refined_events_<photometric_system>_<filt>_<red_law>_companions_rb_mp.fits
    that contains the data for each individual peak for events with multiple peaks.
    (mp stands for multiple peaks).

    """
    start_time = time.time()
    
    if not overwrite and os.path.isfile(output_file):
        raise Exception('That refined_events.fits file name is taken! '
                        'Either delete the .fits file, or pick a new name.')

    if save_phot == True and phot_dir == None:
        raise Exception('phot_dir is "none". Input a directory to save photometry.')
        
    # Error handling/complaining if input types are not right.
    _check_refine_binary_events(events, companions, 
                         photometric_system, filter_name,
                         overwrite, output_file,
                         save_phot, phot_dir)
        
    
    event_table = Table.read(events)
    comp_table = Table.read(companions)
    
    comp_table.add_column( Column(np.zeros(len(comp_table), dtype=float), name='n_peaks') )
    comp_table.add_column( Column(np.zeros(len(comp_table), dtype=float), name='bin_delta_m') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='tE_sys') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='tE_primary') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='primary_t') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='avg_t') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='std_t') )
    comp_table.add_column( Column(np.empty(len(comp_table), dtype=float), name='asymmetry') )
    
    comp_table['bin_delta_m'][:] = np.nan
    comp_table['tE_sys'][:] = np.nan
    comp_table['tE_primary'][:] = np.nan
    comp_table['primary_t'][:] = np.nan
    comp_table['avg_t'][:] = np.nan
    comp_table['std_t'][:] = np.nan
    comp_table['asymmetry'][:] = np.nan
    
    # This table is for events with more than one peak to characterize those peaks
    # comp_id is position of companion in companion table for reference
    # obj_id_L, obj_id_S, and n_peaks are the same as in the companion table, just for reference
    # t is time of peak
    # tE is Einstein crossing time of peak defined by times of 0.5*(max(peak mag) - min(peak mag))
    # delta m is the change in magnitude between the peak and baseline
    # ratio is the magnitude ratio between min peak/max peak
    mult_peaks = Table(names=('comp_id', 'obj_id_L', 'obj_id_S', 'n_peaks', 't', 'tE', 'delta_m', 'ratio'))

    L_idxs = np.where(comp_table['prim_type'] == "L")[0]
    
    for comp_idx in L_idxs:
        name = "{}".format(comp_idx)
        event_id = (np.where(np.logical_and((event_table['obj_id_L'] == comp_table[comp_idx]['obj_id_L']), (event_table['obj_id_S'] == comp_table[comp_idx]['obj_id_S'])))[0])[0]
        L_coords = SkyCoord(l = event_table[event_id]['glat_L']*unit.degree, b = event_table[event_id]['glon_L']*unit.degree, 
                                pm_l_cosb = event_table[event_id]['mu_lcosb_L']*unit.mas/unit.year, 
                                pm_b = event_table[event_id]['mu_b_L']*unit.mas/unit.year, frame ='galactic')
        S_coords = SkyCoord(l = event_table[event_id]['glat_S']*unit.degree, b = event_table[event_id]['glon_S']*unit.degree, 
                                pm_l_cosb = event_table[event_id]['mu_lcosb_S']*unit.mas/unit.year, 
                                  pm_b = event_table[event_id]['mu_b_S']*unit.mas/unit.year, frame ='galactic')


        ##########
        # Calculate binary model and photometry
        ##########
        raL = L_coords.icrs.ra.value # Lens R.A.
        decL = L_coords.icrs.dec.value # Lens dec
        mL1 = event_table[event_id]['mass_L'] # msun (Primary lens mass)
        mL2 = comp_table[comp_idx]['mass'] # msun (Companion lens mass)
        t0 = event_table[event_id]['t0'] # mjd
        xS0 = np.array([0, 0]) #arbitrary offset (arcsec)
        beta = event_table[event_id]['u0']*event_table[event_id]['theta_E']#5.0
        muL = np.array([L_coords.icrs.pm_ra_cosdec.value, L_coords.icrs.pm_dec.value]) #lens proper motion mas/year
        muS = np.array([S_coords.icrs.pm_ra_cosdec.value, S_coords.icrs.pm_dec.value]) #source proper motion mas/year
        dL = event_table[event_id]['rad_L']*10**3 #Distance to lens
        dS = event_table[event_id]['rad_S']*10**3 #Distance to source
        sep = comp_table[comp_idx]['sep'] #mas (separation between primary and companion)
        alpha = comp_table[comp_idx]['alpha']
        mag_src = event_table[event_id]['%s_%s_app_S' % (photometric_system, filter_name)]
        b_sff = event_table[event_id]['f_blend_%s' % filter_name]

        psbl = model.PSBL_PhotAstrom_Par_Param1(mL1, mL2, t0, xS0[0], xS0[1],
                                   beta, muL[0], muL[1], muS[0], muS[1], dL, dS,
                                   sep, alpha, [b_sff], [mag_src], 
                                   raL=raL, decL=decL, 
                                   root_tol = 0.00000001)
        

        # Calculate the photometry 
        duration=1000 # days
        time_steps=5000
        tmin = psbl.t0 - (duration / 2.0)
        tmax = psbl.t0 + (duration / 2.0)
        dt = np.linspace(tmin, tmax, time_steps)
        
        img, amp = psbl.get_all_arrays(dt)
        phot = psbl.get_photometry(dt, amp_arr=amp)
        
        if save_phot == True:
            if not os.path.exists(phot_dir):
                os.makedirs(phot_dir)
            foo = Table((dt, phot), names=['time', 'phot'])
            foo.write(phot_dir + '/' + name + '_phot.fits', overwrite=overwrite)
        
        #because this is magnitudes max(phot) is baseline and min(phot) is peak
        #baseline 2000tE away get_photometry
        comp_table[comp_idx]['bin_delta_m'] = max(phot) - min(phot)
        tenp = np.where(phot < (max(phot) - 0.1*comp_table[comp_idx]['bin_delta_m']))[0]
        if len(tenp) == 0:
            continue
        comp_table[comp_idx]['tE_sys'] = max(dt[tenp]) - min(dt[tenp])
        
        # Find peaks
        peaks, _ = find_peaks(-phot, prominence = 10e-5, width =1) 
                
        if len(peaks) == 0:            
            continue
        
        comp_table[comp_idx]['n_peaks'] = len(peaks)
        comp_table[comp_idx]['primary_t'] = dt[peaks][np.argmin(phot[peaks])]
        comp_table[comp_idx]['avg_t'] = np.average(dt[peaks])
        comp_table[comp_idx]['std_t'] = np.std(dt[peaks])
        
        # Find asymmetry (0 if symmetric, larger if asymmetric)
        # Uses 50 degree chebyshev polynomial (see eq 7 in Night et al. 2010)
        indices = np.arange(0,51)
        cheb_fit = np.polynomial.chebyshev.Chebyshev.fit(dt, phot, 50).coef
        odd_cheb = cheb_fit[indices%2==1]
        even_cheb = cheb_fit[indices%2==0]
        asymm = np.sqrt(np.sum(odd_cheb**2)/np.sum(even_cheb**2))
        comp_table[comp_idx]['asymmetry'] = asymm
        
        # Split up peaks by minima between peaks
        # Note since it's magnitudes all the np.min and such are maxima
        split_data = []
        start_idx = 0
        if len(peaks) > 1:
            for i in range(len(peaks) - 1):
                min_btwn_peaks = np.max(phot[peaks[i]:peaks[i+1]])
                end_idx = np.where(phot[start_idx:] == min_btwn_peaks)[0][0] + start_idx
                split_data.append([dt[start_idx:end_idx], phot[start_idx:end_idx]])
                start_idx = end_idx
        split_data.append([dt[start_idx:], phot[start_idx:]])
        split_data = np.array(split_data, dtype='object')
        
        highest_peak = np.argmin(phot[peaks])
        
        if len(split_data[highest_peak][1]) == 0:
            continue
        
        highest_bump_mag = max(split_data[highest_peak][1]) - min(split_data[highest_peak][1])
        highest_half = np.where(split_data[highest_peak][1] < (max(split_data[highest_peak][1]) - 0.5*highest_bump_mag))[0]
        if len(highest_half) == 0:
            continue
        comp_table[comp_idx]['tE_primary'] = max(dt[highest_half]) - min(dt[highest_half])
        
        # For events with more than one peak, add them to the multi peak table
        if len(peaks) > 1:
            n_peaks = len(peaks)
            obj_id_L = comp_table[comp_idx]['obj_id_L']
            obj_id_S = comp_table[comp_idx]['obj_id_S']
            for i in range(len(peaks)):
                t = dt[peaks[i]]
                delta_m = max(phot) - phot[peaks[i]]
                ratio = np.min(phot[peaks])/phot[peaks[i]]
                
                # Don't log primary peak
                if ratio == 1:
                    continue
                
                if len(split_data[i][1]) == 0:
                    tE = np.nan
                    mult_peaks.add_row([comp_idx, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
                    continue
                    
                split_bump_mag = max(split_data[i][1]) - min(split_data[i][1])
                split_half = np.where(split_data[i][1] < (max(split_data[i][1]) - 0.5*split_bump_mag))[0]
                if len(split_half) == 0:
                    tE = np.nan
                    mult_peaks.add_row([comp_idx, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
                    continue
                    
                tE = max(dt[split_half]) - min(dt[split_half])
                
                mult_peaks.add_row([comp_idx, obj_id_L, obj_id_S, n_peaks, t, tE, delta_m, ratio])
            
        
    # Writes fits file
    if output_file == 'default':
        comp_table.write(companions[:-5] + "_rb.fits", overwrite=overwrite)
        mult_peaks.write(companions[:-5] + "_rb_multi_peaks.fits", overwrite=overwrite)
    else:
        comp_table.write(output_file + ".fits", overwrite=overwrite)
        mult_peaks.write(output_file + "_multi_peaks.fits", overwrite=overwrite)
        
        
        
    ##########
    # Make log file
    ##########
    now = datetime.datetime.now()
    popsycle_path = os.path.dirname(inspect.getfile(perform_pop_syn))
    popstar_path = os.path.dirname(inspect.getfile(imf))
    popsycle_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                             cwd=popsycle_path).decode('ascii').strip()
    popstar_hash = subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                           cwd=popstar_path).decode('ascii').strip()
    
    end_time = time.time()
    
    dash_line = '-----------------------------' + '\n'
    empty_line = '\n'

    line0 = 'FUNCTION INPUT PARAMETERS' + '\n'
    line1 = 'event_file : ' + events + '\n'
    line2 = 'companion_file : ' + companions + '\n'
    line3 = 'save_phot : ' + str(save_phot) + '\n'

    line4 = 'VERSION INFORMATION' + '\n'
    line5 = str(now) + ' : creation date' + '\n'
    line6 = popstar_hash + ' : SPISEA commit' + '\n'
    line7 = popsycle_hash + ' : PopSyCLE commit' + '\n'

    line8 = 'OTHER INFORMATION' + '\n'
    line9 = str(end_time - start_time) + ' : total runtime (s)' + '\n'
    line10 = str(len(L_idxs)) + ' : number of simulated lightcurves' + '\n'

    line11 = 'FILES CREATED' + '\n'
    if output_file == 'default':
        line12 = companions[:-5] + "_rb.fits" + ' : binary refined events' + '\n'
        line13 = companions[:-5] + "_rb_multi_peaks.fits" + ' : multiple peak events table' + '\n'
        log_name = companions[:-5] + "_rb.log"
    else:
        line12 = output_file + ".fits" + ' : refined events' + '\n'
        line13 = output_file + "_multi_peaks.fits" + ' : multiple peak events table' + '\n'
        log_name = output_file + '.log'
     
    line14 = '\n'
    if save_phot == True:
        line14 = phot_dir + ' : directiory of photometry'
       
    
    with open(log_name, 'w') as out:
        out.writelines([line0, dash_line, line1, line2, line3, empty_line,
                        line4, dash_line, line5, line6, line7, empty_line,
                        line8, dash_line, line9, line10, empty_line,
                        line11, dash_line, line12, line13, line14])

    print('refine_binary_events runtime : {0:f} s'.format(end_time - start_time))
    return
