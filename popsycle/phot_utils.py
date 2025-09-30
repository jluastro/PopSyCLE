def make_filt_dict():
    """
    Dictionary for extinction law coefficients f_i, as a function of filter
    
    Damineli values from photometric bands (nm):
    B = 445, V = 551, I = 806, J = 1220, H = 1630, K = 2190, U = 365, R = 658
    
    SDSS photometric bands https://skyserver.sdss.org/dr1/en/proj/advanced/color/sdssfilters.asp (nm):
    u = 354.3, g = 477.0, r = 623.1, i = 762.5, z = 913.4
    
    ZTF photometric bands:
    G = 472.274, R = 633.961, I = 788.613
    
    RUBIN photometric bands https://github.com/lsst-pst/syseng_throughputs (nm):
    u = 372.0, g = 480.3, r = 622.0, i = 755.8, z = 868.0, y = 974.9
    
    Schlegel and Schlafly photometric bands:
    B = 440, V = 543, I = 809, J = 1266, H = 1673, K = 2215, U = 337, R = 651
    Calculated using calc_f
    """
    filt_dict = {}
    filt_dict['ubv_J'] = {'Schlafly11': 0.709, 'Schlegel99': 0.902, 'Damineli16': 0.662}
    filt_dict['ubv_H'] = {'Schlafly11': 0.449, 'Schlegel99': 0.576, 'Damineli16': 0.344}
    filt_dict['ubv_K'] = {'Schlafly11': 0.302, 'Schlegel99': 0.367, 'Damineli16': 0.172}
    filt_dict['ubv_U'] = {'Schlafly11': 4.334, 'Schlegel99': 5.434, 'Damineli16': 5.022}
    filt_dict['ubv_B'] = {'Schlafly11': 3.626, 'Schlegel99': 4.315, 'Damineli16': 3.757}
    filt_dict['ubv_V'] = {'Schlafly11': 2.742, 'Schlegel99': 3.315, 'Damineli16': 2.757}
    filt_dict['ubv_I'] = {'Schlafly11': 1.505, 'Schlegel99': 1.940, 'Damineli16': 1.496}
    filt_dict['ubv_R'] = {'Schlafly11': 2.169, 'Schlegel99': 2.634, 'Damineli16': 2.102}
    filt_dict['ztf_g'] = {'Damineli16': 3.453}
    filt_dict['ztf_r'] = {'Damineli16': 2.228}
    filt_dict['ztf_i'] = {'Damineli16': 1.553}
    filt_dict['sdss_u'] = {'Damineli16': 5.262}
    filt_dict['sdss_g'] = {'Damineli16': 3.401}
    filt_dict['sdss_r'] = {'Damineli16': 2.290}
    filt_dict['sdss_i'] = {'Damineli16': 1.650}
    filt_dict['sdss_z'] = {'Damineli16': 1.192}
    filt_dict['rubin_u'] = {'Damineli16': 4.880}
    filt_dict['rubin_g'] = {'Damineli16': 3.370}
    filt_dict['rubin_r'] = {'Damineli16': 2.296}
    filt_dict['rubin_i'] = {'Damineli16': 1.672}
    filt_dict['rubin_z'] = {'Damineli16': 1.309}
    filt_dict['rubin_y'] = {'Damineli16': 1.051}
    filt_dict['roman_f062'] = {'Damineli16': 2.307}
    filt_dict['roman_f087'] = {'Damineli16': 1.307}
    filt_dict['roman_f106'] = {'Damineli16': 0.889}
    filt_dict['roman_f129'] = {'Damineli16': 0.583}
    filt_dict['roman_f158'] = {'Damineli16': 0.372}
    filt_dict['roman_f146'] = {'Damineli16': 0.441} #same as w146
    filt_dict['roman_f184'] = {'Damineli16': 0.258}
    filt_dict['roman_f213'] = {'Damineli16': 0.184}

    return filt_dict

def make_photometric_system_dict():
    """
    Dictionary for listing out supported photometric systems and filters
    """
    photometric_system_dict = {}
    photometric_system_dict['ubv'] = ['J', 'H', 'K', 'U', 'B', 'V', 'I', 'R']
    photometric_system_dict['ztf'] = ['g', 'r', 'i']
    photometric_system_dict['sdss'] = ['u', 'g', 'r', 'i', 'z']
    photometric_system_dict['rubin'] = ['u','g','r','i','z','y']
    photometric_system_dict['roman'] = ['f062','f087','f106','f129','f158','f146','f184','f213'] #f146 same as w146

    return photometric_system_dict