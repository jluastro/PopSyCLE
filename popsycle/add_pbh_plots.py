import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

def print_plots(output_root, galactic_lin_distance, galactic_lin_b, galactic_lin_l, galactocen_lin_spherical_distance,
                galactocen_lin_spherical_b, galactocen_lin_spherical_l, rho_lin, r_max, n_lin, cdf_los,
                x_cyl, y_cyl, r_cyl, r_proj_los_cyl, n_pbh, d_galac, b_galac, l_galac, area_proj_los_cyl,
                mask_obs_cone, field_of_view_radius, l_radian, b_radian, f_cdf_d, pbh_mass):

    ##################################################################################
    #Line of Sight - Galactic Coordinates
    ##################################################################################
    fig_arr = []
    fig1, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0}, figsize=(8,8))
    fig_arr.append(fig1)
    axs[0].plot(galactic_lin_distance, '.', label='Distance', c='C0')
    axs[1].plot(galactic_lin_b, '.', label='Latitude',c='C1')
    axs[2].plot(galactic_lin_l, '.', label='Longitude', c='C2')


    plt.legend(loc='best')
    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
        ax.legend(loc='upper left')

    axs[0].set_title('Line of Sight Linear Space Galactic Coordinates')
    axs[0].set_ylabel('[kpc]')

    axs[1].set_ylabel('[deg]')

    axs[2].set_ylabel('[deg]')
    axs[2].set_xlabel('Line-of-Sight Linear Space Coordinate Index')

    axs[0].set_ylim(0, axs[0].set_ylim()[1])

    ##################################################################################
    #Line of Sight - Galactocentric Coordinates
    ##################################################################################
    fig2, axs = plt.subplots(3, sharex=True, gridspec_kw={'hspace': 0}, figsize=(8,8))
    fig_arr.append(fig2)
    axs[0].plot(galactic_lin_distance,
                galactocen_lin_spherical_distance, '.', label='Distance', c='C0')
    axs[1].plot(galactic_lin_distance,
                galactocen_lin_spherical_b, '.', label='Latitude', c='C1')
    axs[2].plot(galactic_lin_distance,
                galactocen_lin_spherical_l, '.', label='Longitude',c='C2')

    plt.legend(loc='best')
    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
        ax.legend(loc='upper left')

    axs[0].set_title('Line of Sight Linear Space\nGalactocentric Coordinates')
    axs[0].set_ylabel('[kpc]')

    axs[1].set_ylabel('[deg]')

    axs[2].set_ylabel('[deg]')
    axs[2].set_xlabel('Galactic Coordinate Distance [kpc]')

    axs[0].set_ylim(0, axs[0].set_ylim()[1])

    ##################################################################################
    #Denisty along Line of Sight
    ##################################################################################
    fig3, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4),
                                   sharey=True, gridspec_kw={'wspace': 0})
    fig_arr.append(fig3)
    # 
    ax1.semilogy(galactic_lin_distance,
                 rho_lin, '.',c='k')
    ax1.set_ylabel('DENSITY Msun*pc^-3')
    ax1.set_xlabel('Galactic Coordinate Distance [kpc]')
    ax1.set_xlim(0, ax1.set_xlim()[1])
    ax2.semilogy(galactocen_lin_spherical_distance,
                 rho_lin, '.',c='k')
    ax2.set_xlabel('Galactocentric Distance [kpc]')
    ax2.set_xlim(0, ax2.set_xlim()[1])

    plt.suptitle('Milky Way Density along Line-of-Sight')

    ##################################################################################
    #Cumulative Projected Desnity
    ##################################################################################
    fig4 = plt.figure(figsize=(10,7))
    fig_arr.append(fig4)

    plt.plot(galactic_lin_distance,
             np.cumsum(rho_lin) * r_max / n_lin,
             '.', label='Density',c='k')
    plt.legend(loc='best')
    plt.ylabel('Cumulative Projected Density} [M_sun*pc^{-2}]$')
    plt.xlabel('Galactic Distance [kpc]')
    plt.xlim(0, plt.xlim()[1])

    ##################################################################################
    #Discrete CDF
    ##################################################################################
    fig5 = plt.figure(figsize=(10,7))
    fig_arr.append(fig5)
    plt.plot(cdf_los, lw=3)
    plt.ylabel('CDF')
    plt.xlabel('Line-of-Sight Linear Space Coordinate Index')

    ##################################################################################
    #Interpolated CDF
    ##################################################################################
    fig6 = plt.figure(figsize=(10,7))
    fig_arr.append(fig6)
    x = np.linspace(0,1,100)
    y = f_cdf_d(x)
    plt.plot(y, x, lw=3)
    plt.ylabel('CDF')
    plt.xlabel('Galactic Coordinate Distance along LOS [kpc]')

    ##################################################################################
    #Cylindrical distribution of PBHs
    #Infered Radii Distribution
    #Galactic Distance Distribution
    ##################################################################################
    fig7, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(14,4))
    fig_arr.append(fig7)
    # Create a plot to look at the projected cylindrical distribution of PBH
    ax1.plot(x_cyl * 1000, y_cyl * 1000, ',')
    ax1.set_aspect('equal', 'datalim')
    ax1.set_ylabel('Cylindrical LOS Samples l-axis [pc]')
    ax1.set_xlabel('Cylindrical LOS Samples b-axis [pc]')
    # Create a histogram of the infered radii to make sure that it shows the expected linear trend
    bins = 20
    ax2.hist(r_cyl * 1000, bins=bins, label='Sampled Number')
    ax2.plot((0, r_proj_los_cyl * 1000), (0, 2 * n_pbh / bins), c='C2', label='Expected Number')
    ax2.legend(loc='upper left')
    ax2.set_ylabel('Number of PBH per LOS Cylindrical Annuli')
    ax2.set_xlabel('Cylindrical Radius [pc]')

    # Create a histogram of the galactic distance distribution
    bins=100
    n_per_bin, _, _ = ax3.hist(d_galac, bins=bins, label='Sampled Number')
    ax3.plot(galactic_lin_distance,
             rho_lin * area_proj_los_cyl / pbh_mass * r_max / bins * 1000**3,c='C2', label='Approximate Expectation')
    ax3.legend(loc='best')
    ax3.set_ylabel('Number of PBH')
    ax3.set_xlabel('Galactic Distance [kpc]')

    fig7.tight_layout(pad=2.0)

    plt.suptitle('Distribution of PBHs in Cylindrical Line-of-Sight Tube (i.e. not light cone)')

    ##################################################################################
    #Checking Light Cone Boundaries
    ##################################################################################
    fig8, (ax1, ax2) = plt.subplots(2, 1, figsize=(5,8),
                                   sharex=True, gridspec_kw={'hspace': 0})
    fig_arr.append(fig8)
    ax1.plot(d_galac[~mask_obs_cone], r_cyl[~mask_obs_cone] * 1000,
             '.', alpha=0.1, label='Outside Light Cone')
    ax1.plot(d_galac[mask_obs_cone], r_cyl[mask_obs_cone] * 1000,
             '.', alpha=0.1, label='Inside Light Cone')
    ax1.plot((0,r_max),
             (0, field_of_view_radius * np.pi / 180 * r_max *1000),
             label='Light Cone Boundry')
    ax1.legend(loc='upper left')
    ax1.set_ylabel('Radius on LOS Cylinder [pc]')

    ax2.plot(d_galac[~mask_obs_cone], np.sqrt((l_galac[~mask_obs_cone] - l_radian)**2 + 
                                     (b_galac[~mask_obs_cone] - b_radian)**2) * 180/np.pi,
             '.', alpha=0.1, label='Outside Light Cone' )
    ax2.plot(d_galac[mask_obs_cone], np.sqrt((l_galac[mask_obs_cone] - l_radian)**2 +
                                    (b_galac[mask_obs_cone] - b_radian)**2) * 180/np.pi,
             '.', alpha=0.1, label='Inside Light Cone')
    ax2.plot((0,r_max),
             (field_of_view_radius, field_of_view_radius), label='Light Cone Boundry')
    ax2.set_ylim(0,5*field_of_view_radius)
    ax2.legend(loc='best')
    ax2.set_ylabel('Field of View Radius [deg]')
    ax2.set_xlabel('Galactic Coordinate Distance [kpc]')
    
    ##################################################################################
    #Saving all plots to PDF
    ##################################################################################
    pp = PdfPages(output_root+'_pbh_diagnostic_plots.pdf')
    for fig in fig_arr:
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
