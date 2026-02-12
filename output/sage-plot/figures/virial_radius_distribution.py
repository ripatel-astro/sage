#!/usr/bin/env python

"""
SAGE Virial Radius Distribution

This module generates a histogram showing the distribution of virial radii for SAGE galaxy data.
"""

import os
import random

import matplotlib.pyplot as plt
import numpy as np
from figures import (
    AXIS_LABEL_SIZE,
    IN_FIGURE_TEXT_SIZE,
    LEGEND_FONT_SIZE,
    get_stellar_mass_label,
    setup_legend,
    setup_plot_fonts,
)
from matplotlib.ticker import MultipleLocator


def plot(
    galaxies,
    volume,
    metadata,
    params,
    output_dir="plots",
    output_format=".png",
    verbose=False,
):
    """
    Create a velocity dispersion histogram.

    Args:
        galaxies: Galaxy data as a numpy recarray
        volume: Simulation volume in (Mpc/h)^3
        metadata: Dictionary with additional metadata
        params: Dictionary with SAGE parameters
        output_dir: Output directory for the plot
        output_format: File format for the output

    Returns:
        Path to the saved plot file
    """
    # Set random seed for reproducibility when sampling points
    random.seed(2222)

    # Set up the figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Apply consistent font settings
    setup_plot_fonts(ax)

    # Extract necessary metadata
    hubble_h = metadata["hubble_h"]

    # Maximum number of points to plot (for better performance and readability)
    dilute = 7500

    # Filter for valid galaxies 
    w = np.where((galaxies.StellarMass > 0))[0] 

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for virial radius plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for virial radius plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"VirialRadiusDistribution{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    #if len(w) > dilute:
    #    w = random.sample(list(w), dilute)
    BulgeMass = galaxies.BulgeMass[w]
    TotalMass = galaxies.StellarMass[w]
    BTRatio = BulgeMass/TotalMass

    SFR = galaxies.SfrDisk[w] + galaxies.SfrBulge[w]
    Stellar_Mass = galaxies.StellarMass[w] * 1.0e10 / hubble_h
    SSFR = SFR / Stellar_Mass

    SSFR_new = []
    for galaxy in SSFR:
        if galaxy == 0:
            SSFR_new.append(galaxy+1e-14)
        else:
            SSFR_new.append(galaxy)

    logSSFR = np.log10(SSFR_new)

    # Print some debug information if verbose mode is enabled
    #if verbose:      
    #    print(f"Velocity Dispersion Distribution plot debug:")
    #    print(f"Number of galaxies plotted: {len(w)}")
    #    print(f"  Galaxy Mass range: {min(w):.2f} to {max(w):.2f}")
    #print(BTRatio)

    #    hist, bin_edges = np.histogram(Rvir, 60)
    #    max_counts_index = np.argmax(hist)
    #    bin_start = bin_edges[max_counts_index]
    #    bin_end = bin_edges[max_counts_index + 1]
    #    print("Max count range:", bin_start, bin_end)


    # def normalizing(param_hist, Nbins_global, dmin, dmax):
    
    #     binwidth_cen = float((dmax-dmin)/Nbins_global)
    #     bins_cen = np.linspace(dmin, dmax, Nbins_global+1)

    #     Ntotal_cen = float(len(param_hist))
    #     counts, edge = np.histogram(param_hist, bins_cen)
    #     hist_cen = counts/(binwidth_cen*Ntotal_cen)   #to have the same probability for different width sizes, divide by delta d which is your width length 
    #     dM = edge[1] - edge[0] 
    #     bins = bins_cen[0:-1] + dM/2.

    #     return hist_cen, bins

    # Plot the galaxy data
    ax.hist(
         logSSFR,
         60,
         color = "cornflowerblue",
         edgecolor = "w"
    )

    # Nbins_global = 50
    # dmin = 0.0
    # dmax = 1.0

    # hist, bins = normalizing(BTRatio, Nbins_global, dmin, dmax )

    # ax.step(bins, hist, where='post', color= 'k', label='All')

    # Customize the plot
    ax.set_xlabel(r"Galaxy Virial Radius ($Mpc$)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Virial Radius Count", fontsize=AXIS_LABEL_SIZE)



    # Save the figure, ensuring the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Warning: Could not create output directory {output_dir}: {e}")
        # Try to use a subdirectory of the current directory as fallback
        output_dir = "./plots"
        os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"VirialRadiusDistribution{output_format}")
    if verbose:
        print(f"Saving Virial Radius Distribution plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path

# matplotlib.rcParams['axes.linewidth'] = 1.5
# matplotlib.rcParams['lines.linewidth'] = 1.5
# figsize = (8,8)
# size_label = 30
# size_legend = 15
# size_ticks = 20
# dot_size = 11



# def normalizing(param_hist, Nbins_global, dmin, dmax):
    
#     binwidth_cen = float((dmax-dmin)/Nbins_global)
#     bins_cen = np.linspace(dmin, dmax, Nbins_global+1)

#     Ntotal_cen = float(len(param_hist))
#     counts, edge = np.histogram(param_hist, bins_cen)
#     hist_cen = counts/(binwidth_cen*Ntotal_cen)   #to have the same probability for different width sizes, divide by delta d which is your width length 
#     dM = edge[1] - edge[0] 
#     bins = bins_cen[0:-1] + dM/2.

#     return hist_cen, bins

# Nbins_global = 50
# dmin = 0.0
# dmax = 1.0


# hist, bins = normalizing(BTRatio, Nbins_global, dmin, dmax )


# plt.clf()
# plt.close()
# fig, ax = plt.subplots(1, figsize=figsize, facecolor='white')
# ax.set_facecolor('white')

# ax.step( bins, hist, where='post', color= 'k', label='All' )

# # ax.set_title(r'$f_{\rm Edd} = 0.01-1$', fontsize=size_label)

# ax.set_xlabel(r'$a$', fontsize=size_label)
# ax.set_ylabel( r'$dP(a)/d a$', fontsize=size_label )
# ax.set_yscale('log')

# ax.set_xlim(0.0, 1.0)
# # ax.set_ylim(-0.1, 3.1)
# # ax.set_xticks( [ -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0 ] )
# # ax.set_yticks( [ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 ] )
# # ax.set_xticklabels( [ -12.0, -11.5, -11.0, -10.5, -10.0, -9.5, -9.0 ], fontsize=size_ticks )
# # ax.set_yticklabels( [ 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0 ], fontsize=size_ticks )

# ax.legend(frameon=False, framealpha=1, borderpad=1, labelspacing=0.3, 
#     edgecolor='k', fontsize=size_legend, facecolor='white', loc=2 )

# #plt.legend( loc=1 )
# #plt.tight_layout()

# plt.savefig('/Users/pater32/Documents/Research/SAGEModel/sage/output/results/millennium/plots/BTRatioHist.png')
# plt.show()