#!/usr/bin/env python

"""
SAGE Black Hole occupation fraction

This module generates a plot showing the relationship between galaxy stellar mass vs. black hole occupation fraction for SAGE galaxy data.
"""

## ADD SCATTER INTO THIS PLOT 

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
    Create a galaxy stellar mass vs black hole occupation fraction chart.

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

    # Filter for valid galaxies with non-zero stellar mass and nuclear star clusters
    w = np.where((galaxies.StellarMass > 0))[0]
    has_bh = galaxies.BlackHoleMass[w] > 0
     

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for black hole occupation plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for black hole occupation plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"BHOccupation{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    ## COMMENTED OUT BECAUSE WE WANT ALL GALAXIES
    # If we have too many galaxies, randomly sample a subset
    #if len(w) > dilute:
    #    w = random.sample(list(w), dilute)
    
    
    # Add scatter to NSC galaxies
    StellarMass = galaxies.StellarMass[w]
    logStellarMass = np.log10(StellarMass*(10**10)/hubble_h)
    

    ## Separate galaxies by early-type and late-type (can use either SSFR or B/T)
    # Calculate specific star formation rate
    SFR = galaxies.SfrDisk[w] + galaxies.SfrBulge[w]
    Stellar_Mass = galaxies.StellarMass[w] * 1.0e10 / hubble_h
    SSFR = SFR / Stellar_Mass
    SSFRThreshold = -11.0
    
    #Calculate B/T Ratio
    BulgeMass = galaxies.BulgeMass[w]
    BulgeTotalRatio = BulgeMass/StellarMass
    BulgeRatioThreshold = .5
    
    EarlyGalaxy = logStellarMass[SSFR < 10**SSFRThreshold]
    EarlyBHGalaxy = logStellarMass[(SSFR < 10**SSFRThreshold) & (has_bh)]
    LateGalaxy = logStellarMass[SSFR >= 10**SSFRThreshold]
    LateBHGalaxy = logStellarMass[(SSFR >= 10**SSFRThreshold) & (has_bh)]
    BHGalaxy = logStellarMass[has_bh]

    bins = np.arange(0, 15, .5)
    EarlyGalaxies_count, edges = np.histogram(EarlyGalaxy, bins = bins)
    EarlyBHGalaxies_count, _ = np.histogram(EarlyBHGalaxy, bins = edges)
    LateGalaxies_count, edges_late = np.histogram(LateGalaxy, bins = edges)
    LateBHGalaxies_count, _ = np.histogram(LateBHGalaxy, bins = edges)
    AllGalaxies_count, _ = np.histogram(logStellarMass, bins = edges)
    BHGalaxies_count, _ = np.histogram(BHGalaxy, bins = edges)

    e = np.where(EarlyGalaxies_count > 0)[0]
    l = np.where(LateGalaxies_count > 0)[0]
    g = np.where(AllGalaxies_count > 0)[0]

    f_occ_early = np.zeros_like(EarlyGalaxies_count, dtype=float)
    f_occ_late  = np.zeros_like(LateGalaxies_count,  dtype=float)
    f_occ_all = np.zeros_like(AllGalaxies_count, dtype=float)
    f_occ_early[e] = EarlyBHGalaxies_count[e] / EarlyGalaxies_count[e]
    f_occ_late[l] =  LateBHGalaxies_count[l] / LateGalaxies_count[l]
    f_occ_all[g] = BHGalaxies_count[g] / AllGalaxies_count[g]

    # Print some debug information if verbose mode is enabled
    if verbose:      
        print(f"BH Occupation fraction plot debug:")
        print(f"Number of galaxies plotted: {len(w)}")
        print(f"  Galaxy stellar mass range: {min(logStellarMass):.2f} to {max(logStellarMass):.2f}")
        print(f"  BH Mass range: {min(logBH):.3f} to {max(logBH):.3f}")
        print(f"    Number of late-type galaxies with black holes: {len(EarlyBHGalaxy)}")
        print(f"    Number of early-type galaxies with black holes: {len(LateBHGalaxy)}")

      

    # Plot the galaxy data
    midpoints = .5 * (edges[1:] + edges[:-1])
    
    #ax.plot(
    #    midpoints,
    #    f_occ_early,
    #    marker="o",
    #    color = "firebrick",
    #    label="Early-type galaxies",
    #)
    #ax.plot(
    #    midpoints,
    #    f_occ_late,
    #    marker="o",
    #    color="dodgerblue",
    #    label="Late-type galaxies",
    #)

    ax.plot(
        midpoints,
        f_occ_all,
        marker="o",
        color="black",
        label="All galaxies",
    )
    

    # Customize the plot
    ax.set_xlabel(r"log($M_{\star}$)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Fraction of Galaxies with BHs", fontsize=AXIS_LABEL_SIZE)

    # # Set the x and y axis minor ticks
    # ax.xaxis.set_minor_locator(MultipleLocator(1))
    # ax.yaxis.set_minor_locator(MultipleLocator(1))

    # Set axis limits - matching the original plot
    ax.set_xlim(8, 11.25)
    ax.set_ylim(0, 1)

    # Add consistently styled legend
    setup_legend(ax, loc="upper left")

    # Save the figure, ensuring the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Warning: Could not create output directory {output_dir}: {e}")
        # Try to use a subdirectory of the current directory as fallback
        output_dir = "./plots"
        os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"BHOccupationPlot{output_format}")
    if verbose:
        print(f"Saving BH Occupation plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
