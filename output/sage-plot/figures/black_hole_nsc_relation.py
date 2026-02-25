#!/usr/bin/env python

"""
SAGE Black Hole mass relation

This module generates a plot showing the relationship between galaxy stellar mass vs. black hole mass for SAGE galaxy data.
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
    Create a galaxy stellar mass vs black hole mass chart.

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
    w = np.where((galaxies.StellarMass > 0) & (galaxies.BlackHoleMass > 0) & (galaxies.NSC > 0))[0] 

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for black hole - NSC mass plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for black hole - NSC mass plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"BHNSCRelation{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    if len(w) > dilute:
        w = random.sample(list(w), dilute)
    
    
    # Filter 
    StellarMass = galaxies.StellarMass[w]
    NSCMass = galaxies.NSC[w]
    BHMass = galaxies.BlackHoleMass[w] 
    BHNSCRatio = BHMass/NSCMass

    logStellarMass = np.log10(StellarMass*(10**10)/hubble_h)
    logNSC = np.log10(NSCMass*(10**10)/hubble_h)
    logBH = np.log10(BHMass*(10**10)/hubble_h)
    logBHNSCRatio = np.log10(BHNSCRatio)

    

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
    
    logNSC_early = logNSC[BulgeTotalRatio >= BulgeRatioThreshold]
    logStellarMass_early = logStellarMass[BulgeTotalRatio >= BulgeRatioThreshold]
    logBHNSC_early = logBHNSCRatio[BulgeTotalRatio >= BulgeRatioThreshold]
    logNSC_late = logNSC[BulgeTotalRatio < BulgeRatioThreshold]
    logStellarMass_late = logStellarMass[BulgeTotalRatio < BulgeRatioThreshold]
    logBHNSC_late = logBHNSCRatio[BulgeTotalRatio < BulgeRatioThreshold]
    

    # Print some debug information if verbose mode is enabled
    # if verbose:      
    print(f"BH-NSC Mass plot debug:")
    print(f"Number of galaxies plotted: {len(w)}")
    print(f"  Galaxy stellar mass mass range: {min(logStellarMass):.2f} to {max(logStellarMass):.2f}")
    print(f"  Galaxy NSC mass range: {min(logNSC):.2f} to {max(logNSC):.2f}")
    print(f"  BH range: {min(logBH):.3f} to {max(logBH):.3f}")
    print(f"    BH-NSC Ratio range: {min(logBHNSCRatio):.3f} to {max(logBHNSCRatio):.3f}")


      

    ## Plot the galaxy data
    # Plot a dashed line to show galaxies with equal BH Mass and NSC Mass
    plt.axhline(
        y=0,
        color="black",
        linestyle="dashed",
    )

    ax.scatter(
        logStellarMass_late,
        logBHNSC_late,
        marker="o",
        s=3,
        c="dodgerblue",
        alpha=0.5,
        label="Late-type galaxies",
    )

    ax.scatter(
        logStellarMass_early,
        logBHNSC_early,
        marker="o",
        s=3,
        c="firebrick",
        alpha=0.5,
        label="Early-type galaxies",
    )

    #ax.scatter(
    #    logStellarMass,
    #    logBHNSCRatio,
    #    marker="o",
    #    s=3,
    #    c="k",
    #    alpha=0.5,
    #    label="Model galaxies"
    #)

    

    # Customize the plot
    ax.set_xlabel(r"log($M_{\star}$)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"log($M_{BH}/M_{NSC}$)", fontsize=AXIS_LABEL_SIZE)

    # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(1))

    # Set axis limits - matching the original plot
    ax.set_xlim(8, 12.5)
    ax.set_ylim(-4, 4)

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

    output_path = os.path.join(output_dir, f"BHNSCRelationPlot{output_format}")
    if verbose:
        print(f"Saving BHNSCRelation plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
