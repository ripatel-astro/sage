#!/usr/bin/env python

"""
SAGE Black Hole mass to galaxy type relation

This module generates a plot showing the relationship between galaxy stellar mass vs. bulge-to-total mass ratio for SAGE galaxy data.
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
    verbose = True
    """
    Create a black hole mass vs bulge-to-total mass ratio chart.

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
    w = np.where((galaxies.StellarMass > 0) & (galaxies.BlackHoleMass > 0))[0] 

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for black hole mass - galaxy type plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for black hole mass - galaxy type plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"BHMassGalTypePlot{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    if len(w) > dilute:
        w = random.sample(list(w), dilute)
    
    
    # Calculate log stellar mass and log black hole mass for the selected galaxies
    m = galaxies.BHmergeCount[w] > 0
    e = galaxies.BHejectFlag[w] > 0
    StellarMass = galaxies.StellarMass[w]
    BHMass = galaxies.BlackHoleMass[w] 

    logStellarMass = np.log10(StellarMass*(10**10)/hubble_h)
    logBH = np.log10(BHMass*(10**10)/hubble_h)
    
    #Calculate B/T Ratio
    BulgeMass = galaxies.BulgeMass[w]
    BulgeTotalRatio = BulgeMass/StellarMass
    BulgeRatioThreshold = .7
    

    # Print some debug information if verbose mode is enabled
    if verbose:      
        print(f"BH Mass plot debug:")
        print(f"Number of galaxies plotted: {len(w)}")
        print(f"  Galaxy stellar mass range: {min(logStellarMass):.2f} to {max(logStellarMass):.2f}")
        print(f"  BH range: {min(logBH):.3f} to {max(logBH):.3f}")


      

    # Plot the galaxy data
    ax.scatter(
        logBH,
        BulgeTotalRatio,
        marker="o",
        s=3,
        c="k",
        alpha = 0.5,
        label="Model galaxies"
    )
    ax.scatter(
        logBH[m],
        BulgeTotalRatio[m],
        marker="o",
        s=3,
        c="purple",
        label="Model galaxies with merged BHs"
    )
    ax.scatter(
        logBH[e],
        BulgeTotalRatio[e],
        marker="o",
        s=3,
        c="red",
        label="Model galaxies with ejected BHs"
    )
    

    # Customize the plot
    ax.set_xlabel(r"log$_{10}$ M$_{\rm BH}$ (M$_{\odot}$)", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"B/T", fontsize=AXIS_LABEL_SIZE)

    # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(.05))

    # Set axis limits - matching the original plot
    ax.set_xlim(0, 10)
    ax.set_ylim(-0.05, 1.05)

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

    output_path = os.path.join(output_dir, f"BHMassGalTypePlot{output_format}")
    if verbose:
        print(f"Saving BHMassGalType plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
