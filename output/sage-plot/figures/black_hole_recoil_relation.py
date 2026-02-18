#!/usr/bin/env python

"""
SAGE Black Hole - GW Recoil Relation Plot

This module generates a plot showing the relationship between black hole mass ratio vs. GW recoil velocity for SAGE galaxy data.
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
    Create a black hole mass ratio vs recoil velocity chart.

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

    # Filter for valid galaxies with black hole masses (w) and with merged black holes (m)
    w = np.where((galaxies.BlackHoleMass > 0))[0] 
    m = np.where((galaxies.BHMassRatio > 0))[0]

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for black hole recoil plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for black hole recoil plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"BlackHoleRecoilRelation{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    if len(w) > dilute:
        w = random.sample(list(w), dilute)
    
    # Filter for valid black hole masses and compute black hole mass ratio
    q = galaxies.BHMassRatio[w]
    BHRecoil = galaxies.BHRecoilVmag[w]
    ejected = np.where((galaxies.BHEjectedMass > 0))[0]

    # Print some debug information if verbose mode is enabled
    #if verbose:      
    print(f"Black Hole Recoil Relation plot debug:")
    print(f"Number of black holes plotted: {len(w)}")
    print(f" Number of galaxies with merged black holes: {len(m)}")
    print(f"  Black Hole Mass Ratio range: {min(q):.2f} to {max(q):.2f}")
    print(f"  Recoil range: {min(BHRecoil):.3f} to {max(BHRecoil):.3f}")

    print('Number of ejected black holes', len(ejected))

    BHGal = np.where((galaxies.BlackHoleMass > 0))[0]
    print("Galaxies with black holes", len(BHGal))

    NoBHMR = np.where(galaxies.BHMassRatio == 0)[0]
    print(f" Number galaxies without merged black holes", len(NoBHMR))

    # Plot the galaxy data
    ax.scatter(
        q,
        BHRecoil,
        marker="o",
        s=1,
        c="k",
        alpha=0.5,
        label="Model black holes",
    )

    # Customize the plot
    ax.set_xlabel(r"Black Hole Mass Ratio", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Recoil Velocity (kms$^{-1}$)", fontsize=AXIS_LABEL_SIZE)

    # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(10))

    # Set axis limits - matching the original plot
    ax.set_xlim(0, 1)
    ax.set_ylim(-25, 200)

    # Add consistently styled legend
    setup_legend(ax, loc="upper right")

    # Save the figure, ensuring the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Warning: Could not create output directory {output_dir}: {e}")
        # Try to use a subdirectory of the current directory as fallback
        output_dir = "./plots"
        os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"BlackHoleRecoilRelation{output_format}")
    if verbose:
        print(f"Saving Black Hole - Recoil Magnitude Relation plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
