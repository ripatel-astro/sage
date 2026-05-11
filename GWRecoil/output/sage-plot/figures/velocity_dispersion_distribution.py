#!/usr/bin/env python

"""
SAGE Velocity Dispersion Distribution

This module generates a histogram showing the distribution of velocity dispersions for SAGE galaxy data.
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
    w = np.where((galaxies.Mvir > 0))[0] 

    # Check if we have any galaxies to plot
    if len(w) == 0:
        print("No suitable galaxies found for velocity dispersion plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for velocity dispersion plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"VelocityDispersionDistribution{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    if len(w) > dilute:
        w = random.sample(list(w), dilute)
    VelDisp = galaxies.VelDisp[w]

    # Print some debug information if verbose mode is enabled
    if verbose:      
        print(f"Velocity Dispersion Distribution plot debug:")
        print(f"Number of galaxies plotted: {len(w)}")
        print(f"  Galaxy Mass range: {min(w):.2f} to {max(w):.2f}")
        print(f"  Velocity Dispersion range: {min(VelDisp[w]):.3f} to {max(VelDisp[w]):.3f}")

        hist, bin_edges = np.histogram(VelDisp, 60)
        max_counts_index = np.argmax(hist)
        bin_start = bin_edges[max_counts_index]
        bin_end = bin_edges[max_counts_index + 1]
        print("Max count range:", bin_start, bin_end)
    

    # Plot the galaxy data
    ax.hist(
        VelDisp,
        60,
        color = "cornflowerblue",
        edgecolor = "w"
    )

    # Customize the plot
    ax.set_xlabel(r"Galaxy Velocity Dispersion", fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Velocity Dispersion Count", fontsize=AXIS_LABEL_SIZE)



    # Save the figure, ensuring the output directory exists
    try:
        os.makedirs(output_dir, exist_ok=True)
    except Exception as e:
        print(f"Warning: Could not create output directory {output_dir}: {e}")
        # Try to use a subdirectory of the current directory as fallback
        output_dir = "./plots"
        os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, f"VelocityDispersionDistribution{output_format}")
    if verbose:
        print(f"Saving Velocity Dispersion Distribution plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
