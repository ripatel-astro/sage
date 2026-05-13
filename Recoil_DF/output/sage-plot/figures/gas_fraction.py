#!/usr/bin/env python

"""
SAGE Gas Fraction Plot

This module generates a plot showing the gas fraction vs. stellar mass for SAGE galaxy data.
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
    Create a gas fraction vs. stellar mass plot.

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

    # First filter for valid mass values - avoid division by zero
    valid_mass = ((galaxies.Type == 0) & (galaxies.StellarMass > 0.0) & (galaxies.StellarMass + galaxies.ColdGas > 0.0))
    m = (valid_mass & (galaxies.BHmergeCount > 0))
    r = (valid_mass & (galaxies.BHrefillCount > 0) & (galaxies.BHejectFlag == 0))
    e = (valid_mass & (galaxies.BHejectCount > 0))

    # Calculate ratio safely for valid galaxies
    # bulge_ratio = np.zeros_like(galaxies.StellarMass)
    # bulge_ratio_vm = galaxies.BulgeMass[valid_mass] / galaxies.StellarMass[valid_mass]
    # bulge_ratio_m = galaxies.BulgeMass[m] / galaxies.StellarMass[m]
    # bulge_ratio_r = galaxies.BulgeMass[r] / galaxies.StellarMass[r]
    

    # Now apply all filters
    # w = np.where(valid_mass & (bulge_ratio_vm > 0.1) & (bulge_ratio_vm < 0.5))[0]

    # Check if we have any galaxies to plot
    if len(galaxies[valid_mass]) == 0:
        print("No suitable galaxies found for gas fraction plot")
        # Create an empty plot with a message
        ax.text(
            0.5,
            0.5,
            "No suitable galaxies found for gas fraction plot",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            fontsize=IN_FIGURE_TEXT_SIZE,
        )

        # Save the figure
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"GasFraction{output_format}")
        plt.savefig(output_path)
        plt.close()
        return output_path

    # If we have too many galaxies, randomly sample a subset
    # if len(galaxies[valid_mass]) > dilute:
    #     w = random.sample(list(w), dilute)

    # Calculate gas fraction and convert stellar mass to log scale
    stellar_mass = np.log10(galaxies.StellarMass[valid_mass] * 1.0e10 / hubble_h)
    gas_fraction = galaxies.ColdGas[valid_mass] / (galaxies.StellarMass[valid_mass] + galaxies.ColdGas[valid_mass])

    stellar_mass_merge = np.log10(galaxies.StellarMass[m] * 1.0e10 / hubble_h)
    gas_fraction_merge = galaxies.ColdGas[m] / (galaxies.StellarMass[m] + galaxies.ColdGas[m])

    stellar_mass_refill = np.log10(galaxies.StellarMass[r] * 1.0e10 / hubble_h)
    gas_fraction_refill = galaxies.ColdGas[r] / (galaxies.StellarMass[r] + galaxies.ColdGas[r])

    stellar_mass_eject = np.log10(galaxies.StellarMass[e] * 1.0e10 / hubble_h)
    gas_fraction_eject = galaxies.ColdGas[e] / (galaxies.StellarMass[e] + galaxies.ColdGas[e])



    # Print some debug information if verbose mode is enabled
    if verbose:
        print(f"Gas Fraction plot debug:")
        print(f"  Number of galaxies plotted: {len(w)}")
        print(
            f"  Stellar mass range: {min(stellar_mass):.2f} to {max(stellar_mass):.2f}"
        )
        print(
            f"  Gas fraction range: {min(gas_fraction):.3f} to {max(gas_fraction):.3f}"
        )

    # Plot the galaxy data
    # ax.scatter(
    #     stellar_mass,
    #     gas_fraction,
    #     marker="o",
    #     s=1,
    #     c="k",
    #     alpha=0.5,
    #     label="Model Sb/c galaxies",
    # )

    ax.scatter(
        stellar_mass,
        gas_fraction,
        marker="o",
        s=1,
        c="grey",
        alpha=0.3,
        label="Model galaxies",
    )

    ax.scatter(
        stellar_mass_eject,
        gas_fraction_eject,
        marker="x",
        s=30,
        linewidths=1,
        c="firebrick",
        alpha=0.7,
        label="Model galaxies with ejected BHs",
    )

    ax.scatter(
        stellar_mass_merge,
        gas_fraction_merge,
        marker="o",
        s=5,
        c="royalblue",
        alpha=0.7,
        label="Model galaxies with BH-BH mergers",
    )

    ax.scatter(
        stellar_mass_refill,
        gas_fraction_refill,
        marker="o",
        s=5,
        c="magenta",
        alpha=0.7,
        label="Model galaxies with refilled BHs",
    )

    

    # Customize the plot
    ax.set_xlabel(get_stellar_mass_label(), fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Cold Mass / (Cold+Stellar Mass)", fontsize=AXIS_LABEL_SIZE)

    # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    # Set axis limits - matching the original plot
    ax.set_xlim(9.0, 12.0)
    ax.set_ylim(-0.1, 1.0)

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

    output_path = os.path.join(output_dir, f"GasFraction{output_format}")
    if verbose:
        print(f"Saving Gas Fraction plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
