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
    verbose = True
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
    fig, ax = plt.subplots(figsize=(10, 10))

    # Apply consistent font settings
    setup_plot_fonts(ax)

    # Extract necessary metadata
    hubble_h = metadata["hubble_h"]

    # Maximum number of points to plot (for better performance and readability)
    dilute = 7500

    # Filter for valid galaxies with non-zero stellar mass
    w = np.where((galaxies.StellarMass > 0))[0]
    has_bh = (galaxies.StellarMass > 0) & (galaxies.BlackHoleMass > 0)

    m = (galaxies.StellarMass > 0) & (galaxies.BHmergeCount > 0)
    has_bh_merged = has_bh & (galaxies.BHmergeCount > 0)

     

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

    # Define safe division function to handle cases where the denominator might be zero
    def safe_div(numerator, denominator):
        if denominator > 0:
            return numerator / denominator
        else:
            return 0.0
    
    
    # Calculate required quantities for all selected galaxies - safely handle zeros
    stellar_mass = np.log10(galaxies.StellarMass[w] * 1.0e10 / hubble_h)
    stellar_mass_bh = np.log10(galaxies.StellarMass[has_bh] * 1.0e10 / hubble_h)

    stellar_mass_merged = np.log10(galaxies.StellarMass[m] * 1.0e10 / hubble_h)
    stellar_mass_merged_bh = np.log10(galaxies.StellarMass[has_bh_merged] * 1.0e10 / hubble_h)

    # Define mass bins
    min_range = 9
    max_range = 11.5
    interval = 0.1
    nbins = int((max_range - min_range) / interval)
    mass_bins = np.arange(min_range, max_range, interval)

    # Arrays to store results
    mass = []  # Bin centers
    bh_occ_fraction = []  # BH occupation fraction for all galaxies
    bh_occ_merge_fraction = []  # BH occupation fraction for galaxies with merged BHs

    # Calculate fractions for each mass bin
    for i in range(nbins - 1):
        # All galaxies in this mass bin
        this_bin_mask = (stellar_mass >= mass_bins[i]) & (stellar_mass < mass_bins[i + 1])
        this_bin_count = np.sum(this_bin_mask)

        # Galaxies with black holes in this mass bin
        bh_mask = (stellar_mass_bh >= mass_bins[i]) & (stellar_mass_bh < mass_bins[i + 1])
        bh_count = np.sum(bh_mask)

        #Galaxies with merged black hole in this mass bin
        bh_merge_mask = (stellar_mass_merged >= mass_bins[i]) & (stellar_mass_merged < mass_bins[i + 1])
        bh_merge_count = np.sum(bh_merge_mask)

        #Galaxies with merged black holes that have a black hole in this mass bin
        bh_merge_retain_mask = (stellar_mass_merged_bh >= mass_bins[i]) & (stellar_mass_merged_bh < mass_bins[i + 1])
        bh_merge_retain_count = np.sum(bh_merge_retain_mask)

        #Calculate occupation fraction
        bh_occ_fraction.append(safe_div(bh_count, this_bin_count))
        bh_occ_merge_fraction.append(safe_div(bh_merge_retain_count, bh_merge_count))

        # Store the mass bin center
        mass.append((mass_bins[i] + mass_bins[i + 1]) / 2.0)

    # Convert to numpy arrays
    mass = np.array(mass)
    bh_occ_fraction = np.array(bh_occ_fraction)
    bh_occ_merge_fraction = np.array(bh_occ_merge_fraction)


    if params["BHRecoilOn"] == 0:
        np.save('bhfraction_norecoil.npy', bh_occ_fraction)
        np.save('bhfraction_merge_norecoil.npy', bh_occ_merge_fraction)
    if params["BHRecoilOn"] == 1 or params["BHRecoilOn"] == 2 or params["BHRecoilOn"] == 3 or params["BHRecoilOn"] == 4:
        bh_occ_fraction_norecoil = np.load('bhfraction_norecoil.npy')
        bh_occ_merge_fraction_norecoil = np.load('bhfraction_merge_norecoil.npy')

        f_diff = []
        f_diff_merge = []
        for i in range(len(mass)):
            f_diff_i = safe_div((bh_occ_fraction[i] - bh_occ_fraction_norecoil[i]), bh_occ_fraction_norecoil[i])
            f_diff.append(f_diff_i)
            f_diff_merge_i = safe_div((bh_occ_merge_fraction[i] - bh_occ_merge_fraction_norecoil[i]), bh_occ_merge_fraction_norecoil[i])
            f_diff_merge.append(f_diff_merge_i)
        f_diff = np.array(f_diff)
        f_diff_merge = np.array(f_diff_merge)

    # Print some debug information if verbose mode is enabled
    ejection_current = np.where((galaxies.BHejectFlag[w] > 0) & (galaxies.BlackHoleMass[w]==0))[0]
    if verbose:      
        print(f"BH Occupation fraction plot debug:")
        print(f"Number of galaxies plotted: {len(w)}")
        print(f"Number of galaxies with black holes: {np.sum(has_bh)}")
        print(f"Number of galaxies with merged black holes: {np.sum(has_bh_merged)}")
        print(f"Number of galaxies with ejected black holes (no refill/reseed): {len(ejection_current)}")
        # print(f"  BH Mass range: {min(logBH):.3f} to {max(logBH):.3f}")
        # print(f"    Number of late-type galaxies with black holes: {len(EarlyBHGalaxy)}")
        # print(f"    Number of early-type galaxies with black holes: {len(LateBHGalaxy)}")

      

    # Plot the galaxy data
    if params["BHRecoilOn"] == 0:
        ax.plot(
            mass,
            bh_occ_fraction,
            color="black",
            label="All galaxies",
        )
        ax.plot(
            mass,
            bh_occ_merge_fraction,
            color="purple",
            label="Galaxies with Merged BHs",
            linestyle = "--",
            linewidth = 2
        )
    if params["BHRecoilOn"] == 1 or params["BHRecoilOn"] == 2 or params["BHRecoilOn"] == 3 or params["BHRecoilOn"] == 4:
        ax.plot(
            mass,
            bh_occ_fraction_norecoil,
            color="black",
            linestyle = "--",
            linewidth = 2,
            label="All Galaxies (No Recoil)"
        )
        ax.plot(
            mass,
            bh_occ_merge_fraction_norecoil,
            color="grey",
            linestyle = "--",
            linewidth = 2,
            label="Galaxies with BH-BH Mergers (No Recoil)"
        )  
        ax.plot(
            mass,
            bh_occ_fraction,
            color="plum",
            label="All Galaxies (with Recoil)",
        )
        ax.plot(
            mass,
            bh_occ_merge_fraction,
            color="purple",
            label="Galaxies with BH-BH Mergers (with Recoil)",
        )
        ax.plot(
            mass,
            f_diff,
            color="magenta",
            label="Fractional Difference in BH Occupation (with vs no recoil)"
        )
        # ax.plot(
        #     mass,
        #     f_diff_merge,
        #     color="magenta",
        # )  

    

    # Customize the plot
    ax.set_title("Black Hole Occupation Fraction", fontsize=AXIS_LABEL_SIZE + 2)
    ax.set_xlabel(get_stellar_mass_label(), fontsize=AXIS_LABEL_SIZE)
    ax.set_ylabel(r"Fraction of Galaxies with BHs", fontsize=AXIS_LABEL_SIZE)

    # # Set the x and y axis minor ticks
    ax.xaxis.set_minor_locator(MultipleLocator(1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    # Set axis limits - matching the original plot
    ax.set_xlim(min_range, max_range)
    ax.set_ylim(-0.35, 1.05)

    # Add consistently styled legend
    setup_legend(ax, loc="lower left")

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
