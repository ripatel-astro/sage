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
    w = np.where((galaxies.StellarMass > 0))[0] 
    

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

    # # If we have too many galaxies, randomly sample a subset
    # if len(w) > dilute:
    #     w = random.sample(list(w), dilute)
    
    # Filter for valid black hole masses and compute black hole mass ratio
    q = galaxies.BHMassRatio[w]
    if params['BHRecoilOn'] == 0:
        BHRecoil = np.zeros_like(q)
    BHRecoil = galaxies.BHRecoilVmag[w]
    ejected = galaxies.BHejectCount[w]
    merged = galaxies.BHmergeCount[w]

    m = np.where((galaxies.BHmergeType == 1) | (galaxies.BHmergeType == 4))[0]

    gal_BH = np.where((galaxies.StellarMass > 0) & (galaxies.BlackHoleMass > 0))[0]

    # Print some debug information if verbose mode is enabled
    verbose = True
    if verbose:      
        print(f"Black Hole Recoil Relation plot debug:")
        print(f"Number of galaxies with black holes plotted: {len(w)}")
        print(f" Number of galaxies with merged black holes: {len(np.where(merged > 0)[0])}")
        print(f" Number of total merged black holes (all history): {np.sum(merged)}")
        print(f"  Number of galaxies with ejected black holes: {len(np.where(ejected > 0)[0])}")
        print(f"  Number of ejected black holes (all history): {np.sum(ejected)}")
        print(f"  Max number of ejected black holes for one galaxy+progenitors: {np.max(ejected)}")
        print(f"  Black Hole Mass Ratio range: {min(q):.2f} to {max(q):.2f}")
        print(f"  Recoil range: {min(BHRecoil):.3f} to {max(BHRecoil):.3f}")
        print(f"  Position of subhalos: {galaxies.Pos[w][0:5]}")


    # Plot recoil lines
    A = 12000
    B = -0.93
    H = 7300
    K = 60000

    def calc_recoil_velocity(q, e, a_1, a_2, xi, phi_1, phi_2, theta):
        p = (q**2)/((1+q)**5)
        a_1_perp = a_1*np.sin(phi_1)
        a_2_perp = a_2*np.sin(phi_2)
        a_1_par = a_1*np.cos(phi_1)
        a_2_par = a_2*np.cos(phi_2)
        v_m = (A*p*(1-q))*(1+(B*(q/((1+q)**2)))) #velocity based on unequal mass
        v_perp =  (H*p)*(a_1_par - q*a_2_par) #velocity perpendicular to orbital angular momentum direction
        v_par = (K*p*np.cos(theta))*(a_1_perp - q*a_2_perp) #velocity parallel to orbital angular momentum direction
        #v_rec_vector = [v_m + v_perp*np.cos(xi), v_perp*np.sin(xi), v_par]
        v_rec = (1+e)*np.sqrt((v_m**2)+(v_perp**2)+(v_par**2)+2*v_m*v_perp*np.cos(xi))
        return v_rec

    num_sim = 500
    x = np.linspace(0, 1, 100)
    
    recoil_no_spin = calc_recoil_velocity(x, 0, 0, 0, 0, 0, 0, 0)
    
    a_1 = np.random.uniform(0, 1, num_sim)
    a_2 = np.random.uniform(0, 1, num_sim)
    xi = np.random.uniform(0, 2*np.pi, num_sim)
    phi_1 = np.random.uniform(0, np.pi, num_sim)
    phi_2 = np.random.uniform(0, np.pi, num_sim)

    recoil_max = []


    if params['BHRecoilOn'] == 1:
        recoil_max = recoil_no_spin
    if params['BHRecoilOn'] == 2:
        for ratio in x:
            recoils = []
            for n in range(num_sim):
                recoil = calc_recoil_velocity(ratio, 0, a_1[n], a_2[n], 0, 0, 0, 0)
                recoils.append(recoil)
            recoil_max.append(max(recoils))
    if params['BHRecoilOn'] == 3:
        for ratio in x:
            recoils = []
            for n in range(num_sim):
                recoil = calc_recoil_velocity(ratio, 0, -a_1[n], -a_2[n], 0, 0, 0, 0)
                recoils.append(recoil)
            recoil_max.append(max(recoils))
    if params['BHRecoilOn'] == 4:
        for ratio in x:
            recoils = []
            for n in range(num_sim):
                recoil = calc_recoil_velocity(ratio, 0, a_1[n], a_2[n], xi[n], phi_1[n], phi_2[n], 0)
                recoils.append(recoil)
            recoil_max.append(max(recoils))

    if params['BHRecoilOn'] == 1 or params['BHRecoilOn'] == 2 or params['BHRecoilOn'] == 3 or params['BHRecoilOn'] == 4:
        plt.plot(x, 
                recoil_no_spin, 
                color='black', 
                label='No Spin Recoil Relation',
                linewidth = 1
                )
        
        plt.plot(x, 
                recoil_max, 
                color='black', 
                label='Max Recoil Relation',
                linewidth = 1,
                linestyle='dashed'
                )

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

    output_path = os.path.join(output_dir, f"BlackHoleRecoilRelation{output_format}")
    if verbose:
        print(f"Saving Black Hole - Recoil Magnitude Relation plot to: {output_path}")
    plt.savefig(output_path)
    plt.close()

    return output_path
