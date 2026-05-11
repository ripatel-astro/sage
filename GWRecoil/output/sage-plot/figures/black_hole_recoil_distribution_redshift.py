# #!/usr/bin/env python

# """
# SAGE Black Hole Recoil Distribution

# This module generates a histogram showing the distribution of ejected and refilled black holes for SAGE galaxy data.
# """

# import os
# import random

# import matplotlib.pyplot as plt
# import numpy as np
# from figures import (
#     AXIS_LABEL_SIZE,
#     IN_FIGURE_TEXT_SIZE,
#     LEGEND_FONT_SIZE,
#     get_stellar_mass_label,
#     setup_legend,
#     setup_plot_fonts,
# )
# from matplotlib.ticker import MultipleLocator


# def plot(
#     galaxies,
#     volume,
#     metadata,
#     params,
#     output_dir="plots",
#     output_format=".png",
#     verbose=False,
# ):
#     verbose = True
#     """
#     Create a black hole recoil histogram.

#     Args:
#         galaxies: Galaxy data as a numpy recarray
#         volume: Simulation volume in (Mpc/h)^3
#         metadata: Dictionary with additional metadata
#         params: Dictionary with SAGE parameters
#         output_dir: Output directory for the plot
#         output_format: File format for the output

#     Returns:
#         Path to the saved plot file
#     """
#     # Set random seed for reproducibility when sampling points
#     random.seed(2222)

#     # Set up the figure
#     fig, ax = plt.subplots(figsize=(8, 6))

#     # Apply consistent font settings
#     setup_plot_fonts(ax)

#     # Extract necessary metadata
#     hubble_h = metadata["hubble_h"]

#     # Maximum number of points to plot (for better performance and readability)
#     dilute = 7500

#     # Filter for valid galaxies 
#     w = np.where((galaxies.StellarMass > 0))   

#     # Check if we have any galaxies to plot
#     if len(w) == 0:
#         print("No suitable galaxies found for virial radius plot")
#         # Create an empty plot with a message
#         ax.text(
#             0.5,
#             0.5,
#             "No suitable galaxies found for virial radius plot",
#             horizontalalignment="center",
#             verticalalignment="center",
#             transform=ax.transAxes,
#             fontsize=IN_FIGURE_TEXT_SIZE,
#         )

#         # Save the figure
#         os.makedirs(output_dir, exist_ok=True)
#         output_path = os.path.join(output_dir, f"BHRecoilDistribution_redshift{output_format}")
#         plt.savefig(output_path)
#         plt.close()
#         return output_path

#     # Create galaxy stellar mass bins
#     StellarMass = galaxies.StellarMass[w]
#     logStellarMass = np.log10((StellarMass*1e10)/hubble_h)

#     bin_min = 8
#     bin_max = 12
#     bin_width = 0.1
#     bin_count = int((bin_max - bin_min) / bin_width)


    
#     # Get relevant black hole values
#     mergeCount = galaxies.BHmergeCount[w]
#     ejectedCount = galaxies.BHejectCount[w]
#     refillCount = galaxies.BHrefillCount[w]
#     reseedCount = galaxies.BHreseedCount[w]
    


#     # Print some debug information if verbose mode is enabled
#     if verbose:      
#         print(f"Recoil Distribution plot debug:")
#         print(f"Number of black holes plotted: {len(w)}")
#         print(f"Number of merged black holes: {sum(mergeCount)}")
#         print(f"Number of ejected black holes: {sum(ejectedCount)}")
#         print(f"Number of refilled black holes: {sum(refillCount)}")  
#         print(f"Number of reseeded black holes: {sum(reseedCount)}")  
    

    

#     # Plot the galaxy data
#     ax.hist(
#          logStellarMass,
#          weights = mergeCount,
#          bins = bin_count,
#          range = (bin_min, bin_max),
#          alpha = 0.7,
#          color = "black",
#          edgecolor = "w",
#          label = "Merged Black Holes"
#     )
#     if params["BHRecoilOn"] > 0:
#         ax.hist(
#             logStellarMass,
#             weights = ejectedCount,
#             bins = bin_count,
#             range = (bin_min, bin_max),
#             alpha = 0.7,
#             color = "firebrick",
#             edgecolor = "w",
#             label = "Ejected Black Holes"
#         )
#         ax.hist(
#             logStellarMass,
#             weights = refillCount,
#             bins = bin_count,
#             range = (bin_min, bin_max),
#             alpha = 0.7,
#             color = "magenta",
#             edgecolor = "w",
#             label = "Refilled Black Holes"
#         )
#         ax.hist(
#             logStellarMass,
#             weights = reseedCount,
#             bins = bin_count,
#             range = (bin_min, bin_max),
#             alpha = 0.7,
#             color = "blue",
#             edgecolor = "w",
#             label = "Reseed Black Holes"
#         )



#     # Customize the plot
#     ax.set_xlabel(get_stellar_mass_label(), fontsize=AXIS_LABEL_SIZE)
#     ax.set_ylabel(r"Black Hole Count", fontsize=AXIS_LABEL_SIZE)
#     ax.set_title(r"Black Hole Recoil Distribution (Redshift)", fontsize=AXIS_LABEL_SIZE)
#     setup_legend(ax, loc="upper left")



#     # Save the figure, ensuring the output directory exists
#     try:
#         os.makedirs(output_dir, exist_ok=True)
#     except Exception as e:
#         print(f"Warning: Could not create output directory {output_dir}: {e}")
#         # Try to use a subdirectory of the current directory as fallback
#         output_dir = "./plots"
#         os.makedirs(output_dir, exist_ok=True)

#     output_path = os.path.join(output_dir, f"BHRecoilDistribution_redshift{output_format}")
#     if verbose:
#         print(f"Saving Black Hole Recoil Distribution (redshift) plot to: {output_path}")
#     plt.savefig(output_path)
#     plt.close()

#     return output_path






