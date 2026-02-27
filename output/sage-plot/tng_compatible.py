import h5py
import numpy as np
import os

# ---------- CONFIG ----------
input_dir = "./input/data/illustristng/tng50_4"      # folder with raw TNG50 tree files
output_dir = "./input/data/illustristng/tng50_4_sage"  # folder for SAGE-ready trees
h = 0.6774  # Hubble parameter in TNG
kpc_to_Mpc = 1.0 / 1000.0  # kpc -> Mpc

os.makedirs(output_dir, exist_ok=True)

# List all TNG HDF5 tree files
tree_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".hdf5")])

for fname in tree_files:
    infile = os.path.join(input_dir, fname)
    outfile = os.path.join(output_dir, fname.replace(".hdf5", "_sage.hdf5"))

    print(f"Converting {fname} → {os.path.basename(outfile)}")

    with h5py.File(infile, "r") as f_in, h5py.File(outfile, "w") as f_out:
        for tree_name in f_in.keys():
            g_in = f_in[tree_name]

            # Skip empty trees
            if g_in["SubhaloLen"].size == 0:
                continue

            g_out = f_out.create_group(tree_name)

            # ---------- Map datasets ----------
            g_out.create_dataset("Mvir", data=g_in["SubhaloMassType"][:,0])  # DM mass
            g_out.create_dataset("Len", data=g_in["SubhaloLen"][:])          # particle count
            g_out.create_dataset("Pos", data=g_in["SubhaloPos"][:] * kpc_to_Mpc * h)
            g_out.create_dataset("Vel", data=g_in["SubhaloVel"][:])
            g_out.create_dataset("Vmax", data=g_in["SubhaloVMax"][:])
            g_out.create_dataset("MostBoundID", data=g_in["SubhaloIDMostBound"][:])
            g_out.create_dataset("SnapNum", data=g_in["SnapNum"][:])
            g_out.create_dataset("FirstProgenitor", data=g_in["FirstProgenitor"][:])
            g_out.create_dataset("NextProgenitor", data=g_in["NextProgenitor"][:])
            g_out.create_dataset("Descendant", data=g_in["Descendant"][:])