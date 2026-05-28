


import importlib.util

path = "/Users/pater32/Documents/GitHub/sage/GWRecoil/output/sage-plot/sage-plot.py"

spec = importlib.util.spec_from_file_location("sage_plot", path)
sage_plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sage_plot)

import os
from snapshot_redshift_mapper import SnapshotRedshiftMapper


class DummyArgs:
    verbose = False

sage_plot.args = DummyArgs()

def get_params(param_file):
    """
    Returns parameter inputs from simulation .par file
    """
    params = sage_plot.SAGEParameters(param_file)
    return params


def load_snapshot(param_file, snapshot, output_dir = None):
    """
    Return (galaxies, volume, metadata) for one snapshot.
    """

    # Parse params once
    params = sage_plot.SAGEParameters(param_file)

    # Make output directory 
    output_directory = output_dir if output_dir is not None else params['OutputDir']

    # Build snapshot -> redshift mapper
    mapper = sage_plot.SnapshotRedshiftMapper(
        param_file,
        params.params,
        output_directory,
    )

    # Get redshift info
    redshift = mapper.get_redshift(snapshot)

    # Construct same model path used in sage-plot.py
    model_base = os.path.join(
        output_directory,
        f"{params['FileNameGalaxies']}{mapper.get_redshift_str(snapshot)}"
    )

    # Load galaxies
    try: 
        galaxies, volume, metadata = sage_plot.read_galaxies(
            model_path=model_base,
            first_file=params["FirstFile"],
            last_file=params["LastFile"],
            params=params.params,
        )
        metadata["redshift"] = redshift

        return galaxies, volume, metadata
    
    except Exception as e:
        print(f'Error generating galaxies: {e}')


def load_variable(param_file, snapshot, variable, output_dir=None):
    """
    Returns array of one galaxy variable for one snapshot.
    """

    galaxies, volume, metadata = load_snapshot(param_file, snapshot, output_dir)

    if not hasattr(galaxies, variable):
        print(f"No attribute found: {variable}")
        return None

    return getattr(galaxies, variable)
   

def load_all_snapshots(param_file):
    """
    Returns dict: snapshot → (galaxies, volume, metadata)
    """

    params = sage_plot.SAGEParameters(param_file)

    mapper = sage_plot.SnapshotRedshiftMapper(
        param_file,
        params.params,
        params["OutputDir"],
    )

    data = {}

    for snap in mapper.get_all_snapshots():
        data[snap] = load_snapshot(param_file, snap)

    return data


def load_variable_all_snapshots(param_file, variable, output_dir=None):
    """
    Returns array of one galaxy variable for all snapshot.
    """

    params = sage_plot.SAGEParameters(param_file)

    mapper = sage_plot.SnapshotRedshiftMapper(
        param_file,
        params.params,
        params["OutputDir"],
    )

    data = {}
    for snap in mapper.get_all_snapshots():
        data[snap] = load_variable(param_file, snap, variable, output_dir)
    
    return data