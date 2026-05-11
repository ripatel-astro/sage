


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


def load_snapshot(param_file, snapshot):
    """
    Return (galaxies, volume, metadata) for one snapshot.
    """

    # Parse params once
    params = sage_plot.SAGEParameters(param_file)

    # Build snapshot -> redshift mapper
    mapper = sage_plot.SnapshotRedshiftMapper(
        param_file,
        params.params,
        params["OutputDir"],
    )

    # Get redshift info
    redshift = mapper.get_redshift(snapshot)

    # Construct same model path used in sage-plot.py
    model_base = os.path.join(
        params["OutputDir"],
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