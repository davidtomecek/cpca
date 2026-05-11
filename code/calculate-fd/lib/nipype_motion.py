from nipype.algorithms.confounds import FramewiseDisplacement
import os
import pandas as pd

def fd(motion_params_file, output_file):

    # Initialize the interface
    fd = FramewiseDisplacement()

    # Inputs
    fd.inputs.in_file = motion_params_file  # Path to your SPM motion file
    fd.inputs.parameter_source = 'SPM'   # Tells Nipype to expect 6 columns (mm, rad)
    fd.inputs.save_plot = False          # Set to True if you want a visual .png as well

    # Manually specify the path and name of the output text file
    fd.inputs.out_file = output_file

    # Run it
    res = fd.run()

    # The FD time series is stored in 'out_file'
    print(f"FD time series saved at: {res.outputs.out_file}")

def create_completion_table(records):
    """
    Creates a summary table from a list of results.
    
    Args:
        records (list): A list of dictionaries, e.g.:
                        [{'subject': 'sub-01', 'visit': 'V1', 'status': 1}, ...]
                        Only include records where the test was completed.
    
    Returns:
        pd.DataFrame: A table with subjects as rows and visits as columns.
    """
    # Create the DataFrame
    df = pd.DataFrame(records)
    
    # Pivot the table:
    # - index: unique subjects
    # - columns: unique visits
    # - values: the completion flag (1)
    return df.pivot(index='subject', columns='visit', values='status')
