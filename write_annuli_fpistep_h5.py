"""
This script writes the detector table for each FPI step in a separate HDF5 file.
The detector table is read from the eorspec_dettable.h5 file and processed for each
FPI step using the process_annuli_fpistep function. The resulting detector table is
then written to an HDF5 file for each FPI step in the fpisteps_h5 directory.

The purpose of this script is to store the detector data for each FPI step in a
separate HDF5 file for further analysis and visualization.
"""
#Imports
from astropy.table import QTable
import os
from prepare_annuli import process_annuli_fpistep
from fpi_step import process_steps

# Base directory for storing HDF5 files
annuli_dir = './fpi_data/annuli_data/'
fpi_data_dir = './fpi_data/' 


# Get the list of fpi steps
fpi_steps = process_steps()

for i, step in enumerate(fpi_steps):
    # if i == 2:
    #     break
    print("-"*50)
    print(f"Processing FPI {step} ...")
    fpistep_infotxt = os.path.join(annuli_dir, f"annulus_results_{step}.txt")

    # # Load the detector table
    dets_table_full = QTable.read('./eorspec_dettable.h5', path='dettable_stack')

    dets_table = process_annuli_fpistep(dets_table_full, fpistep_infotxt, step)
    
    fpi_step_h5file = os.path.join(fpi_data_dir, "fpisteps_h5")

    if not os.path.exists(fpi_step_h5file):
        os.makedirs(fpi_step_h5file)
        print(f"Created directory: {fpi_step_h5file}")
    
    h5file_path = os.path.join(fpi_step_h5file, f'{step}_dettable.h5')
    dets_table.write(h5file_path, path=f'{step}', 
                                    serialize_meta=True,overwrite=True)
    print(f"Wrote h5 file for all annuli and freq for FPI {step} in {h5file_path}...", "\n")

print("Done!")