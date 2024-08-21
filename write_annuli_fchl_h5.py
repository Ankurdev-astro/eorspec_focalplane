"""
This script processes the annuli data for each channel and each step and writes the results to HDF5 files.
It iterates over frequency channels and FPI steps to generate separate directories and files
for each channel. The script reads a detector table from an HDF5 file, 
processes the annuli data using the `process_annuli_fchl` function, 
and writes the results to separate HDF5 files for each frequency channel and FPI step.

Purpose:
The purpose of this script is to organize and store the processed annuli data in a 
structured manner for further analysis and visualization in the EoR-Spec project.
"""

#Imports
from astropy.table import QTable
import os
from prepare_annuli import process_annuli_fchl
from fpi_step import process_steps, getall_freq_chl

# Base directory for storing HDF5 files
annuli_dir = './fpi_data/annuli_data/'
fpi_data_dir = './fpi_data/' 

# Get the list of fpi steps
fpi_steps = process_steps()
freq_channels = getall_freq_chl()
freq_channels.sort()

freq_channels_h5 = os.path.join(fpi_data_dir, "fchl_h5")
if not os.path.exists(freq_channels_h5):
    os.makedirs(freq_channels_h5)
    print(f"Created directory: {freq_channels_h5}")
        
for i, freq_chl in enumerate(freq_channels):
    # if i == 2:
    #     break
    print("\n","="*50,"\n")
    print(f"Processing freq channel {freq_chl} GHz...")
    
    fchl_h5 = os.path.join(freq_channels_h5, f"chnl_{freq_chl}")
    if not os.path.exists(fchl_h5):
        os.makedirs(fchl_h5)
        print(f"Created directory: {fchl_h5}")
    

    for j, step in enumerate(fpi_steps):           
        print("\t"*2,"-"*20)
        print("\t"*2,f"on FPI {step} ...")
        fpistep_infotxt = os.path.join(annuli_dir, f"annulus_results_{step}.txt")   
        
        # Load the detector table
        dets_table_full = QTable.read('./eorspec_dettable.h5', path='dettable_stack')
        
        dets_table = process_annuli_fchl(dets_table_full, 
                                        fpistep_infotxt, freq_chl, step)
                
        if dets_table is not None:
            ndets = len(dets_table)
            if ndets == 0:
                continue
            fchl_step_h5 = os.path.join(fchl_h5, f"{step}")
            if not os.path.exists(fchl_step_h5):
                os.makedirs(fchl_step_h5)
                print(f"Created directory: {fchl_step_h5}")
            
            h5file_path = os.path.join(fchl_step_h5, f'f{freq_chl}_{step}_d{ndets}_dettable.h5')
            dets_table.write(h5file_path, path=f'{freq_chl}_{step}', 
                                        serialize_meta=True,overwrite=True)
            print(f"Wrote h5 file for freq channel {freq_chl} GHz", 
                f"and FPI {step} in {h5file_path}...", 
                "\n")
            
print("Done!")