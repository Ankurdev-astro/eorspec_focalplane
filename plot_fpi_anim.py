# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation
# from matplotlib.colors import BoundaryNorm
# import numpy as np

# import toast.qarray as qa


import os, h5py
from astropy.table import QTable
import astropy.units as u
from toast.instrument import Focalplane
from plotting_func import animate_eorspec_annuli

fpi_data_dir = './fpi_data/' 
target_dir = os.path.join(fpi_data_dir, "fpisteps_h5")
plots_dir = os.path.join(fpi_data_dir, "fpi_plots")
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

fp_list = []

#Loop through all the h5 files in the target directory
for h5file in os.listdir(target_dir):
    if not h5file.endswith('.h5'):
        continue
    h5file_path = os.path.join(target_dir, h5file)
    print(f"Reading h5 file: {h5file} ...")
    # Open the HDF5 file
    with h5py.File(h5file_path, 'r') as dets_h5file:
        # Get the first path (assuming there's only one)
        path = list(dets_h5file.keys())[0]
        print(f"Processing Path: {path}")
        
    # Read the detector table from the HDF5 file        
    dets_table = QTable.read(h5file_path, path)

    # Plotting the annuli
    sample_rate = 244 * u.Hz
    width= 1.3 * u.degree
    #width of plot set to 1.3 degrees fixed in plot_eorspec_annuli
    
    fp_fpi_step =  Focalplane(
                    detector_data=dets_table,
                    sample_rate=sample_rate,
                    field_of_view=1.1 * width,
                    )
    fp_list.append(fp_fpi_step)

anim_file = os.path.join(plots_dir, "EoR-Spec_anim_FPI_01fps.gif")
animate_eorspec_annuli(focalplane_list=fp_list, outfile=anim_file, 
                       label_step=True, interval=1000)

    