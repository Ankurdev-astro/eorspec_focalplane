"""
This script reads detector data from an HDF5 file into an astropy QTable,
initializes a Focalplane object with specific parameters such as sample rate
and field of view, and then plots the focal plane annuli using the 
plot_eorspec_annuli function. The purpose of this script is to visualize the 
arrangement of detectors in the focal plane based on the provided data.
"""
import os, h5py
from astropy.table import QTable
from toast.instrument import Focalplane
import matplotlib.pyplot as plt
import astropy.units as u
from plotting_func import plot_eorspec_annuli


fpi_data_dir = './fpi_data/' 
target_dir = os.path.join(fpi_data_dir, "fpisteps_h5")
plots_dir = os.path.join("./fpi_data", "fpi_plots")
if not os.path.exists(plots_dir):
    os.makedirs(plots_dir)

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
    
    outfile_path = os.path.join(plots_dir, f"fp_eorspec_{path}.png")
    # outfile_path = os.path.join(plots_dir, f"fp_eorspec_{path}.pdf")
    fp_test =  Focalplane(
                    detector_data=dets_table,
                    sample_rate=sample_rate,
                    field_of_view=1.1 * width,
                    )

    plot_eorspec_annuli(
        focalplane=fp_test,
        label_step=True,
        outfile=outfile_path)

    plt.show()
    # break



