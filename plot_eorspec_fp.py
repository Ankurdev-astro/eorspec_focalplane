"""
This script reads detector data from an HDF5 file into an astropy QTable,
initializes a Focalplane object with specific parameters such as sample rate
and field of view, and then plots the focal plane using the 
plot_focalplane_eorspec function. The purpose of this script is to visualize the 
arrangement of detectors in the focal plane based on the provided data.
"""
from astropy.table import QTable
from plotting_func import plot_focalplane_eorspec
from toast.instrument import Focalplane
import matplotlib.pyplot as plt
import astropy.units as u

# Read the Astropy Qtable from the HDF5 file
# hf_fulltable_file = './test_dir/eorspec_dettable.h5'
hf_fulltable_file = './eorspec_dettable.h5'
# QTable includes the Quantities as well
dettable_stack = QTable.read(hf_fulltable_file, path='dettable_stack')

#print(dettable_stack[:10])

sample_rate = 244 * u.Hz
width= 1.3 * u.degree
fp_test =  Focalplane(
                detector_data=dettable_stack,
                sample_rate=sample_rate,
                field_of_view=1.1 * width,
                )

plot_focalplane_eorspec(
    focalplane=fp_test,
    width=width,
    height=width,
    show_labels=False,
    outfile="test_dir/fp_eorspec_rot2_2024_v4.pdf")

plt.show()

