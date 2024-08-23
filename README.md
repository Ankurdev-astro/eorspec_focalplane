# eorspec_focalplane

This set of scripts makes detector tables for EoR-Spec focalplane simulation in the TOASTv3 format.

The two cases are:
* We want to view the full focalplane with different frequencies at a single FPI step
* We want to simulate 1 frequency channel that are illuminated at multiple FPI steps

Final output(s) is/are Astropy Detector table(s) with Detector parameters which can provided as argument
for a TOAST FocalPlane Class.


![Example EoR-Spec Focal Plane simulation for an FPI step](fpi_data/fpi_plots/fp_eorspec_step.png)


EoR-Spec Frequency and Annulus data taken from:
https://github.com/ccatobs/eor_spec_mapping_simulations/blob/main/params/annulus_radii.csv





