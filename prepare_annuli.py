"""
This script contains functions to prepare and mask detectors for annuli in the EoR-Spec focal plane.
It includes functionality to filter detectors based on radius and other parameters for TOAST simulations.
"""
#Imports
from astropy.table import Column
import toast.qarray as qa
import numpy as np

def mask_dets_annulus(
        dets_table, r_min_rad, r_max_rad, wtype, 
        freq_channel, fpi_step, annulus_name,
        remove_dets=False
        ):
    
    """
    Returns a detector table with dets to be considered for the toast simulation of a particular
    annulus.
    
    Arguments:
    dets_table = astropy table with detector information
    r_min_rad = minimum radius value in radians
    r_max_rad = maximum radius value in radians
    wtype = "lfa" or "hfa"
    freq_channel = frequency channel of EoR-Spec
    fpi_step = FPI step e.g. "step210"
    annulus_name = name of the annuli e.g. "A1"
    remove_dets = remove detectors not in the annulus if True   
    
    Returns:
    dets_table = astropy table with selected detectors
    """    
    zaxis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    dets_selected = 0
    dets_idx_delete = []
    wtype = wtype.lower()
    
    n_det = len(dets_table)
    # dets_table.add_column(Column(name="freq_channel", length=n_det, dtype='i4'))
    # dets_table.add_column(Column(name="fpi_step", length=n_det, dtype='S8'))
    # dets_table.add_column(Column(name="annuli_name", length=n_det, dtype='S8'))
    
    # Add columns only if they don't already exist
    if "freq_channel" not in dets_table.colnames:
        dets_table.add_column(Column(name="freq_channel", length=n_det, dtype='i4', data=[0]*n_det))
    if "fpi_step" not in dets_table.colnames:
        dets_table.add_column(Column(name="fpi_step", length=n_det, dtype='S8', data=['']*n_det))
    if "annuli_name" not in dets_table.colnames:
        dets_table.add_column(Column(name="annuli_name", length=n_det, dtype='S8', data=['']*n_det))
    
    for idet, quat in enumerate(dets_table["quat"]):
        det_wtype = dets_table[idet]["wtype"]
        
        # rotation from boresight
        rdir = qa.rotate(quat, zaxis).flatten()
        ang = np.arctan2(rdir[1], rdir[0])

        # get x, y positions, in radians
        # angular distance
        mag = np.arccos(rdir[2])
        xpos = mag * np.cos(ang)
        ypos = mag * np.sin(ang)

        r = np.sqrt(xpos**2 + ypos**2)  # distance from the origin
        if r > r_min_rad and r < r_max_rad and det_wtype == wtype:
            dets_table[idet]["freq_channel"] = freq_channel
            dets_table[idet]["fpi_step"] = fpi_step
            dets_table[idet]["annuli_name"] = annulus_name
            dets_selected += 1
        else:
            dets_idx_delete.append(idet)

    if remove_dets:
        dets_table.remove_rows(dets_idx_delete)
    print(f"Number of Detectors selcted: {dets_selected}")

    return dets_table


def process_annuli_fpistep(dets_table, fpistep_infotxt, fpi_step):
    """
    Returns a detector table with dets to be considered for the toast simulation of a particular
    FPI step
    
    Arguments:
    dets_table = astropy table with detector information
    fpistep_infotxt = text file with annulus information for each FPI step
    fpi_step = FPI step
    
    Returns:
    dets_table = astropy table with updated detectors    
    """
    
    # Read text file and skip the header
    with open(fpistep_infotxt, 'r') as f:
        lines = f.readlines()[1:]
    
    for line in lines:
        col_data = line.strip().split('\t')
       
        wtype = col_data[0]
        annulus_num = int(col_data[1])
        freq_channel = int(col_data[2])
        r_min = float(col_data[3])
        r_max = float(col_data[4])
        # freq_min = float(col_data[5])
        # freq_max = float(col_data[6])
        freq_center = float(col_data[7])
        freq_delta = float(col_data[8])
    
        # Create annulus_name
        annulus_name = f"A{annulus_num}"
        print(wtype, annulus_name, freq_channel)

        dets_table = mask_dets_annulus(dets_table, r_min, r_max,
                                        wtype, freq_channel, fpi_step, 
                                        annulus_name)
        
    return dets_table

def process_annuli_fchl(dets_table, fpistep_infotxt, target_fchl, fpi_step):
    """
    Returns a detector table with dets to be considered for the toast simulation of a particular
    FPI step and frequency channel
    
    Arguments:
    dets_table = astropy table with detector information
    fpistep_infotxt = text file with annulus information for each FPI step
    target_fchl = target frequency channel
    fpi_step = FPI step
    
    Returns:
    dets_table = astropy table with updated detectors 
    """
    
    # Read text file and skip the header
    with open(fpistep_infotxt, 'r') as f:
        lines = f.readlines()[1:]
    
    # Flag to track if the target frequency channel is found
    found_fchl = False
    
    for line in lines:
        col_data = line.strip().split('\t')
       
        wtype = col_data[0]
        annulus_num = int(col_data[1])
        freq_channel = int(col_data[2])
        r_min = float(col_data[3])
        r_max = float(col_data[4])
        # freq_min = float(col_data[5])
        # freq_max = float(col_data[6])
        freq_center = float(col_data[7])
        freq_delta = float(col_data[8])
    
        if freq_channel != target_fchl:
            continue
        elif freq_channel == target_fchl:
            found_fchl = True
            # Create annulus_name
            annulus_name = f"A{annulus_num}"
            print(wtype, annulus_name, freq_channel)
            # print("Debug point 1")

            dets_table = mask_dets_annulus(dets_table, r_min, r_max,
                                            wtype, freq_channel, fpi_step, 
                                            annulus_name, remove_dets=True)
            break
         
    if not found_fchl:
        # print(f"Frequency channel {target_fchl} GHz not found in {fpi_step} !")
        return None
    return dets_table