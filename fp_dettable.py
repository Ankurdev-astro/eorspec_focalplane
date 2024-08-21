"""
This script handles the creation of detector Qtables for EoR-Spec in TOAST3 format.
It includes functions to simulate wafers and manage detector data for the focal plane.
"""
## File handles making detector Qtables for EoR-Spec in TOAST3 format

#imports
import numpy as np
import toast.qarray as qa
from toast.instrument_sim import rhomb_dim, rhombus_hex_layout
from toast.instrument_coords import xieta_to_quat
import astropy.units as u
from astropy.table import Column, QTable, vstack
from copy import deepcopy

#Functions
def sim_wafer(wafertype: str, wafername: str) -> dict: 
    """
    Function to simulate a wafer for EoR-Spec focalplane
    Args:
    wafertype (str): 'lfa' or 'hfa'
    wafername (str): name of the wafer
    
    Returns:
    dict: dictionary containing the wafer information
    """ 
    
    # The plate scale in degrees / mm.
    platescale = 0.00495 #SO line #917
        
    if wafertype == 'lfa':
        pix_mm = 2.75 #CCAT Parameter Specs Doc
        npix=1728
    elif wafertype == 'hfa':
        pix_mm = 2.09 #CCAT Parameter Specs Doc
        npix=3072
    else: raise ValueError(f'Wafertype {wafertype} not supported. Must be lfa or hfa')
    
    # The number of pixel locations in one rhombus.
    nrhombus = npix // 3
    # This dim is the number of pixels along the short axis.
    dim = rhomb_dim(nrhombus)
    
    # Pixel separation in degrees
    # https://github.com/simonsobs/sotodlib/blob/master/sotodlib/toast/sim_focalplane.py#L82
    pixsep = platescale * pix_mm * u.degree #SO line 82
    
    # https://github.com/simonsobs/sotodlib/blob/master/sotodlib/toast/sim_focalplane.py#L97
    # This is the center-center distance along the short axis
    # The angle subtended by pixel center-to-center distance 
    # along the short dimension of one rhombus.
    width = (dim - 1) * pixsep #SO Line 97
    
    # Setting pol to zero
    # pol = u.Quantity(np.zeros(nrhombus, dtype=np.float64), u.degree)
    pol = None
    
    # By default, the rhombi are aligned such that the center spacing across
    # the gap is the same as the center spacing within a rhombus.
    # https://github.com/simonsobs/sotodlib/blob/master/sotodlib/toast/sim_focalplane.py#L87
    gap = 0.0 * u.degree #additional gap between the edges of the rhombi.
    suffix = "_" + wafername
    wafer = rhombus_hex_layout(
            nrhombus, width, "eorspec_", suffix, gap=gap, pol=pol
        )
    
    return wafer

def make_det_table(focalplane_dict):
    """
    Function to make a full EoR-Spec detector table
    Three EoR-Spec wafers are simulated and arranged in a tube
        
    Args:
    focalplane_dict (dict): dictionary containing the focalplane information
    
    Returns:
    astropy.table.QTable: detector table
    """
    
    #Setting arbitrary values for now
    bandcenter = 280 * u.GHz
    bandwidth = 2 * u.GHz
    sample_rate = 244 * u.Hz
    psd_net = 1 * u.K * np.sqrt(1 * u.second)
    psd_fmin = 0.1 * u.Hz
    psd_alpha = 1
    psd_fknee = 1 * u.Hz
    # width=1 * u.degree
    fwhm_lfa = 0.8 * u.arcmin #280 GHz
    fwhm_hfa = 0.62 * u.arcmin #350 HHz    

    det_data = {x: focalplane_dict[x] for x in sorted(focalplane_dict.keys())}
    n_det = len(det_data)
    
    nominal_freq = str(int(bandcenter.to_value(u.GHz)))
    det_names = [f"{x}-{nominal_freq}" for x in det_data.keys()]
    det_gamma = u.Quantity([det_data[x]["gamma"] for x in det_data.keys()], u.radian)
    
    det_table = QTable(
    [
        Column(name="name", data=det_names),
        Column(name="wname", length=n_det, dtype='S4'),
        Column(name="wtype", length=n_det, dtype='S4'),
        Column(name="quat", data=[det_data[x]["quat"] for x in det_data.keys()]),
        Column(name="psi_pol", length=n_det, unit=u.rad),
        Column(name="gamma", length=n_det, unit=u.rad),
        Column(name="fwhm", length=n_det, unit=u.arcmin),
        Column(name="psd_fmin", length=n_det, unit=u.Hz),
        Column(name="psd_fknee", length=n_det, unit=u.Hz),
        Column(name="psd_alpha", length=n_det, unit=None),
        Column(name="psd_net", length=n_det, unit=(u.K * np.sqrt(1.0 * u.second))),
        Column(name="bandcenter", length=n_det, unit=u.GHz),
        Column(name="bandwidth", length=n_det, unit=u.GHz),
    ]
    )
    # det_table.add_column(Column(name="wtype", length=n_det, dtype='S4'))
        
    for idet, det in enumerate(det_data.keys()):
        det_table[idet]["wname"] = det.split("_")[-1]
        if "lfa" in det_table[idet]["wname"]:
            det_table[idet]["wtype"] = str("lfa")
        elif "hfa" in det_table[idet]["wname"]:
            det_table[idet]["wtype"] = str("hfa")
        else: raise ValueError('Keys must include lfa or hfa')            
        
        # psi_pol is the rotation from the PXX beam frame to the polarization
        # sensitive direction.
        det_table[idet]["psi_pol"] = 0 * u.rad
        det_table[idet]["gamma"] = det_gamma[idet]
        if "lfa" in str(det):
            det_table[idet]["fwhm"] = fwhm_lfa
        elif "hfa" in str(det):
            det_table[idet]["fwhm"] = fwhm_hfa
        else: raise ValueError('Keys must include lfa or hfa')
        det_table[idet]["bandcenter"] = bandcenter
        det_table[idet]["bandwidth"] = bandwidth 
        det_table[idet]["psd_fmin"] = psd_fmin
        det_table[idet]["psd_fknee"] = psd_fknee
        det_table[idet]["psd_alpha"] = psd_alpha
        det_table[idet]["psd_net"] = psd_net
    
    return det_table


##Creating EoR-Spec detectors and tube
#Focal plane wafer dicts for eorspec wafers
eorspec_lfa1 = sim_wafer(wafertype="lfa", wafername="lfa1")
eorspec_lfa2 = sim_wafer(wafertype="lfa", wafername="lfa2")
eorspec_hfa = sim_wafer(wafertype="hfa", wafername="hfa")

#Make det tables from eorspec wafers
det_table_lfa1 = make_det_table(eorspec_lfa1)
det_table_lfa2 = make_det_table(eorspec_lfa2)
det_table_hfa = make_det_table(eorspec_hfa)

#Defining quaternion wafer centers
#https://github.com/simonsobs/sotodlib/blob/master/sotodlib/sim_hardware.py#L478
## This tube spacing in mm corresponds to 1.78 degrees projected on
## the sky at a plate scale of 0.00495 deg/mm.
platescale = 0.00495 #SO line #478
# https://github.com/simonsobs/sotodlib/blob/master/sotodlib/sim_hardware.py#L385
waferspace = 128.4 #SO line #385

#https://github.com/simonsobs/sotodlib/blob/master/sotodlib/toast/sim_focalplane.py#L690   
wradius = 0.5 * (waferspace * platescale * np.pi / 180.0) #SO Line 690

#https://github.com/simonsobs/sotodlib/blob/master/sotodlib/toast/sim_focalplane.py#L696
qwcenters = [
                xieta_to_quat(
                    -wradius, wradius / np.sqrt(3.0), 0.0
                ),
                xieta_to_quat(
                    wradius, wradius / np.sqrt(3.0), 0.0
                ),
                xieta_to_quat(
                    0.0, -2.0 * wradius / np.sqrt(3.0), 0.0
                ),
            ]

thirty = np.pi / 6.0  # 30 degrees in radians
sixty = np.pi / 3.0
zaxis = np.array([0, 0, 1], dtype=np.float64)

#Rotating wafers
cp_dettable_lfa1 = deepcopy(det_table_lfa1)
cp_dettable_lfa2 = deepcopy(det_table_lfa2)
cp_dettable_hfa = deepcopy(det_table_hfa)
rotation_quat = qa.from_axisangle(zaxis, sixty)

#HFA is at position 1
#LFA takes position 0,2
# This is arbitrary. Self definition

for det,det_quat in enumerate(cp_dettable_lfa1["quat"]):
    rotated_quat = qa.mult(rotation_quat, det_quat)
    quat_shift = qa.mult(qwcenters[0], rotated_quat)
    cp_dettable_lfa1["quat"][det] = quat_shift

for det,det_quat in enumerate(cp_dettable_hfa["quat"]):
    rotated_quat = qa.mult(rotation_quat, det_quat)
    quat_shift = qa.mult(qwcenters[1], rotated_quat)
    cp_dettable_hfa["quat"][det] = quat_shift

for det,det_quat in enumerate(cp_dettable_lfa2["quat"]):
    rotated_quat = qa.mult(rotation_quat, det_quat)
    quat_shift = qa.mult(qwcenters[2], rotated_quat)
    cp_dettable_lfa2["quat"][det] = quat_shift 

dettable_stack = vstack([
                    cp_dettable_lfa1,
                    cp_dettable_lfa2,
                    cp_dettable_hfa
                    ])


#No plotting here
# plot_focalplane_helper(
#     focalplane=fp_test_edit8,
#     width=width,
#     height=width,
#     show_labels=False,
#     outfile="fp_eorspec_rot2_2024.pdf")

# plt.show()


# save dettable_stack to hdf5 file using astropy write
print(" Writing detector table to HDF5 file ...")
dettable_stack.write('./eorspec_dettable.h5', path='dettable_stack', 
                                serialize_meta=True,overwrite=True)


