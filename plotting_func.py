"""
This module contains helper functions for plotting the EoR-Spec full focal plane and annuli.
It includes functions to visualize the arrangement of detectors and their properties.
"""
### Helper EoR-Spec FP plotting functions
### Modified from: 
# https://github.com/hpc4cmb/toast/blob/cc7a91abb8ddf86a31272c8c49de633dcf4d2d06/src/toast/instrument_sim.py#L1048

#Imports
import numpy as np
import astropy.units as u
from toast import qarray as qa
from toast.instrument_coords import quat_to_xieta
from toast.vis import set_matplotlib_backend
from matplotlib.colors import BoundaryNorm

def plot_focalplane_eorspec(
    focalplane=None,
    width=None,
    height=None,
    outfile=None,
    show_labels=False,
    face_color=None,
    pol_color=None,
    xieta=False,
    show_centers=False,
    show_gamma=False,
):
    """Visualize a projected EoRSpec Focalplane.

    This makes a simple plot of the detector positions on the projected focalplane.
    By default, this plots the focalplane in the boresight X / Y / Z frame, as seen
    by incoming photons.  If `xieta` is set to True, the focalplane is plotted in
    Xi / Eta / Gamma coordinates as seen from the observer looking out at the sky.

    To avoid python overhead in large MPI jobs, we place the matplotlib import inside
    this function, so that it is only imported when the function is actually called.

    Args:
        focalplane (Focalplane):  The focalplane to plot
        width (Quantity):  Width of plot.
        height (Quantity):  Height of plot.
        outfile (str):  Output PDF path.  If None, then matplotlib will be
            used for inline plotting.
        show_labels (bool):  If True, plot detector names.
        face_color (dict): dictionary of color values for the face of each
            detector circle.
        pol_color (dict): dictionary of color values for the polarization
            arrows.
        xieta (bool):  Plot in observer xi/eta/gamma coordinates rather than
            boresight X/Y/Z.
        show_centers (bool):  If True, label the pixel centers.
        show_gamma (bool):  If True, show gamma angle (for debugging).

    Returns:
        (Figure):  The figure.

    """
    if focalplane is None:
        raise RuntimeError("You must specify a Focalplane instance")

    if outfile is not None:
        set_matplotlib_backend(backend="pdf")

    import matplotlib.pyplot as plt

    if width is None:
        width = 10.0 * u.degree

    if height is None:
        height = 10.0 * u.degree

    width_deg = width.to_value(u.degree)
    height_deg = height.to_value(u.degree)

    # xfigsize = int(width_deg) + 1
    # yfigsize = int(height_deg) + 1
    xfigsize = 6
    yfigsize = 6
    figdpi = 100

    # Compute the font size to use for detector labels
    fontpix = 0.05 * figdpi
    fontpt = int(0.75 * fontpix)

    fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
    # fig = plt.figure(figsize=(5, 5), dpi=figdpi)
    ax = fig.add_subplot(1, 1, 1)

    half_width = 0.6 * width_deg
    half_height = 0.6 * height_deg
    if xieta:
        ax.set_xlabel(r"Boresight $\xi$ Degrees", fontsize="medium")
        ax.set_ylabel(r"Boresight $\eta$ Degrees", fontsize="medium")
    else:
        ax.set_xlabel("Boresight X Degrees", fontsize="medium")
        ax.set_ylabel("Boresight Y Degrees", fontsize="medium")
    ax.set_xlim([-half_width, half_width])
    ax.set_ylim([-half_height, half_height])

    xaxis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    yaxis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    zaxis = np.array([0.0, 0.0, 1.0], dtype=np.float64)

    for d in focalplane.detectors:
        quat = focalplane[d]["quat"]
        fwhm = focalplane[d]["fwhm"].to_value(u.arcmin)

        # radius in degrees
        detradius = 0.5 * 5.0 / 60.0
        if fwhm is not None:
            detradius = 0.5 * fwhm / 60.0
            # detradius = 0.4 * fwhm / 60.0
            
        if xieta:
            xi, eta, gamma = quat_to_xieta(quat)
            xpos = xi * 180.0 / np.pi
            ypos = eta * 180.0 / np.pi
            # Polang is plotted relative to visualization x/y coords
            polang = 1.5 * np.pi - gamma
            plot_gamma = polang
        else:
            # rotation from boresight
            rdir = qa.rotate(quat, zaxis).flatten()
            mag = np.arccos(rdir[2]) * 180.0 / np.pi
            ang = np.arctan2(rdir[1], rdir[0])
            orient = qa.rotate(quat, xaxis).flatten()
            polang = np.arctan2(orient[1], orient[0])
            xpos = mag * np.cos(ang)
            ypos = mag * np.sin(ang)
            xi, eta, gamma = quat_to_xieta(quat)
            plot_gamma = gamma

        detface = "gray"
        if face_color is not None:
            detface = face_color[d]

        # circ = plt.Circle((xpos, ypos), radius=detradius, fc=detface, ec="gray")
        circ = plt.Circle((xpos, ypos), radius=detradius, fc=detface)
        ax.add_artist(circ)

        ascale = 1.5

        xtail = xpos - ascale * detradius * np.cos(polang)
        ytail = ypos - ascale * detradius * np.sin(polang)
        dx = ascale * 2.0 * detradius * np.cos(polang)
        dy = ascale * 2.0 * detradius * np.sin(polang)

        detcolor = "black"
        if pol_color is not None:
            detcolor = pol_color[d]

        if show_centers:
            ysgn = -1.0
            if dx < 0.0:
                ysgn = 1.0
            ax.text(
                (xpos + 0.1 * dx),
                (ypos + 0.1 * ysgn * dy),
                f"({xpos:0.4f}, {ypos:0.4f})",
                color="green",
                fontsize=fontpt,
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(fc="w", ec="none", pad=1, alpha=0.0),
            )

        if show_labels:
            xsgn = 1.0
            if dx < 0.0:
                xsgn = -1.0
            labeloff = 0.05 * xsgn * fontpix * len(d) / figdpi
            ax.text(
                (xtail + 1.3 * dx + labeloff),
                (ytail + 1.2 * dy),
                d,
                color="k",
                fontsize=fontpt,
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(fc="w", ec="none", pad=1, alpha=0.0),
            )

        if show_gamma:
            ax.arrow(
                xtail,
                ytail,
                1.3 * dx,
                1.3 * dy,
                width=0.1 * detradius,
                head_width=0.2 * detradius,
                head_length=0.2 * detradius,
                fc="gray",
                ec="gray",
                length_includes_head=True,
            )


    # Draw a "mini" coordinate axes for reference
    # xmini = -0.8 * half_width
    # ymini = -0.8 * half_height
    # xlen = 0.1 * half_width
    # ylen = 0.1 * half_height
    # mini_width = 0.005 * half_width
    # mini_head_width = 3 * mini_width
    # mini_head_len = 3 * mini_width
    # if xieta:
    #     aprops = [
    #         (xlen, 0, "-", r"$\xi$"),
    #         (0, ylen, "-", r"$\eta$"),
    #         (-xlen, 0, "--", "Y"),
    #         (0, -ylen, "--", "X"),
    #     ]
    # else:
    #     aprops = [
    #         (xlen, 0, "-", "X"),
    #         (0, ylen, "-", "Y"),
    #         (-xlen, 0, "--", r"$\eta$"),
    #         (0, -ylen, "--", r"$\xi$"),
    #     ]
    # for ap in aprops:
    #     lx = xmini + 1.5 * ap[0]
    #     ly = ymini + 1.5 * ap[1]
    #     lw = figdpi / 200.0
    #     ax.arrow(
    #         xmini,
    #         ymini,
    #         ap[0],
    #         ap[1],
    #         width=mini_width,
    #         head_width=mini_head_width,
    #         head_length=mini_head_len,
    #         fc="k",
    #         ec="k",
    #         linestyle=ap[2],
    #         linewidth=lw,
    #         length_includes_head=True,
    #     )
    #     ax.text(
    #         lx,
    #         ly,
    #         ap[3],
    #         color="k",
    #         fontsize=int(figdpi / 10),
    #         horizontalalignment="center",
    #         verticalalignment="center",
    #     )

    # st = "Focalplane Looking Towards Observer"
    # if xieta:
    #     st = "Focalplane on Sky From Observer"
    # fig.suptitle(st)

    if outfile is None:
        output_plt = plt.show();
        print()
    else:
        plt.savefig(outfile, dpi=figdpi, bbox_inches="tight", format="pdf")
        plt.close()
    return fig


def plot_eorspec_annuli(
    focalplane=None,
    outfile=None,
    label_step=False,
):
    """Visualize different annuli of EoR-Spec

    This makes a simple plot of the detector positions on the projected
    focalplane.  The size of detector circles are controlled by the detector
    "fwhm" key, which is in arcminutes.
    
    Arguments:
    focalplane (Focalplane): The focalplane to plot
    outfile (str): Output PDF path. If None, then matplotlib will be used for inline plotting.
    show_labels (bool): If True, plot EoR-Spec wafers labels.
    label_step (bool): If True, includes FPI step info in label
    
    Returns:
    None

    """
    
    if focalplane is None:
        raise RuntimeError("You must specify a Focalplane instance")

    if outfile is not None:
        set_matplotlib_backend(backend="pdf")

    import matplotlib.pyplot as plt

    width = 1.3 * u.degree 
    height = 1.3 * u.degree

    width_deg = width.to_value(u.degree)
    height_deg = height.to_value(u.degree)

    xfigsize = 8
    yfigsize = 8
    figdpi = 100

    # Compute the font size to use for detector labels
    # fontpix = 0.05 * figdpi
    # fontpt = int(0.75 * fontpix)

    fig = plt.figure(figsize=(xfigsize, yfigsize), dpi=figdpi)
    ax = fig.add_subplot(1, 1, 1)
    fig.patch.set_facecolor('whitesmoke')

    half_width = 0.6 * width_deg
    half_height = 0.6 * height_deg

    xaxis = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    yaxis = np.array([0.0, 1.0, 0.0], dtype=np.float64)
    zaxis = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    
    # Initialize the colormap and the Normalize object
    cmap = plt.cm.get_cmap('nipy_spectral_r')
    # Define the boundaries and normalize
    boundaries = np.arange(210, 425, 2)
    norm = BoundaryNorm(boundaries, cmap.N)
    
    if label_step:
        for d in focalplane.detectors:
            try:
                fpi_step = focalplane[d]["fpi_step"]
                ax.text(0.4,0.6, 
                fr"FPI step $f_0$: {fpi_step.split('step')[1]} GHz", 
                color='black', fontsize="large", 
                horizontalalignment='center',
                verticalalignment='center'
                )
                
                print(f"Plotting FPI Step: {fpi_step}")
                break
            except:
                pass
            
        #{fpi_steps[step_FPI]}  
  

    #Loop over all detectors
    for d in focalplane.detectors:
        quat = focalplane[d]["quat"]
        fwhm = focalplane[d]["fwhm"].to_value(u.arcmin)
        wtype  = focalplane[d]["wtype"]
        annuli_name = focalplane[d]["annuli_name"]
        freq_channel = focalplane[d]["freq_channel"]

        # radius in degrees
        detradius = 0.5 * fwhm / 60.0

        # rotation from boresight
        rdir = qa.rotate(quat, zaxis).flatten()
        mag = np.arccos(rdir[2]) * 180.0 / np.pi
        ang = np.arctan2(rdir[1], rdir[0])
        orient = qa.rotate(quat, xaxis).flatten()
        polang = np.arctan2(orient[1], orient[0])
        xpos = mag * np.cos(ang)
        ypos = mag * np.sin(ang)
        xi, eta, gamma = quat_to_xieta(quat)

        #EoR-Spec Det Face Color
        
        
        if len(annuli_name) != 0:
            detface = cmap(norm(freq_channel))
            # print(detface)
        else:
            if wtype == "lfa":
                detface = "dimgrey"
            elif wtype == "hfa":
                detface = "slategrey"
            else:
                detface = "black"
        

        # circ = plt.Circle((xpos, ypos), radius=detradius, fc=detface, ec="gray")
        circ = plt.Circle((xpos, ypos), radius=detradius, fc=detface)
        ax.add_artist(circ)

        ascale = 1.5

    #Plotting frame details
    plt.title("EoR-Spec projected Focal Plane Arrays", fontsize="x-large")
    ax.set_xlabel("Boresight X Degrees", fontsize="medium")
    ax.set_ylabel("Boresight Y Degrees", fontsize="medium")
    ax.set_xlim([-half_width, half_width])
    ax.set_ylim([-half_height, half_height])
    
    # Create a new axis for the colorbar
    #[left, bottom, width, height] weights[0:1]
    cb_ax = fig.add_axes([0.91, 0.124, 0.03, 0.754])
    
    # Mapping the norm values to the colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Placeholder

    reduced_ticks = boundaries[::10]
    
    # Add the colorbar
    cbar = plt.colorbar(sm, boundaries=boundaries, 
                        ticks=reduced_ticks, 
                        orientation='vertical', cax=cb_ax)
    cbar.set_label('EoR-Spec Band Frequencies [GHz]', size="x-large")

    if outfile is None:
        output_plt = plt.show();
        print()
    elif outfile.endswith(".png"):
        plt.savefig(outfile, dpi=figdpi, bbox_inches="tight", format="png")
        plt.close()
    else:
        plt.savefig(outfile, dpi=figdpi, bbox_inches="tight", format="pdf")
        plt.close()
    return fig
    
