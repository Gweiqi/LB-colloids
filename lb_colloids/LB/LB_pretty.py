"""
LB_pretty contains methods to format and generate matplotlib images of
lattice Boltzmann velocity results. These methods are called internally
in LB2DModel if the user specifies plotting within the output dictionary.

example code to save an image of fluid magnitude within a domain from hdf5 output

>>> from lb_colloids.LB.LB_pretty import velocity_image
>>> from lb_colloids import ColloidOutput
>>>
>>> model = ColloidOutput.Hdf5Reader("LB_model.hdf5")
>>> u_y = model.get_data("lb_velocity_y")
>>> u_x = model.get_data("lb_velocity_x")
>>> image = model.get_data("image")
>>> u = np.array([u_y, u_x])
>>> velocity_image(u, image, "LB", 1000, vel=False, vmin=0.0001, vmax=0.1)
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys


def velocity_image(u, img, name, numit, vel, vmin, vmax):
    """
    Lattice Boltmann fluid velocity image generation routine. Can return
    plot velocity or fluid magnitude.

    Parameters:
    ----------
    :param np.ndarray u: Lattice Boltzmann macroscopic velocity array
    :param np.ndarray img: LBImage binarized image array
    :param str name: base name to save figure to
    :param int numit: time step that image is produced at
    :param bool vel: velocity flag, if false magnitude is plotted
    :param float vmin: Matplotlib vmin
    :param float vmax: Matplotlib vmax
    """
    temp = np.ones((len(img), len(img[0])))
    imgi = np.invert(img.astype(np.bool))
    if vel:
        uxy = u[0]
        uxy = np.array([uxy[i] * imgi[i] for i in range(len(uxy))])
        img = np.array([temp[i] * img[i] for i in range(len(uxy))])
        img = np.ma.masked_where(img == 0, img)
        uxy = np.ma.masked_where(uxy == 0, uxy)
    else:
        uxy = np.sqrt(u[0] * u[0] + u[1] * u[1])
        uxy = np.array([uxy[i] * imgi[i] for i in range(len(uxy))])
        img = np.array([temp[i] * img[i] for i in range(len(uxy))])
        img = np.ma.masked_where(img == 0, img)
        uxy = np.ma.masked_where(uxy == 0, uxy)
    
    plt.imshow(img, cmap=mpl.cm.Dark2_r, interpolation='nearest')
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    labelbottom='off', right='off', left='off', labelleft='off')

    cmap = set_colormap(vel)
    # if we plot y-velocity do not use log-normal distribution of colors
    if vel:
        plt.imshow(uxy, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        plt.imshow(uxy, cmap=mpl.cm.jet, norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
    plt.colorbar()
    numit = str(numit)
    if len(numit) < 5:
        x = 5-len(numit)
        numit = '0' * x+numit
    figname = name[:-5] + '_' + numit + '.png'
    plt.savefig(figname)
    plt.close()


def set_colormap(vel):
    """
    Method to set the matplotlib colormap option based
    on python version and datatype

    Parameters:
    ----------
    :param bool vel: Are we plotting a yvelocity image

    Returns:
    -------
    :return: (object) Matplotlib colormap object
    """
    whichpython = sys.version_info
    if vel:
        if whichpython[0] > 2 or whichpython[2] >= 10:
            cmap = mpl.cm.viridis_r
        else:
            cmap = mpl.cm.jet_r
    else:
        # set the color map based on python version
        if whichpython[0] > 2 or whichpython[2] >= 10:
            cmap = mpl.cm.viridis
        else:
            cmap = mpl.cm.jet
    return cmap
