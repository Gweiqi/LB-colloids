import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys

def velocity_image(u, img, name, numit, vel, vmin, vmax):    
    temp = np.ones((len(img),len(img[0])))
    imgi = np.invert(img.astype(np.bool))
    if vel == True:
	uxy = u[0]
        uxy = np.array([uxy[i]*imgi[i] for i in range(len(uxy))])
	img = np.array([temp[i]*img[i] for i in range(len(uxy))])
        img = np.ma.masked_where(img == 0, img)
	uxy = np.ma.masked_where(uxy == 0, uxy)
    else:
	uxy = np.sqrt(u[0]*u[0]+u[1]*u[1])
    	uxy = np.array([uxy[i]*imgi[i] for i in range(len(uxy))])
    	img = np.array([temp[i]*img[i] for i in range(len(uxy))])
    	img = np.ma.masked_where(img == 0, img)
    	uxy = np.ma.masked_where(uxy == 0, uxy)

    
    plt.imshow(img, cmap=mpl.cm.Dark2_r, interpolation='nearest')
    plt.tick_params(axis='both', which='both',bottom='off',top='off',
                    labelbottom='off',right='off',left='off', labelleft='off')

    
    cmap = set_colormap(vel)
    # if we plot y-velocity do not use log-normal distribution of colors
    if vel == True:	    
        plt.imshow(uxy, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        plt.imshow(uxy, cmap=mpl.cm.jet, norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
    plt.colorbar()
    numit = str(numit)
    if len(numit) < 5:
        x = 5-len(numit)
        numit = '0'*x+numit
    Figname = name[:-5]+'_'+numit+'.png'
    plt.savefig(Figname)
    plt.close()

def set_colormap(yvel):
    """
    Set colormap based on python version and datatype

    Input:
    -----
    yvel: (bool) Are we plotting a yvelocity image

    Returns:
    --------
    cmap: (object) Matplotlib colormap object
    """
    whichpython = sys.version_info
    if yvel == True:
        if whichpython[0] > 2 or whichpython[2] >=10:
	    cmap = mpl.cm.viridis_r
	else:
	    cmap = mpl.cm.jet_r
    else:
        # set the color map based on python version
        if whichpython[0] > 2 or whichpython[2] >=10:
            cmap = mpl.cm.viridis
        else:
            cmap = mpl.cm.jet
    return cmap
    
