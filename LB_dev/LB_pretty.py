import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys

def velocity_image(u, img, name, numit, vel, vmin, vmax):    
    temp = np.ones((len(img),len(img[0])))
    imgi = np.invert(img)
    if vel == 'uy':
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

    whichpython = sys.version_info
    
    plt.imshow(img, cmap=mpl.cm.Dark2_r, interpolation='nearest')
    plt.tick_params(axis='both', which='both',bottom='off',top='off',
                    labelbottom='off',right='off',left='off', labelleft='off')
    if whichpython[0] > 2 or whichpython[2] >=10:
        plt.imshow(uxy, cmap=mpl.cm.viridis, norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax))
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
    
