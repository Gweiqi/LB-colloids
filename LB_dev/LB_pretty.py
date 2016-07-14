import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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
    
    plt.imshow(img, cmap=mpl.cm.Accent_r, interpolation='nearest')
    plt.tick_params(axis='both', which='both',bottom='off',top='off',
                    labelbottom='off',right='off',left='off', labelleft='off')
    plt.imshow(uxy, cmap=mpl.cm.nipy_spectral, norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax))#cmap=mpl.cm.nipy_spectral,
    plt.colorbar()
    numit = str(numit)
    if len(numit) < 5:
        x = 5-len(numit)
        numit = '0'*x+numit
    Figname = name[:-5]+'_'+numit+'.png'
    plt.savefig(Figname)
    plt.close()
    
'''
u=np.arange(64)
u.shape = (2,8,4)

image = np.array([[True,False,True,True],[False,True,True,True],
               [True,False,True,True],[False,True,True,True],
               [True,False,True,True],[False,True,True,True],
               [True,False,True,True],[False,True,True,True]])
i = 100
velocity_image(u,image, 'testr.hdf5', i)
'''
