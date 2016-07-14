import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate
import h5py as H
import optparse
import sys

class HDF5_reader:
    def __init__(self, HDF_name):
        hdf = H.File(HDF_name, 'r+')
        self.imarray = hdf['Binary_image'][()]
        self.uarray = hdf['results/uarray'][()]
        self.yu = hdf['results/uarray'][()][0]
        self.xu = hdf['results/uarray'][()][1]

class gridarray:
    def __init__(self, arr, gridres, gridsplit):
        self.yarr = np.copy(arr).T
        self.gridx = self.arrx(arr, gridres, gridsplit)
        self.gridy = self.arry(self.yarr, gridres, gridsplit).T

    def arrx(self, arr, gridres, gridsplit):
        for line in range(len(arr)):
            arr[line] = self.distance_gridx(arr[line], gridres, gridsplit)
        return arr

    def arry(self, yarr, gridres, gridsplit):
        for line in range(len(yarr)):
            ylen = len(yarr[0])
            yarr[line] = self.distance_gridy(yarr[line], gridres, ylen, gridsplit)
            #print yarr[line]
        return yarr
    
    def distance_gridx(self,line,gridres, gridsplit):
        rplc = np.linspace(0,1,gridsplit)
        rplc = rplc[gridsplit/2:]
        for i in rplc:
            line[i] = 1.    # pore space == 0, solid is still 1.
        line[line < 1] = 0. # linear interpolation also interpolates pore
        line[line > 1] = 0. # boundaries, this block corrects that error

        # create an array of pore boundaries from binary system
        boundary= np.where(np.abs(np.diff(line)) >= 1)[0]
        try:
            for i in range(1,len(boundary),2):
                rbound = boundary[i]+1
                lbound = boundary[i-1]+1
                gap = rbound - lbound

                if gap % 2 == 0:
                    gap = gap//2
                    left = np.arange(1,gap+1)*-1*gridres
                    right = left[::-1]*-1
                    line[lbound:rbound] = np.append(left,right)

                else:
                    gap = gap//2
                    left = np.arange(1,gap+1)*-1*gridres
                    right = left[::-1]*-1
                    adjust = (len(left)+1)*-1*gridres
                    left = np.append(left,adjust)
                    line[lbound:rbound] = np.append(left,right)
                    
        except IndexError:
            print 'volume does not percolate'
        return line

    def distance_gridy(self,line,gridres,ylen,gridsplit):
        rplc = np.linspace(0,1,gridsplit)
        rplc = rplc[gridsplit/2:]
        for i in rplc:
            line[i] = 1.    # pore space == 0, solid is still 1.
        line[line > 1.] = 0. # linear interpolation also interpolates pore
        line[line < 1.] = 0. # boundaries, this block corrects that error

        # create an array of pore boundaries from binary system
        boundary= np.where(np.abs(np.diff(line)) >= 1)[0]

        if len(boundary) > 0:
            rbound = boundary[0]+1
            lbound = 0
            top = np.arange(rbound,lbound,-1)*gridres*-1
            line[lbound:rbound] = top
               
            for i in range(1,len(boundary),2):
                rbound = boundary[i]+1
                lbound = boundary[i-1]+1
                gap = rbound - lbound
                if gap % 2 == 0:
                    gap = gap//2
                    left = np.arange(1,gap+1)*gridres
                    right = left[::-1]*-1
                    line[lbound:rbound] = np.append(left,right)
                else:
                    gap = gap//2
                    left = np.arange(1,gap+1)*gridres
                    right = left[::-1]*-1
                    adjust = (len(left)+1)*gridres
                    left = np.append(left,adjust)
                    line[lbound:rbound] = np.append(left,right)
            rbound = ylen
            lbound = boundary[-1]+1
            gap = rbound - lbound
            bottom =np.arange(1,gap+1)*gridres
            line[lbound:rbound] = bottom
        else:
            pass
        return line

def LB_varray(LBv, img):
    vel = np.zeros((len(LB.yu),len(LB.yu[0])))
    invert = np.invert(img)
    vel = np.array([LBv[i]*invert[i] for i in range(len(LBv))])
    print vel.shape
    return vel

def interp_v(LBv, gridsplit):
    ylen = len(LBv)
    xlen = len(LBv[0])
    xindex = np.arange(0,xlen)
    yindex = np.arange(0,ylen)
    xnew = np.arange(0, xlen-1+(1./gridsplit), 1./gridsplit)
    ynew = np.arange(0, ylen-1+(1./gridsplit), 1./gridsplit)
    f = interpolate.interp2d(xindex, yindex, LBv, kind='linear')
    znew = f(xnew, ynew)
    return znew

LB = HDF5_reader('Synthetic255.hdf5')
gridsplit = 10
gridres = 110/gridsplit

LBy = LB_varray(LB.yu, LB.imarray)
LBx = LB_varray(LB.xu, LB.imarray)
LBv = np.sqrt(LB.yu*LB.yu+LB.xu*LB.xu)
LBvp = LB_varray(LBv, LB.imarray)

plt.imshow(LBvp, interpolation='nearest')
plt.colorbar()
plt.show()

#### interpolate over grid array and interpolate veloctity profiles ####
Col_vp = interp_v(LBvp, gridsplit)
Col_img = interp_v(LB.imarray, gridsplit)

### Use grids function to measure distance from pore space and correct
### boundaries for interpolation effects. 
Grids = gridarray(Col_img, gridres, gridsplit)
xArr = Grids.gridx
yArr = Grids.gridy


xArr = np.ma.masked_where(yArr == 1., xArr)
plt.imshow(xArr, interpolation = 'nearest')
plt.show()

yArr = np.ma.masked_where(yArr == 1., yArr)
plt.imshow(yArr, interpolation = 'nearest')
plt.show()


plt.imshow(Col_img, cmap=mpl.cm.Accent, interpolation='nearest')
Col_vp = np.ma.masked_where(Col_vp == 0, Col_vp)
plt.imshow(Col_vp, interpolation = 'nearest')
plt.colorbar()
plt.show()

