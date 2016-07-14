import numpy as np
import matplotlib.pyplot as plt
import colloid_math as cm

check = np.zeros((18,18))
check[::4, 1::4] = True
check[1::2, ::4] = True
check[0,:] = 0.
check[-1,:] = 0.
check[:,0] = True
check[:,-1] = True

checky = np.copy(check.T)

gridres = 1e-6

################################################################
### Definitions to set grid distances based on vector forces ###
### Pos ++ Up and Right, Neg. -- Down and Left               ###
################################################################

def distance_gridx(line, gridres):
    boundary= np.where(np.abs(np.diff(line)) >= 1)[0]
    try:
        for i in range(1,len(boundary),2):
            rbound = boundary[i] + 1
            lbound = boundary[i-1] + 1
            gap = rbound - lbound

            if gap % 2 == 0:
                gap = gap//2
                left = np.arange(1, gap + 1)*-1*gridres
                right = left[::-1]*-1
                line[lbound:rbound] = np.append(left, right)

            else:
                gap = gap//2
                left = np.arange(1, gap + 1)*-1*gridres
                right = left[::-1]*-1
                adjust = (len(left) + 1)*-1*gridres
                left = np.append(left, adjust)
                line[lbound:rbound] = np.append(left, right)
                    
    except IndexError:
        print 'volume does not percolate'
    return line


def distance_gridy(line, gridres, ylen):
    boundary= np.where(np.abs(np.diff(line)) >= 1)[0]

    if len(boundary) > 0:
        rbound = boundary[0] + 1
        lbound = 0
        top = np.arange(rbound, lbound, -1)*gridres*-1
        line[lbound:rbound] = top
               
        for i in range(2, len(boundary), 2):
            rbound = boundary[i] + 1
            lbound = boundary[i-1] + 1
            gap = rbound - lbound
            if gap % 2 == 0:
                gap = gap//2
                left = np.arange(1, gap+1)*gridres
                right = left[::-1]*-1
                line[lbound:rbound] = np.append(left, right)
            else:
                gap = gap//2
                left = np.arange(1, gap + 1)*gridres
                right = left[::-1]*-1
                adjust = (len(left) + 1)*gridres
                left = np.append(left, adjust)
                line[lbound:rbound] = np.append(left, right)
        rbound = ylen
        lbound = boundary[-1] + 1
        gap = rbound - lbound
        bottom =np.arange(1, gap + 1)*gridres
        line[lbound:rbound] = bottom
    else:
        pass
    return line

##################################################################
##################################################################
##################################################################

for line in range(len(check)):
    check[line] = distance_gridx(check[line], gridres)

for line in range(len(checky)):
    ylen = len(checky[0])
    checky[line] = distance_gridy(checky[line], gridres, ylen)

#plt.imshow(check, interpolation='nearest')
#plt.show()

checky = checky.T
#plt.imshow(checky, interpolation='nearest')
#plt.show()

check = np.ma.masked_where(check == 1., check)
checky = np.ma.masked_where(checky == 1., checky)

gaps = cm.gap(check, checky)
testr = cm.brownian(check,checky, gaps.f1, gaps.f4)
dtestr = cm.drag(check, checky, check, checky, gaps.f1, gaps.f2, gaps.f3, gaps.f4)
tdlvo = cm.DLVO(check, checky)

plt.imshow(tdlvo.EDLy, interpolation='nearest')
plt.colorbar()
plt.show()

plt.imshow(dtestr.drag_x, interpolation='nearest')
plt.colorbar()
plt.show()

plt.imshow(dtestr.drag_y, interpolation='nearest')
plt.colorbar()
plt.show()
