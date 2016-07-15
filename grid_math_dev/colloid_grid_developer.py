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

gridres = 1e-8

################################################################
### Definitions to set grid distances based on vector forces ###
### Pos ++ Up and Right, Neg. -- Down and Left               ###
################################################################

def distance_gridx(line, vline, gridres, solid):
    boundary= np.where(np.abs(np.diff(line)) >= 1)[0]
    
    try:
        for i in range(1,len(boundary),2):
            rbound = boundary[i] + 1
            lbound = boundary[i-1] + 1
            gap = rbound - lbound

            if gap % 2 == 0:
                gap = gap//2
                left = np.arange(1, gap + 1)*gridres#*-1*gridres
                right = left[::-1]#*-1
                line[lbound:rbound] = np.append(left, right)

                left = np.ones(gap)*-1
                right = np.ones(gap)
                vline[lbound:rbound] = np.append(left, right)

            else:
                gap = gap//2
                left = np.arange(1, gap + 1)*gridres#-1*gridres
                right = left[::-1]#*-1
                adjust = (len(left) + 1)*gridres#-1*gridres
                left = np.append(left, adjust)
                line[lbound:rbound] = np.append(left, right)

                left = np.ones(gap + 1)*-1
                right = np.ones(gap)
                vline[lbound:rbound] = np.append(left, right)
                    
    except IndexError:
        print 'volume does not percolate'
    line[line == solid] = np.nan
    
    return line, vline


def distance_gridy(line, vline, gridres, ylen, solid):
    boundary= np.where(np.abs(np.diff(line)) >= 1)[0]

    if len(boundary) > 0:
        rbound = boundary[0] + 1
        lbound = 0
        top = np.arange(rbound, lbound, -1)*gridres#*-1
        line[lbound:rbound] = top

        vtop = np.ones(len(top))*-1
        vline[lbound:rbound] = vtop
               
        for i in range(2, len(boundary), 2):
            rbound = boundary[i] + 1
            lbound = boundary[i-1] + 1
            gap = rbound - lbound
            if gap % 2 == 0:
                gap = gap//2
                left = np.arange(1, gap+1)*gridres
                right = left[::-1]#*-1
                line[lbound:rbound] = np.append(left, right)
                # create our vector array!
                left = np.ones(gap)*-1
                right = np.ones(gap)
                vline[lbound:rbound] = np.append(left, right)
                
                
            else:
                gap = gap//2
                left = np.arange(1, gap + 1)*gridres
                right = left[::-1]#*-1
                adjust = (len(left) + 1)*gridres
                left = np.append(left, adjust)
                line[lbound:rbound] = np.append(left, right)

                left = np.ones(gap + 1)*-1
                right = np.ones(gap)
                vline[lbound:rbound] = np.append(left, right)
                
        rbound = ylen
        lbound = boundary[-1] + 1
        gap = rbound - lbound
        bottom =np.arange(1, gap + 1)*gridres
        line[lbound:rbound] = bottom

        bottom = np.ones(gap)
        vline[lbound:rbound] = bottom
    else:
        pass
    line[line == solid] = np.nan
    line[line == 0.] = np.nan
    return line, vline
    
def create_vector_array(img, solid):
    '''
    creates an x-dir copy and y-dir copy of the image domain for vector
    directions to populate.
    '''
    vimgx = np.copy(img)
    vimgy = np.copy(img.T)
    vimgx[vimgx == solid] = np.nan
    vimgy[vimgy == solid] = np.nan
    return vimgx, vimgy

##################################################################
##################################################################
##################################################################

vector_x, vector_y = create_vector_array(check, True)

for line in range(len(check)):
    check[line], vector_x[line] = distance_gridx(check[line], vector_x[line], gridres, solid=True)

for line in range(len(checky)):
    ylen = len(checky[0])
    checky[line], vector_y[line] = distance_gridy(checky[line], vector_y[line], gridres, ylen, solid=True)

#plt.imshow(check, interpolation='nearest')
#plt.show()

checky = checky.T
vector_y = vector_y.T
#plt.imshow(checky, interpolation='nearest')
#plt.show()

plt.imshow(vector_x, interpolation='nearest')
plt.colorbar()
plt.show()

plt.imshow(vector_y, interpolation='nearest')
plt.colorbar()
plt.show()


'''
plt.imshow(checky, interpolation='nearest')
plt.colorbar()
plt.show()

gaps = cm.Gap(check, checky)
testr = cm.Brownian(check,checky, gaps.f1, gaps.f4)
dtestr = cm.Drag(check, checky, check, checky, gaps.f1, gaps.f2, gaps.f3, gaps.f4)
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
'''
