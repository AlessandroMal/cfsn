import subRoutines
import matplotlib.pyplot as plt
import scipy.ndimage.morphology as mph
import numpy as np

plt.close('all')

##### TO DO:
### how to deal with superposition of dots?
### add lengthscale (what is the length of one pixel? important for height profile!)

##################
### user input ###
##################

N_resolution = 200 # pixel in x and y direction (in total: N_resolution² pixel)
N_dots = int(input('N_dots: ')) # number of dots on surface
tipSize = 50 # size of tip (=size of structuring element)


######################################
### generate surface, tip and image###
######################################

#surface = subRoutines.genUnifSph(N_resolution, N_dots, 0.05*N_resolution, 0.25*N_resolution)
surface = subRoutines.genNormSph(N_resolution, N_dots, 0.1*N_resolution, 0.05*N_resolution)
#surface = subRoutines.generateSurface(N_resolution, N_dots)
tip = subRoutines.generateTip(tipSize)
image = mph.grey_dilation(surface, structure=-tip)


###################################
### plot surface, tip and image ###
###################################

plotDict = {'surface': surface, 'tip': tip, 'image': image}
for i in plotDict:
    
    # false color plot
    plt.figure()
    plt.title(i)
    plt.axis('equal')
    plt.pcolor(plotDict[i])
    plt.colorbar()
    
    # 3D surface plot
    if i=='tip':
        x = np.linspace(0, tipSize-1, tipSize) * np.ones([tipSize, tipSize])
    else:
        x = np.linspace(0, N_resolution-1, N_resolution) * np.ones([N_resolution, N_resolution])
    y = x.T
    fig = plt.figure()
    ax= fig.gca(projection='3d')
    ax.set_title(i)
    ax.plot_surface(x, y, plotDict[i])

plt.show()

'''
# plot 
plt.figure()
plt.plot(surface[int(N_resolution/2)])
plt.plot(image[int(N_resolution/2)])
plt.plot(tip[int(tipSize/2)])
'''