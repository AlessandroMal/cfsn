import subRoutines
import matplotlib.pyplot as plt
import scipy.ndimage.morphology as mph
import numpy as np

plt.close('all')


##################
### user input ###
##################

N_resolutionSurface = 200 # number of pixel in x and y direction of surface (in total: N_resolutionSurface² pixel)
sideLength = 10 # side length of surface (in certain unit)
N_dots = 5 # number of dots on surface
N_resolutionTip = 50 # numberof pixel in x and y direction of tip (in total: N_resolutionTip² pixel)

pixelSideLength = sideLength/N_resolutionSurface # side length of one pixel 

######################################
### generate surface, tip and image###
######################################

surface = subRoutines.generateSurface(N_resolutionSurface, N_dots, sideLength)
tip = subRoutines.generateTip(N_resolutionTip, pixelSideLength)
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
        x = np.linspace(0, pixelSideLength*N_resolutionTip, N_resolutionTip) * np.ones([N_resolutionTip, N_resolutionTip])
    else:
        x = np.linspace(0, sideLength, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
    y = x.T
    fig = plt.figure()
    ax= fig.gca(projection='3d')
    ax.set_title(i)
    ax.plot_surface(x, y, plotDict[i])

plt.show()

'''
# plot 
plt.figure()
plt.plot(surface[int(N_resolutionSurface/2)])
plt.plot(image[int(N_resolutionSurface/2)])
plt.plot(tip[int(N_resolutionTip/2)])
'''