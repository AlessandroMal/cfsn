import sys
import subRoutines
import matplotlib.pyplot as plt
import scipy.ndimage.morphology as mph
import numpy as np

plt.close('all')


##################
### user input ###
##################

argsInput = True if len(sys.argv) == 5 else False # checks if an input is given (attention: "main.py" is counted as an input)

N_resolutionSurface = int(sys.argv[1]) if argsInput == True else 200 # number of pixel in x and y direction of surface (in total: N_resolutionSurface² pixel)
sideLengthSurface = int(sys.argv[2]) if argsInput == True else 10 # side length of surface (in certain unit)
N_dots = int(sys.argv[3]) if argsInput == True else 10 # number of dots on surface
N_resolutionTip = int(sys.argv[4]) if argsInput == True else 50 # number of pixel in x and y direction of tip (in total: N_resolutionTip² pixel)

pixelLength = sideLengthSurface/(N_resolutionSurface-1) # side length of one pixel (in certain unit)

######################################
### generate surface, tip and image###
######################################

#surface = subRoutines.generateSurface(N_resolutionSurface, N_dots, sideLengthSurface)
#surface = subRoutines.genUnifSph(N_resolutionSurface, N_dots, sideLengthSurface, pixelLength, 0.05*sideLengthSurface, 0.25*sideLengthSurface)
surface = subRoutines.genNormSph(N_resolutionSurface, N_dots, sideLengthSurface, pixelLength, 0.1*sideLengthSurface, 0)

tip = subRoutines.generateTip(N_resolutionTip, pixelLength)

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
    plt.xlabel('x [a.u.]')
    plt.ylabel('y [a.u.]')
    plt.pcolor(plotDict[i])
    plt.colorbar()
    
    # 3D surface plot
    if i=='tip':
        x = np.linspace(0, pixelLength*N_resolutionTip, N_resolutionTip) * np.ones([N_resolutionTip, N_resolutionTip])
    else:
        x = np.linspace(0, sideLengthSurface, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
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