import sys
import subRoutines
import matplotlib.pyplot as plt
import scipy.ndimage.morphology as mph
import numpy as np

plt.close('all')


##################
### user input ###
##################

argsInput = True if len(sys.argv) == 8 else False # checks if an input is given (attention: "main.py" is counted as an input)

N_resolutionSurface = int(sys.argv[1]) if argsInput == True else 200 # number of pixel in x and y direction of surface (in total: N_resolutionSurface² pixel)
sideLengthSurface = int(sys.argv[2]) if argsInput == True else 10 # side length of surface (in certain unit)
N_dots = int(sys.argv[3]) if argsInput == True else 15 # number of dots on surface
N_resolutionTip = int(sys.argv[4]) if argsInput == True else 100 # number of pixel in x and y direction of tip (in total: N_resolutionTip² pixel)
typeTip = sys.argv[5] if argsInput == True else 'parabolic' # type of tip
R_tip_max = int(sys.argv[6]) if argsInput == True else 15 # maximum radius of curvature of tip
m_tip_max = int(sys.argv[7]) if argsInput == True else 0.5 # maximum asymptotic slope of tip (only for hyperbolic tip: theta = 2m)

pixelLength = sideLengthSurface/(N_resolutionSurface-1) # side length of one pixel (in certain unit)

########################
### generate surface ###
########################

surface = subRoutines.genNormSph(N_resolutionSurface, N_dots, sideLengthSurface, pixelLength, 0.06*sideLengthSurface, 0.02*sideLengthSurface, seperateDots=False)

# 3D surface plot
x = np.linspace(0, sideLengthSurface, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
y = x.T
fig = plt.figure()
ax= fig.gca(projection='3d')
ax.set_title('surface')
ax.set_xlabel('x [a.u.]')
ax.set_ylabel('y [a.u.]')
ax.set_zlabel('z [a.u.]')
ax.plot_surface(x, y, surface)

plt.show()

####################################
### calculate surface parameters ###
####################################

# average height of surface
surfaceAverageHeight = np.mean(surface)
imageAverageHeight = []

# standard deviation of height of surface
surfaceStd = np.std(surface)
imageStd = []

# skewness of surface
surfaceSkewness = np.mean((surface - np.mean(surface))**3/ np.std(surface)**3)
imageSkewness = []

# total volume of surface
surfaceVolume = np.trapz(np.trapz(surface, dx=pixelLength), dx=pixelLength)
imageVolume = []

# total surface area
surfaceArea = sideLengthSurface**2
imageMeanArea= []

R_tip = np.linspace(0, R_tip_max, 30) # array with different tip radi

for R_tip_i in R_tip:
    print(R_tip_i)
    
    tip = subRoutines.generateTip(N_resolutionTip, pixelLength, typeTip, R_tip_i, m_tip_max)
    image = mph.grey_dilation(surface, structure=-tip)
    
    imageAverageHeight.append(np.mean(image))
    imageStd.append(np.std(image))
    imageSkewness.append(np.mean((image - np.mean(image))**3/ np.std(image)**3))
    imageVolume.append(np.trapz(np.trapz(image, dx=pixelLength), dx=pixelLength))

####################
### plot results ###
####################

print('Skewness: {:.3f} a.u.'.format(surfaceSkewness))

plt.figure()
plt.title('average height')
plt.xlabel(r'$R_{tip}$')
plt.ylabel(r'$\bar z_{image}/ \bar z_{surface}$')
plt.grid()
plt.plot(R_tip, imageAverageHeight/ surfaceAverageHeight, marker='.')

plt.figure()
plt.title('standard deviation')
plt.xlabel(r'$R_{tip}$')
plt.ylabel(r'$\bar \sigma_{image}/ \bar \sigma_{surface}$')
plt.grid()
plt.plot(R_tip, imageStd/ surfaceStd, marker='.')

plt.figure()
plt.title('skewness')
plt.xlabel(r'$R_{tip}$')
plt.ylabel(r'$R_{SK, image}/ R_{SK, surface}$')
plt.grid()
plt.plot(R_tip, imageSkewness/ surfaceSkewness, marker='.')

plt.figure()
plt.title('volume')
plt.xlabel(r'$R_{tip}$')
plt.ylabel(r'$V_{image}/ V_{surface}$')
plt.grid()
plt.plot(R_tip, imageVolume/ surfaceVolume, marker='.')
