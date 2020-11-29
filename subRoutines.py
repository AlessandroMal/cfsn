import numpy as np
from random import uniform
import scipy.ndimage

def generateSurface(N_resolutionSurface, N_dots, sideLengthSurface):
    
    x = np.linspace(0, sideLengthSurface, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
    y = x.T
    z = np.zeros([N_resolutionSurface, N_resolutionSurface])

    for i in range(N_dots):
        
        R = uniform(0.05*sideLengthSurface, 0.25*sideLengthSurface)
        x0 = uniform(0, sideLengthSurface)
        y0 = uniform(0, sideLengthSurface)
        
        dot = np.sqrt(R**2 - (x - x0)**2 - (y - y0)**2) + R
        dot = np.nan_to_num(dot)

        z = z + dot
    
    return z
   

def generateTip(N_resolutionTip, pixelLength, typeTip, R, m):
    
    tipSideLength = pixelLength*N_resolutionTip
    
    x = np.linspace(0, tipSideLength, N_resolutionTip) * np.ones([N_resolutionTip, N_resolutionTip])
    y = x.T
    
    if typeTip == 'parabolic':
        tip = ((x - tipSideLength/2)**2 + (y - tipSideLength/2)**2)/ (2*R) # t(x) = 1/2 * x²/R
    if typeTip == 'hyperbolic':
        tip = 1/ m * np.sqrt(R**2/m**2 + (x - tipSideLength/2)**2 + (y - tipSideLength/2)**2) - R/ m**2 # t(x) = 1/m * sqrt(R²/m² + x²) - R/m²

    return tip


def genUnifSph(N_res, N_dots, l_side, pixelLength,  rmin, rmax):
#questa genera la mappa senza creare due matrici in più, usa già la griglia che z ha nella sua definizione naturale
#può essere utile e più rapido se ho tanti pixel (N_res alto)
#mi pare più efficace di definire una doppia griglia, ditemi voi
#in più passo l'intervallo di generazione del raggio dal main
#Alessandro
    z = np.zeros([N_res, N_res])

    for i in range(N_dots):
        
        R = uniform(rmin, rmax)
        x0 = uniform(0, l_side)
        y0 = uniform(0, l_side)
        dot = np.zeros([N_res, N_res])
        
        for x in range(N_res):
            for y in range(N_res):
                if R**2 - (x*pixelLength - x0)**2 - (y*pixelLength - y0)**2 > 0:
                    dot[x,y] += np.sqrt(R**2 - (x*pixelLength - x0)**2 - (y*pixelLength - y0)**2) + R
#l'if è solo per togliere il warning del nan di python mi dava fastidio            
        z = z + dot
    return z


def genNormSph(N_res, N_dots, l_side, pixelLength, av, var):
#come prima ma raggi con distribuzione normale
#Alessandro
    
    while 1>0:
        R = np.random.normal(av, var)
        if R>0:
            break
#uso una tecnica di rigetto per evitare R negativi se metto
#media e varianza pericolose
#Ale
    
    x0 = []
    y0 = []
    
    # find sufficiently seperated points on surface to prevent overlap of spheres
    while len(x0) < N_dots:
        x0_new = uniform(0, l_side)
        y0_new = uniform(0, l_side)
        
        validPoint = True
        for j in range(len(x0)):
            if np.sqrt((x0_new - x0[j])**2 + (y0_new - y0[j])**2) < 3*R: # SPECIFY DEGREE OF SEPERATION
                validPoint=False
        
        if validPoint == True:
            x0.append(x0_new)
            y0.append(y0_new)
    
    z = np.zeros([N_res, N_res])
    
    # place spheres on surface
    for i in range(N_dots):

        dot = np.zeros([N_res, N_res])
        
        for x in range(N_res):
            for y in range(N_res):
                if R**2 - (x*pixelLength - x0[i])**2 - (y*pixelLength - y0[i])**2>0:
                    dot[x,y] += np.sqrt(R**2 - (x*pixelLength - x0[i])**2 - (y*pixelLength - y0[i])**2) + R
            
        z = z + dot
        
    return z


def identSph(image, N_dots):

    imageNonZero = (image>0) # true if element is not zero
    labeledObjects = scipy.ndimage.label(imageNonZero)[0] # label single objects with 1, 2, 3,...
    singleObjectsIndexes = scipy.ndimage.find_objects(labeledObjects) # find indexes of single objects
    
    if len(singleObjectsIndexes) != N_dots:
        raise Exception('identSph couldn\'t find all objects!')
    
    singleObjects = [image[i] for i in singleObjectsIndexes] # put subarrays of single objects in one list

    return singleObjects