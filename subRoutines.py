import numpy as np
from random import uniform

def generateSurface(N_resolutionSurface, N_dots, sideLength):
    
    x = np.linspace(0, sideLength, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
    y = x.T
    z = np.zeros([N_resolutionSurface, N_resolutionSurface])
    
    for i in range(N_dots):
        
        R = uniform(0.05*sideLength, 0.25*sideLength)
        x0 = uniform(0, sideLength)
        y0 = uniform(0, sideLength)
        
        dot = np.sqrt(R**2 - (x - x0)**2 - (y - y0)**2) + R
        dot = np.nan_to_num(dot)

        z = z + dot
    
    return z
   

def generateTip(N_resolutionTip, pixelSideLength):
    
    tipSideLength = pixelSideLength*N_resolutionTip
    
    x = np.linspace(0, tipSideLength, N_resolutionTip) * np.ones([N_resolutionTip, N_resolutionTip])
    y = x.T
    
    tip = 2*(x - tipSideLength/2)**2 + 2*(y - tipSideLength/2)**2

    return tip