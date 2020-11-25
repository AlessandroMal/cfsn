import numpy as np
from random import uniform

def generateSurface(N_resolution, N_dots):
    
    x = np.linspace(0, N_resolution-1, N_resolution) * np.ones([N_resolution, N_resolution])
    y = x.T
    z = np.zeros([N_resolution, N_resolution])
    
    for i in range(N_dots):
        
        R = uniform(0.05*N_resolution, 0.25*N_resolution)
        x0 = uniform(0, N_resolution)
        y0 = uniform(0, N_resolution)
        
        dot = np.sqrt(R**2 - (x - x0)**2 - (y - y0)**2) + R
        dot = np.nan_to_num(dot)

        z = z + dot
    
    return z
   

def generateTip(tipSize):

    x = np.linspace(0, tipSize-1, tipSize) * np.ones([tipSize, tipSize])
    y = x.T
    
    tip = 0.3*(x - tipSize/2)**2 + 0.3*(y - tipSize/2)**2

    return tip