import numpy as np
from random import uniform

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
   

def generateTip(N_resolutionTip, pixelLength):
    
    tipSideLength = pixelLength*N_resolutionTip
    
    x = np.linspace(0, tipSideLength, N_resolutionTip) * np.ones([N_resolutionTip, N_resolutionTip])
    y = x.T
    
    tip = 2*(x - tipSideLength/2)**2 + 2*(y - tipSideLength/2)**2

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
    z = np.zeros([N_res, N_res])

    for i in range(N_dots):
#uso una tecnica di rigetto per evitare R negativi se metto
#media e varianza pericolose
#Ale
        while 1>0:
            R = np.random.normal(av, var)
            if R>0:
                break
        
        x0 = uniform(0, l_side)
        y0 = uniform(0, l_side)
        dot = np.zeros([N_res, N_res])
        
        for x in range(N_res):
            for y in range(N_res):
                if R**2 - (x*pixelLength - x0)**2 - (y*pixelLength - y0)**2>0:
                    dot[x,y] += np.sqrt(R**2 - (x*pixelLength - x0)**2 - (y*pixelLength - y0)**2) + R
            
        z = z + dot
    return z