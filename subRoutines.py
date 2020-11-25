import numpy as np
from random import uniform

def generateSurface(N_resolutionSurface, N_dots, sideLength):
    
    x = np.linspace(0, sideLength, N_resolutionSurface) * np.ones([N_resolutionSurface, N_resolutionSurface])
    y = x.T
<<<<<<< HEAD
    z = np.zeros([N_resolution, N_resolution])

=======
    z = np.zeros([N_resolutionSurface, N_resolutionSurface])
    
>>>>>>> 36ff4f129e582131ec365675ff64dc5c5e17a944
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

def genUnifSph(N_res, N_dots, rmin, rmax):
#questa genera la mappa senza creare due matrici in più, usa già la griglia che z hanella sua definizione naturale
#può essere utile e più rapido se ho tanti pixel (N_res alto)
#mi pare più efficace di definire una doppia griglia, ditemi voi
#in più passo l'intervallo di generazione del raggio dal main
#Alessandro
    z = np.zeros([N_res, N_res])

    for i in range(N_dots):
        
        R = uniform(rmin, rmax)
        x0 = uniform(0, N_res)
        y0 = uniform(0, N_res)
        dot=np.zeros([N_res, N_res])
        
        for x in range(N_res):
            for y in range(N_res):
                if R**2 - (x - x0)**2 - (y - y0)**2>0:
                    dot[x,y] += np.sqrt(R**2 - (x - x0)**2 - (y - y0)**2) + R
#l'if è solo per togliere il warning del nan di python mi dava fastidio            
        z = z + dot
    return z

def genNormSph(N_res, N_dots, av, var):
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
        
        x0 = uniform(0, N_res)
        y0 = uniform(0, N_res)
        dot=np.zeros([N_res, N_res])
        
        for x in range(N_res):
            for y in range(N_res):
                if R**2 - (x - x0)**2 - (y - y0)**2>0:
                    dot[x,y] += np.sqrt(R**2 - (x - x0)**2 - (y - y0)**2) + R
            
        z = z + dot
    return z