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