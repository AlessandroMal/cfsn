import numpy as np
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt

def skew(z,pxlen):
    std=np.std(z*pxlen)
    if std!=0: z_skew = np.mean((z*pxlen - np.mean(z*pxlen))**3/ std**3)
    else: z_skew=0
    return z_skew

def V(z,pxlen):
    somma=0
    
    for x in range(len(z)):
        for y in range(len(z)):
            somma+=z[y,x]
    
    return pxlen*pxlen*somma

def specArea(z,pxlen):
#first order specific area
    A=0
    for x in range(len(z)-1):
        for y in range(len(z)-1):
            A+=pxlen*np.linalg.norm(np.array([-z[y,x+1]+z[y,x], -z[y+1,x]+z[y,x], pxlen]))/2
            A+=pxlen*np.linalg.norm(np.array([z[y+1,x]-z[y+1,x+1], z[y,x+1]-z[y+1,x+1], pxlen]))/2
            '''
            if y==0: a=np.array([pxlen,0,z[y,x+1]-z[y,x]])
            else: a=a*(-1)
            b=np.array([0,pxlen,z[y+1,x]-z[y,x]])
            A+=np.linalg.norm(np.cross(a,b))/2 #fare il prodotto vettore non conviene computazionalmente
            
            a=np.array([-pxlen,0,z[y+1,x]-z[y+1,x+1]])
            b=np.array([0,-pxlen,z[y,x+1]-z[y+1,x+1]])
            A+=np.linalg.norm(np.cross(a,b))/2
            '''
    return A/( pxlen*(len(z)-1) )**2

def coverage(z, thres):
    N_tot = len(z)**2
    N = np.sum(z>thres)
    cov = N/N_tot
    return cov

def calcParams(z,pxlen,thres):
    params = {'mean': np.mean(z*pxlen),
              'std': np.std(z*pxlen),
              'skew': skew(z,pxlen),
              'V': V(z,pxlen),
              'specArea': specArea(z,pxlen),
              'coverage': coverage(z,thres)}
    return params

def paramvsTip(surf, pxlen, thres, tipType, h, aspectratio_min, aspectratio_max, step):
    aspectratio = np.linspace(aspectratio_min, aspectratio_max, step)
    
    surfParams = calcParams(surf,pxlen,thres)
    imgParams = []
    
    for i in aspectratio:
        tip = globals()[tipType](pxlen,h,i)
        img = mph.grey_dilation(surf, structure=-tip)
        
        imgParams.append(calcParams(img,pxlen,thres))
    
    plt.figure()
    
    for i in imgParams[0]:
        param = []
        for j in range(len(imgParams)):
            param.append(imgParams[j][i]/ surfParams[i])
        plt.plot(aspectratio, param, marker='.', label=i)
        
    plt.xlabel(r'$aspectratio \propto 1/ R_{tip}$')
    plt.ylabel('rel. values (image/surface)')
    plt.grid()   
    plt.legend()