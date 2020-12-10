import numpy as np
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
import mapsfunction as mf

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

def paramvsTip(surf, pxlen, thres, tipFunc, h, aspectratio_min, aspectratio_max, aspectratio_step):
    aspectratio = np.linspace(aspectratio_min, aspectratio_max, aspectratio_step)
    
    surfParams = calcParams(surf, pxlen, thres)
    imgParams = []
    
    for ar in aspectratio:
        tip = tipFunc(pxlen, h, ar)
        img = mph.grey_dilation(surf, structure=-tip)
        
        imgParams.append(calcParams(img,pxlen,thres))
        
        print('progress: ' + str((np.where(aspectratio==ar)[0][0] + 1)/len(aspectratio) * 100) + '%')
        
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
    
def paramvsNpart(Npx, pxlen, rmin, rmax, N_part_min, N_part_max,
                 tipFunc, h, aspectratio_min, aspectratio_max, aspectratio_step,
                 N_sample, paramFunc, y_label):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    
    z = mf.genFlat(Npx)
    N_part = np.arange(N_part_min, N_part_max+1, 1)
    aspectratio = np.linspace(aspectratio_min, aspectratio_max, aspectratio_step)
    
    plt.figure()
    np.random.seed(123)
    plt_colors = [np.random.random(3) for _ in range(len(aspectratio) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        
        z_param = []
        img_param = []
        
        for N in N_part:
            print('N_part = ' + str(N))
            z_N = mf.genUnifIsolSph(z,pxlen,N,rmin,rmax)
            z_param.append(paramFunc(z_N*pxlen))
                
            for ar in aspectratio:
                tip_ar = tipFunc(pxlen,h,ar)
                img_ar = mph.grey_dilation(z_N, structure=-tip_ar)
                img_param.append(paramFunc(img_ar*pxlen)) 
        
        plt_label = 'surface' if i==0 else '' # visualizza label solo una volta
        plt.plot(N_part, z_param, marker='.', color=plt_colors[-1], label=plt_label)
        for j in range(len(aspectratio)):
            plt_label = 'a.r. = '+str(aspectratio[j])  if i==0 else '' # visualizza label solo una volta
            plt.plot(N_part, img_param[j::len(aspectratio)], marker='.', color=plt_colors[j], label = plt_label)
        
    plt.xlabel(r'$N_{part}$')
    plt.ylabel(y_label)
    plt.grid()   
    plt.legend()
    plt.tight_layout()

def paramvsRpart(Npx, pxlen, R_part_min, R_part_max, R_part_step, N_part,
                 tipFunc, h, aspectratio_min, aspectratio_max, aspectratio_step,
                 N_sample, paramFunc, y_label):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    
    z = mf.genFlat(Npx)
    R_part = np.linspace(R_part_min, R_part_max, R_part_step)
    aspectratio = np.linspace(aspectratio_min, aspectratio_max, aspectratio_step)
    
    plt.figure()
    np.random.seed(123)
    plt_colors = [np.random.random(3) for _ in range(len(aspectratio) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        
        z_param = []
        img_param = []
        
        for R in R_part:
            print('R_part = ' + str(R))
            z_N = mf.genUnifIsolSph(z,pxlen,N_part,R,R)
            z_param.append(paramFunc(z_N*pxlen))
                
            for ar in aspectratio:
                tip_ar = tipFunc(pxlen,h,ar)
                img_ar = mph.grey_dilation(z_N, structure=-tip_ar)
                img_param.append(paramFunc(img_ar*pxlen)) 
        
        plt_label = 'surface' if i==0 else '' # visualizza label solo una volta
        plt.plot(R_part, z_param, marker='.', color=plt_colors[-1], label=plt_label)
        for j in range(len(aspectratio)):
            plt_label = 'a.r. = '+str(aspectratio[j])  if i==0 else '' # visualizza label solo una volta
            plt.plot(R_part, img_param[j::len(aspectratio)], marker='.', color=plt_colors[j], label = plt_label)
        
    plt.xlabel(r'$R_{part} [nm]$')
    plt.ylabel(y_label)
    plt.grid()   
    plt.legend()
    plt.tight_layout()