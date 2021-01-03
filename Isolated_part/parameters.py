import numpy as np
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
import mapsfunction as mf

def skew(z,pxlen):
    std=h_std(z,pxlen)
    if std!=0: z_skew = (h_av(z**3,pxlen)-3*h_av(z,pxlen)*std**2-h_av(z,pxlen)**3)/ std**3
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
    return A/( pxlen*(len(z)-1) )**2

def coverage(z, thres):
    N_tot = len(z)**2
    N = np.sum(z>thres)
    cov = N/N_tot
    return cov

def h_av(z,pxlen):
    summ=0
    for x in range(len(z)-1):
        for y in range(len(z)-1):
            summ+=z[y][x]
    return summ/len(z)**2

def h_std(z,pxlen):
    summ=0
    sum2=0
    for x in range(len(z)-1):
        for y in range(len(z)-1):
            summ+=z[y][x]
            sum2+=z[y][x]**2
    summ/=len(z)**2
    sum2/=len(z)**2
    return np.sqrt(sum2-summ**2)

def h_max(z,n):
    top=[]
    for x in range(len(z)-1):
        for y in range(len(z)-1):
            if y<n and x==0:
                top.append(z[y,x])
                if y==n-1: 
                    top=sorted(top)
            else:
                if z[y,x]>top[0]:
                    top[0]=z[y,x]
                    top=sorted(top)
    if n%2==0: return (top[int(n/2)]+top[int(n/2)-1])/2
    else: return top[int(n/2)]

def calcParams(z,pxlen,thres): 
    params = {'mean': h_av(z,pxlen),
              'std': h_std(z,pxlen),
              'skew': skew(z,pxlen),
              'V': V(z,pxlen),
              'specArea': specArea(z,pxlen),
              'coverage': coverage(z,thres)}
    return params
    
def paramvsTip(surf, pxlen, thres, tipFunc, h, rtipmin, rtipmax, rtipstep):
    R = np.linspace(rtipmin, rtipmax, rtipstep)
    surfParams = calcParams(surf, pxlen, thres)
    imgParams = []
    
    for rtip in R:
        tip = tipFunc(pxlen, h, r=rtip)
        img = mph.grey_dilation(surf, structure=-tip)
        imgParams.append(calcParams(img,pxlen,thres))
        for i in imgParams[0]:
            imgParams[-1][i]=imgParams[-1][i]/surfParams[i]
        print('progress: ' + str((np.where(R==rtip)[0][0] + 1)/len(R) * 100) + '%')
        
    imgParams.append('image val/surface val')
    Rlist = []
    Rlist.append(R)
    Rlist.append(r'$R_{tip}$')
    return Rlist, imgParams

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
            z_param.append(paramFunc(z_N))
                
            for ar in aspectratio:
                tip_ar = tipFunc(pxlen,h,r=ar)
                img_ar = mph.grey_dilation(z_N, structure=-tip_ar)
                img_param.append(paramFunc(img_ar)) 
        
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
                 tipFunc, h, rtipmin, rtipmax, rtipstep,
                 N_sample, paramFunc, y_label):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    
    z = mf.genFlat(Npx)
    R_part = np.linspace(R_part_min, R_part_max, R_part_step)
    Rtip = np.linspace(rtipmin, rtipmax, rtipstep)
    
    plt.figure()
    np.random.seed(123)
    plt_colors = [np.random.random(3) for _ in range(len(Rtip) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        
        z_param = []
        img_param = []
        
        for R in R_part:
            print('R_part = ' + str(R), '    Ncp = '+str(mf.Ncp(pxlen,Npx,R)))
            z_N = mf.genUnifIsolSph(z,pxlen,N_part,R,R)
            z_param.append(paramFunc(z_N*pxlen))
                
            for radius in Rtip:
                tip_ar = tipFunc(pxlen,h,r=radius)
                img_ar = mph.grey_dilation(z_N, structure=-tip_ar)
                img_param.append(paramFunc(img_ar*pxlen)) 
        
        plt_label = 'surface' if i==0 else '' # visualizza label solo una volta
        plt.plot(R_part, z_param, marker='.', color=plt_colors[-1], label=plt_label)
        for j in range(len(Rtip)):
            plt_label = r'$R_{tip}=$'+str(Rtip[j])+r'$nm$' if i==0 else '' # visualizza label solo una volta
            plt.plot(R_part, img_param[j::len(Rtip)], marker='.', color=plt_colors[j], label = plt_label)
    
    plt.title(r'$res=$'+str(pxlen)+r'$nm, n=$'+str(N_part/pxlen/Npx)+r'$nm^{-2}$')
    plt.xlabel(r'$R_{part} [nm]$')
    plt.ylabel(y_label)
    plt.grid()   
    plt.legend()
    plt.tight_layout()    

def plotParams(x,y):
    xlabel=x[len(x)-1]
    x.pop(len(x)-1)
    ylabel=y[len(y)-1]
    y.pop(len(y)-1)
    
    plt.figure()
    
    for i in y[0]:
        param = []
        for j in range(len(y)):
            param.append(y[j][i])
        plt.plot(x[0], param, marker='.', label=i)
    
#    plt.title()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid()   
    plt.legend()
