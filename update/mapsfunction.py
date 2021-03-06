import numpy as np
from random import uniform
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.ndimage

def genFlat(Npx):
    z = np.zeros([Npx, Npx])
    return z


def genNormNoise(z,pxlen,var,Ltile):
    l=int(Ltile/pxlen)
    ntiles=int(len(z)/l)
    
    if len(z)%l!=0:
        ntiles+=1
    
    for i in range(ntiles):
        for j in range(ntiles):
            r=np.random.normal(0, var)
            for x in range(l):
                for y in range(l):
                    if j*l+y<len(z) and i*l+x<len(z):
                        z[j*l+y,i*l+x]+=r
    return z

'''
    if len(z)%l!=0:
        for i in range(ntiles):
            r=np.random.normal(0, var)
            for x in range(l):
                for y in range(l-len(z)%l):
                    z[ntiles*l+y,i*l+x]+=r
        
        for i in range(ntiles):
            r=np.random.normal(0, var)
            for y in range(l):
                for x in range(l-len(z)%l):
                    z[i*l+y,ntiles*l+x]+=r
                    
        for x in range(l-len(z)%l):
            for y in range(l-len(z)%l):
 '''                          
    
def plotview(z,pxlen,theta,phi):
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    X = np.linspace(0,pxlen*len(z),num=len(z))
    Y = np.linspace(0,pxlen*len(z),num=len(z))
    X, Y = np.meshgrid(X, Y)

    ax.plot_surface(X,Y,z, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax.view_init(theta,phi)
    
    #ax.set_xlim(0, 10)
    #ax.set_ylim(0, 10)
    #ax.set_zlim(0, 10)
    ax.set_xlabel('X (nm)')
    ax.set_ylabel('Y (nm)')
    ax.set_zlabel('Z (nm)')
    
def plotfalsecol(z,pxlen):
    plt.figure()
    #plt.title()
    plt.axis('equal')
    
    X = np.linspace(0,pxlen*len(z),num=len(z))
    Y = np.linspace(0,pxlen*len(z),num=len(z))
    X, Y = np.meshgrid(X, Y)
    
    plt.pcolormesh(z)
    clb = plt.colorbar()    
    
    clb.set_label('Z (nm)')
    plt.xlabel('X (nm)')
    plt.ylabel('Y (nm)')
    
def genSphere(z,pxlen,centre,R):
    if(len(centre)!=len(R)*2):
        print("input genSphere non validi")
        return 0 #controllo input
    
    for i in range(len(R)):
        
        x0 = centre[2*i]
        y0 = centre[2*i+1]
        
        for x in range(len(z)):
            for y in range(len(z)):
                if R[i]**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x] += np.sqrt(R[i]**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R[i]
    return z

def genUnifSph(z,pxlen,Nsph,rmin,rmax, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
        
    for i in range(Nsph):
        R  = uniform(rmin, rmax)
        x0 = uniform(xmin, xmax)
        y0 = uniform(ymin, ymax)
        
        for x in range(len(z)):
            for y in range(len(z)):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

    return z

def genUnifIsolSph(z,pxlen,Nsph,rmin,rmax, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
    
    n=0
    xyr=[]

    while n<Nsph:
        R  = uniform(rmin, rmax)
        x0 = uniform(xmin, xmax)
        y0 = uniform(ymin, ymax)
        
        b=True
        for i in range(int(len(xyr)/3)):
           if R+xyr[3*i+2]>np.sqrt((x0-xyr[3*i])**2+(y0-xyr[3*i+1])**2):
               b=False
        
        if b==True:
            n+=1
            xyr.append(x0)
            xyr.append(y0)
            xyr.append(R)
            
            for x in range(len(z)):
                for y in range(len(z)):
                    if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                        z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

    return z

def genNormSph(z,pxlen,Nsph,av,var, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))

    for i in range(Nsph):
#uso una tecnica di rigetto per evitare R negativi se metto
#media e varianza pericolose
#Ale
        while 1>0:
            R = np.random.normal(av, var)
            if R>0:
                break
        
        x0 = uniform(xmin,xmax)
        y0 = uniform(ymin,ymax)
        
        for x in range(len(z)):
            for y in range(len(z)):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

    return z

def genParabolicTip(pxlen,h,aspectratio):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    l=h/aspectratio
    px=int(l/pxlen)
    z=np.zeros([px,px])
    
    for x in range(len(z)):
        for y in range(len(z)):
            if (l/2)**2 - (x*pxlen - px/2*pxlen)**2 - (y*pxlen - px/2*pxlen)**2 > 0:
                z[y,x]=4*aspectratio/l*((x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2)
            else:
                z[y,x]=h
    return z

def genPyramidTip(pxlen,h,aspectratio):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    l=h/aspectratio
    px=int(l/pxlen)
    z=np.zeros([px,px])
    
    for x in range(len(z)):
        for y in range(len(z)):
                z[y,x]=2*aspectratio*max(abs(x*pxlen-px/2*pxlen), abs(y*pxlen-px/2*pxlen))
            
    return z

def genSemisphTip(pxlen,h,aspectratio):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    l=h/aspectratio
    px=int(l/pxlen)
    z=np.ones([px,px])*h
    R=l/2
    
    for x in range(len(z)):
        for y in range(len(z)):
            r2= (x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2
            if r2 < R**2:
                z[y,x]=R*(1-np.sqrt(1-r2/R**2))

    return z

def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    objInd = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj = [z[i] for i in objInd]
    print('\nidentObj ha trovato ' + str(len(obj)) + ' particelle\n')

    return obj

def skew(z,pxlen):
    z_skew = np.mean((z*pxlen - np.mean(z*pxlen))**3/ np.std(z*pxlen)**3)
    return z_skew

def V(z,pxlen):
    somma=0
    
    for x in range(len(z)):
        for y in range(len(z)):
            somma+=z[y,x]
    
    return pxlen*pxlen*somma

def specArea(z,pxlen):
    A=0
    for x in range(len(z)-1):
        for y in range(len(z)-1):
            if x==0: a=np.array([pxlen,0,z[y,x+1]-z[y,x]])
            else: a=a*(-1)
            b=np.array([0,pxlen,z[y+1,x]-z[y,x]])
            A+=np.linalg.norm(np.cross(a,b))/2
            
            a=np.array([-pxlen,0,z[y+1,x]-z[y+1,x+1]])
            b=np.array([0,-pxlen,z[y,x+1]-z[y+1,x+1]])
            A+=np.linalg.norm(np.cross(a,b))/2
            
    return A/(pxlen*pxlen*len(z)*len(z))

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

def paramTipDepend(surf, pxlen, thres, tipType, h, aspectratio_min, aspectratio_max, step):
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

def paramDepend(Npx, pxlen, rmin, rmax, N_part_min, N_part_max,
                tipType, h, aspectratio_min, aspectratio_max, aspectratio_step,
                N_sample, calcParam, y_label):
    
    z = genFlat(Npx)
    N_part = np.arange(N_part_min, N_part_max+1, 1)
    aspectratio = np.linspace(aspectratio_min, aspectratio_max, aspectratio_step)
    
    plt.figure()
    plt_colors = [np.random.random(3) for _ in range(len(aspectratio) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        
        z_param = []
        img_param = []
        
        for N in N_part:
            print('N = ', N)
            z_N = genUnifIsolSph(z,pxlen,N,rmin,rmax)
            z_param.append(calcParam(z_N*pxlen))
                
            for ar in aspectratio:
                print('ar = ', ar)
                tip_ar = tipType(pxlen,h,ar)
                img_ar = mph.grey_dilation(z_N, structure=-tip_ar)
                img_param.append(calcParam(img_ar*pxlen)) 
        
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
 