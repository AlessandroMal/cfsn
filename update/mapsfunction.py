import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from matplotlib import cm

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
    #plt.axis('equal')
    
    X = np.linspace(0,pxlen*len(z),num=len(z))
    Y = np.linspace(0,pxlen*len(z),num=len(z))
    X, Y = np.meshgrid(X, Y)
    
    plt.pcolormesh(z)
    plt.colorbar()    
    
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
           if  R+xyr[3*i+2]>np.sqrt((x0-xyr[3*i])**2+(y0-xyr[3*i+1])**2):
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