import numpy as np
from random import uniform
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
    
def genSphere(z,pxlen,centre,R):
    if(len(centre)!=len(R)*2):
        print("input genSphere non validi")
        return 0 #controllo input
    
    for i in range(len(R)):
        
        x0 = centre[2*i]
        y0 = centre[2*i+1]
        
        lwrx=int((x0-R[i])/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
        if lwrx<0: lwrx=0            #tutta la mappa, importante se ho
        uprx=int((x0+R[i])/pxlen)+1  #alta risoluzione
        if uprx>len(z): uprx=len(z)
        lwry=int((y0-R[i])/pxlen)-1
        if lwry<0: lwry=0
        upry=int((y0+R[i])/pxlen)+1
        if upry>len(z): upry=len(z)
        
        for x in range(lwrx,uprx):
            for y in range(lwry,upry):
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
        
        lwrx=int((x0-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
        if lwrx<0: lwrx=0            #tutta la mappa, importante se ho
        uprx=int((x0+R)/pxlen)+1  #alta risoluzione
        if uprx>len(z): uprx=len(z)
        lwry=int((y0-R)/pxlen)-1
        if lwry<0: lwry=0
        upry=int((y0+R)/pxlen)+1
        if upry>len(z): upry=len(z)
        
        for x in range(lwrx,uprx):
            for y in range(lwry,upry):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

    return z

def genHexSpikes(z, pxlen, h, aspectratio, d, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
    
    R=h/aspectratio/2
    b=True
    while xmin<xmax:
        if b: y0=ymin
        else: y0=ymin+d/2/pxlen
        while y0<ymax:
            lwrx=int((xmin-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
            if lwrx<0: lwrx=0            #tutta la mappa, importante se ho
            uprx=int((xmin+R)/pxlen)+1  #alta risoluzione
            if uprx>len(z): uprx=len(z)
            lwry=int((y0-R)/pxlen)-1
            if lwry<0: lwry=0
            upry=int((y0+R)/pxlen)+1
            if upry>len(z): upry=len(z)
        
            for x in range(lwrx,uprx):
                for y in range(lwry,upry):
                    if R**2 - (x*pxlen - xmin)**2 - (y*pxlen - y0)**2 > 0:
                        z[y,x]+=2*aspectratio*(R-np.sqrt( (x*pxlen - xmin)**2 + (y*pxlen - y0)**2))
            
            y0+=d/pxlen
        xmin+=d*np.sqrt(3)/2/pxlen
        b=not(b)
    
    return z

def genUnifIsolSph(z,pxlen,Nsph,rmin,rmax, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
    
    extraline=0 #calcolo numero di sfere per close packing
    d=rmin+rmax
    l=(int(pxlen*len(z)/d)+1)*2
    col=int(pxlen*len(z)/(np.sqrt(3)*d))+1
    if pxlen*len(z)%d<d/2: l=l-1
    if pxlen*len(z)%(np.sqrt(3)*d)>np.sqrt(3)/2*d: extraline=int(l/2)
    #print("number of particle with r=",d/2,"for hexagonal cp: ", l*col+extraline)
    
    if Nsph>=l*col+extraline:
        print("Error: too many particles for this map")
        return 0
        
    n=0
    xyr=[]
    while n<Nsph:
        if n<0.5*l*col+extraline:
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
                
                lwrx=int((x0-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
                if lwrx<0: lwrx=0         #tutta la mappa, importante se ho
                uprx=int((x0+R)/pxlen)+1  #alta risoluzione
                if uprx>len(z): uprx=len(z)
                lwry=int((y0-R)/pxlen)-1
                if lwry<0: lwry=0
                upry=int((y0+R)/pxlen)+1
                if upry>len(z): upry=len(z)
            
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                            z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R
        else:
            remaining=[]
            for x in range(len(z)):
                for y in range(len(z)): #ciclo sui punti della mappa
                    b=False
                    for i in range(int(len(xyr)/3)): #ciclo su sfere già messe
                        if rmax+xyr[3*i+2]>np.sqrt((x*pxlen-xyr[3*i])**2+(y*pxlen-xyr[3*i+1])**2):
                            b=True
                            break
                    if b:
                        break
                    else:
                        remaining.append(x)
                        remaining.append(y)
                        
            if len(remaining)>0:
                n+=1
                rndp = int(uniform(0,len(remaining)/2))
                x0 = remaining[rndp*2]
                y0 = remaining[rndp*2+1]
                remaining.clear()
                R  = uniform(rmin, rmax)
                
                lwrx=int((x0-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
                if lwrx<0: lwrx=0         #tutta la mappa, importante se ho
                uprx=int((x0+R)/pxlen)+1  #alta risoluzione
                if uprx>len(z): uprx=len(z)
                lwry=int((y0-R)/pxlen)-1
                if lwry<0: lwry=0
                upry=int((y0+R)/pxlen)+1
                if upry>len(z): upry=len(z)
                
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                            z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R
            else: break
    if(n!=Nsph): print("space over: ",n," particles generated")
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
        
        lwrx=int((x0-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
        if lwrx<0: lwrx=0         #tutta la mappa, importante se ho
        uprx=int((x0+R)/pxlen)+1  #alta risoluzione
        if uprx>len(z): uprx=len(z)
        lwry=int((y0-R)/pxlen)-1
        if lwry<0: lwry=0
        upry=int((y0+R)/pxlen)+1
        if upry>len(z): upry=len(z)
        
        for x in range(lwrx,uprx):
            for y in range(lwry,upry):
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