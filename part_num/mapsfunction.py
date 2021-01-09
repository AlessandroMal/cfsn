from operator import index
import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from matplotlib import cm, markers
import scipy.ndimage.morphology as mph
import scipy.ndimage
from scipy.stats import lognorm

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

def genLogNormSolidSph(z,pxlen,Nsph,av,var, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
        
    for i in range(Nsph):
        R = np.random.lognormal(np.log(av / np.sqrt(1 + var**2 / av**2)), np.sqrt(np.log(1 + (var/av)**2)))
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
        
        xtop=int((lwrx+uprx)/2)
        ytop=int((lwry+upry)/2)
        ztop=z[ytop,xtop]
        for x in range(lwrx+1,uprx): #mi appoggio sul punto più alto nel quadrato
            for y in range(lwry+1,upry): #in cui poggia la sfera
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0 and z[y,x] > ztop:
                    ztop=z[y,x]
                    xtop=x
                    ytop=y
        
        r_top=np.sqrt((xtop*pxlen-x0)**2+(ytop*pxlen-y0)**2)
        z0=ztop+np.sqrt(R**2 - r_top**2) #pongo il centro della sfera
        
        lwrxcap=int((x0-r_top)/pxlen)-1  #vedo se ci sono appoggi migliori sotto
        if lwrxcap<0: lwrxcap=0
        uprxcap=int((x0+r_top)/pxlen)+1
        if uprxcap>len(z): uprxcap=len(z)
        lwrycap=int((y0-r_top)/pxlen)-1
        if lwrycap<0: lwrycap=0
        uprycap=int((y0+r_top)/pxlen)+1
        if uprycap>len(z): uprycap=len(z)        
        
        r_top=R
        for x in range(lwrxcap,uprxcap): #trovo contect point
            for y in range(lwrycap,uprycap):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0 and np.sqrt( (z[y,x]-z0)**2+(x*pxlen-x0)**2+(y*pxlen-y0)**2) <r_top: 
                    xtop=x
                    ytop=y
                    r_top=np.sqrt( (z[y,x]-z0)**2+(x*pxlen-x0)**2+(y*pxlen-y0)**2)
        
        r_top=np.sqrt((xtop*pxlen-x0)**2+(ytop*pxlen-y0)**2)
        z0=z[ytop,xtop]- (R - np.sqrt(R**2 - r_top**2)) #punto piu basso della nuova sfera
    
        for x in range(lwrx,uprx): #finalemente deposito
            for y in range(lwry,upry):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x]=z0 + np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

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

def genLogNormSph(z,pxlen,Nsph,av,var, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))

    for i in range(Nsph):

        R = np.random.lognormal(np.log(av / np.sqrt(1 + var**2 / av**2)), np.sqrt(np.log(1 + (var/av)**2)))
        
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

def genIsolLogNormSph(z,pxlen,Nsph,av,var, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))

    i = 0
    while i < Nsph:

        R = np.random.lognormal(np.log(av / np.sqrt(1 + var**2 / av**2)), np.sqrt(np.log(1 + (var/av)**2)))

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
        
        valid_position = True
        for x in range(lwrx,uprx):
            for y in range(lwry,upry):
                if z[y,x]!=0: 
                    valid_position = False
        
        if valid_position==False:
            continue
        
        for x in range(lwrx,uprx):
            for y in range(lwry,upry):
                if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                    z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R
        
        i+=1
        
    return z

def genParabolicTip(pxlen,h, **kwargs):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    ar=kwargs.get('ar',-1)
    r=kwargs.get('r',-1)
    curv=kwargs.get('curv',-1)
    
    if ar*r*curv<0:
        print('Error: tip parameters conflict: priority ar>r>curv')
    if ar>0 and r>0 and curv>0:
        print('Error: tip parameters conflict: used ar')
    
    if ar>0:
        l=h/ar
        a=4*ar/l
    else:
        if curv>0:
            l=2*np.sqrt(2*h/curv)
            a=curv/2
        if r>0:
            l=2*np.sqrt(2*r*h)
            a=1/2/r
    
    px=int(l/pxlen)
    z=np.zeros([px,px])
    for x in range(len(z)):
        for y in range(len(z)):
            if (l/2)**2 - (x*pxlen - px/2*pxlen)**2 - (y*pxlen - px/2*pxlen)**2 > 0:
                z[y,x]=a*((x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2)
            else:
                z[y,x]=h

    z = z -np.amin(z)
    return z
        
def genPyramidTip(pxlen,h, **kwargs):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    ar=kwargs.get('ar',-1)
    sideangle=kwargs.get('angle',-1)
    r=kwargs.get('r',-1)
    
    if ar*sideangle*r<0:
        print('Error: tip parameters conflict: priority ar>sideangle>sideL')
    if ar>0 and sideangle>0 and r>0:
        print('Error: tip parameters conflict: used ar')
    
    if ar>0:
        l=h/ar
        m=2*ar
    else:
        if r>0:
            l=r*2
        if sideangle>0:
            l=2*h*np.tan(sideangle/2)
        m=2*h/l
    
    px=int(l/pxlen)
    z=np.zeros([px,px])
    for x in range(len(z)):
        for y in range(len(z)):
                z[y,x]=m*max(abs(x*pxlen-px/2*pxlen), abs(y*pxlen-px/2*pxlen))

    z = z -np.amin(z)
    return z

def genSemisphTip(pxlen,h, **kwargs):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    ar=kwargs.get('ar',-1)
    r=kwargs.get('r',-1)
    
    if ar*r>0:
        print('Error: tip parameters conflict: used ar')
    
    if ar>0:
        l=h/ar
        R=l/2
    else:
        l=2*r
        R=r
    
    px=int(l/pxlen)
    z=np.ones([px,px])*h
    for x in range(len(z)):
        for y in range(len(z)):
            r2= (x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2
            if r2 < R**2:
                z[y,x]=R*(1-np.sqrt(1-r2/R**2))

    z = z -np.amin(z)
    return z

def genConeTip(pxlen,h, **kwargs):
    #OCCHIO: PER AVERE UNA BUONA FIGURA MI SERVE l>>pxlen
    r=kwargs.get('r',-1)
    angle=kwargs.get('angle',-1)
    
    if r*angle>0:
        print('Error: tip parameters conflict: used ar')
    
    if r>0:
        l=r*2
        m=h/r
    else:
        l=2*h/np.tan(angle/2)
        m=2*h/l
    
    px=int(l/pxlen)
    z=np.zeros([px,px])
    for x in range(len(z)):
        for y in range(len(z)):
            if (l/2)**2 - (x*pxlen - px/2*pxlen)**2 - (y*pxlen - px/2*pxlen)**2 > 0:
                z[y,x]+=m*np.sqrt( (x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2)
            else:
                z[y,x]=h

    z = z -np.amin(z)
    return z

def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    obj_ind = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj_list = [z[i] for i in obj_ind]    

    return obj_list, z_labeled

#def optThres(z, prec):
#    thres = 0
#    while :
#        obj_list, z_thres = identObj(z, thres)
#        N_obj = len(obj_list)
#        z_mean = np.mean(z_thres)

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

def plotThres(z, z_labeled, pxlen, title):
    plotfalsecol(z, pxlen)
    for i in  range(1, np.max(z_labeled) + 1):
        obj_i_edge = (z_labeled==i) & mph.binary_dilation(z_labeled!=i, structure=np.ones([3,3])) # edge is part of structure
        index_x = np.where(obj_i_edge==1)[1]
        index_y = np.where(obj_i_edge==1)[0]
        plt.scatter(index_x*pxlen + pxlen/2, index_y*pxlen + pxlen/2, color='r', s=0.25)
        plt.title(title)