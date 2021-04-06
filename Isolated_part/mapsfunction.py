import numpy as np
from random import uniform
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.ndimage
import scipy.ndimage.morphology as mph

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

def genUnifSolidSph(z,pxlen,Nsph,rmin,rmax, **kwargs):
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

def genHexSpikes(z, pxlen, hmin, hmax, dist, parmin, parmax, emax, pbroken, **kwargs):
    par=kwargs.get('par','r')
    if par=='r': R=parmax
    if par=='angle': R=hmax*np.tan(parmax/2)
    xmin=kwargs.get('xmin',R)
    ymin=kwargs.get('ymin',R)
    xmax=kwargs.get('xmax',pxlen*len(z)-R)
    ymax=kwargs.get('ymax',pxlen*len(z)-R)

    counter=0
    bol=True
    while xmin<xmax:
        if bol: y0=ymin
        else: y0=ymin+dist/2
        while y0<ymax:
            counter+=1
            
            h=uniform(hmin, hmax)
            if par=='r': R=uniform(parmin, parmax)
            if par=='angle': R=h*np.tan(uniform(parmin, parmax)/2)
            e=uniform(0,emax)
            if uniform(0,1)>0.5:
                r=R*np.sqrt(1-e**2)
            else:
                r=R
                R=r*np.sqrt(1-e**2)
            coeff2=r/R #calcolo coefficiente per altezza spike con base ellittica
            
            lwrx=int((xmin-max(R,r))/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
            if lwrx<0: lwrx=0            #tutta la mappa, importante se ho
            uprx=int((xmin+max(R,r))/pxlen)+1  #alta risoluzione
            if uprx>len(z): uprx=len(z)
            lwry=int((y0-max(R,r))/pxlen)-1
            if lwry<0: lwry=0
            upry=int((y0+max(R,r))/pxlen)+1
            if upry>len(z): upry=len(z)
            
            if uniform(0,1)>pbroken: #spike non rotta
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                         if ((x*pxlen-xmin)/R)**2 + ((y*pxlen-y0)/r)**2 <1:
                            if x*pxlen-xmin!=0:
                                coeff1=(y*pxlen-y0)/(x*pxlen-xmin)
                                x_ell=xmin + r/np.sqrt(coeff1**2 + coeff2**2)
                                y_ell=y0 + coeff1*(x_ell-xmin)
                            else:
                                x_ell=xmin
                                y_ell=y0 + r
                                
                            z[y,x]+=h*(1 - np.sqrt( ((x*pxlen-xmin)**2+(y*pxlen-y0)**2)/((x_ell-xmin)**2+(y_ell-y0)**2)  ) )
            
            else: #spike rotta
                p=[]
                for i in range(3): #prendo 3 punti (centro in 0)
                    xrnd=xmin+R #rejection method for point in the ellipse
                    yrnd=y0+r
                    x_ell = xmin
                    y_ell = y0
                    
                    while (xrnd-xmin)**2+(yrnd-y0)**2> (x_ell-xmin)**2+(y_ell-y0)**2:
                        xrnd=uniform(xmin-R,xmin+R)
                        yrnd=uniform(y0-r  ,y0+r)
                        if xrnd-xmin!=0:
                            coeff1=(yrnd-y0)/(xrnd-xmin)
                            x_ell=xmin + r/np.sqrt(coeff1**2 + coeff2**2)
                            y_ell=y0 + coeff1*(x_ell-xmin)
                        else:
                            x_ell=xmin
                            y_ell=y0 + r

                    p.append(xrnd)
                    p.append(yrnd)
                    p.append(h*(1 - np.sqrt( ((xrnd-xmin)**2+(yrnd-y0)**2)/((x_ell-xmin)**2+(y_ell-y0)**2)  ) ))
                a=(p[4]-p[1])*(p[8]-p[2])-(p[5]-p[2])*(p[7]-p[1]) #build plane
                b=-(p[3]-p[0])*(p[8]-p[2])+(p[5]-p[2])*(p[6]-p[0])
                c=(p[3]-p[0])*(p[7]-p[1])-(p[4]-p[1])*(p[6]-p[0])
                d=-(a*p[0]+b*p[1]+c*p[2])
                p.clear()
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if ((x*pxlen-xmin)/R)**2 + ((y*pxlen-y0)/r)**2 <1:
                            if x*pxlen-xmin!=0:
                                coeff1=(y*pxlen-y0)/(x*pxlen-xmin)
                                if x*pxlen>xmin: x_ell=xmin + r/np.sqrt(coeff1**2 + coeff2**2)
                                else: x_ell=xmin - r/np.sqrt(coeff1**2 + coeff2**2)
                                y_ell=y0 + coeff1*(x_ell-xmin)
                            else:
                                x_ell=xmin
                                y_ell=y0 + r
                            
                            if -(a*x*pxlen+b*y*pxlen+d)/c>0:
                                z[y,x]+=min(h*(1 - np.sqrt( ((x*pxlen-xmin)**2+(y*pxlen-y0)**2)/((x_ell-xmin)**2+(y_ell-y0)**2)  ) ), -(a*x*pxlen+b*y*pxlen+d)/c)

            y0+=dist
        xmin+=dist*np.sqrt(3)/2
        bol=not(bol)
    print(counter, 'particles generated')
    return z

def Ncp(pxlen,Npx,r):
    extraline=0 #calcolo numero di sfere per close packing
    d=2*r
    lin= int( pxlen*Npx/d +1)*2
    col= int( pxlen*Npx/(2*d*np.cos(np.pi/6)) ) +1
    if pxlen*Npx%d < d/2: lin=lin-1 #tolgo l'estremo se non c'è spazio
    if pxlen*Npx%(2*d*np.cos(np.pi/6) ) < d*np.cos(np.pi/6): extraline=int(lin/2)
    return lin*col-extraline

def genUnifIsolSph(z,pxlen,Nsph,rmin,rmax, xyr=[], centres=False, **kwargs):
    xmin=kwargs.get('xmin',0)
    ymin=kwargs.get('ymin',0)
    xmax=kwargs.get('xmax',pxlen*len(z))
    ymax=kwargs.get('ymax',pxlen*len(z))
    
    ncp=Ncp(pxlen,len(z), (rmin+rmax)/2)
    if Nsph>=ncp:
        print("Error: too many particles for this map")
        return 0

    n=len(xyr)/3
    while n<Nsph:
        if n<0.5*ncp:
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
                
                lwrx=int((x0-R)/pxlen)-2  #ottimizza l'algoritmo, non ciclo su
                if lwrx<0: lwrx=0         #tutta la mappa, importante se ho
                uprx=int((x0+R)/pxlen)+2  #alta risoluzione
                if uprx>len(z): uprx=len(z)
                lwry=int((y0-R)/pxlen)-2
                if lwry<0: lwry=0
                upry=int((y0+R)/pxlen)+2
                if upry>len(z): upry=len(z)
            
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                            z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R

        else:
            remaining=[]
            for x in range(len(z)):
                for y in range(len(z)): #ciclo sui punti della mappa
                    b=True
                    for i in range(int(len(xyr)/3)): #ciclo su sfere già messe
                        if rmax+xyr[3*i+2]>np.sqrt((x*pxlen-xyr[3*i])**2+(y*pxlen-xyr[3*i+1])**2):
                            b=False
                            break
                    if b:
                        remaining.append(x)
                        remaining.append(y)
                        
          #  print(len(remaining),n)
            if len(remaining)>0:
                n+=1
                rndp = int(uniform(0,len(remaining)/2))
                x0 = remaining[rndp*2]*pxlen
                y0 = remaining[rndp*2+1]*pxlen
                remaining.clear()
                R  = uniform(rmin, rmax)
                xyr.append(x0)
                xyr.append(y0)
                xyr.append(R)
                
                lwrx=int((x0-R)/pxlen)-2  #ottimizza l'algoritmo, non ciclo su
                if lwrx<0: lwrx=0         #tutta la mappa, importante se ho
                uprx=int((x0+R)/pxlen)+2  #alta risoluzione
                if uprx>len(z): uprx=len(z)
                lwry=int((y0-R)/pxlen)-2
                if lwry<0: lwry=0
                upry=int((y0+R)/pxlen)+2
                if upry>len(z): upry=len(z)
                
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2 > 0:
                            z[y,x]+=np.sqrt(R**2 - (x*pxlen - x0)**2 - (y*pxlen - y0)**2) + R
            else: break

    if n!=Nsph: print("space over: ",n," particles generated")
    if centres: return z, xyr
    else:
        xyr.clear()
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
 #   if px<10: print('Warning: low tip resolution')
    z=np.zeros([px,px])
    for x in range(len(z)):
        for y in range(len(z)):
            z[y,x]=a*((x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2)

    z = z -np.amin(z)
    return z
        
def genPyramidTip(pxlen,h, **kwargs):
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
        R=h/ar/2
        if h<R: l=2*np.sqrt(R**2 - (R-h)**2)
        else: l=h/ar
    else:
        R=r
        if h<R: l=2*np.sqrt(R**2 - (R-h)**2)
        else: l=2*r
    
    px=int(l/pxlen)
    if R*2/pxlen<10: print('Warning: low tip resolution')
    if h<R: z=np.ones([px,px])*R
    else: z=np.ones([px,px])*h
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
        l=2*h*np.tan(angle/2)
        m=2*h/l
    
    px=int(l/pxlen)
    z=np.zeros([px,px])
    for x in range(len(z)):
        for y in range(len(z)):
            z[y,x]+=m*np.sqrt( (x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2)

    z = z -np.amin(z)
    return z

def genRndSemisphTip(pxlen,h,low,up):
    R=uniform(low, up)
    print('generated rnd Rtip:',R)
    if h<R:
        l=2*np.sqrt(R**2 - (R-h)**2)
        px=int(l/pxlen)
        z=np.ones([px,px])*R
    else:
        l=2*R
        px=int(l/pxlen)
        z=np.ones([px,px])*h
    if R*2/pxlen<10: print('Warning: low tip resolution')
    print('npx tip:',px)
    for x in range(len(z)):
        for y in range(len(z)):
            r2= (x*pxlen - px/2*pxlen)**2 + (y*pxlen - px/2*pxlen)**2
            if r2 < R**2:
                z[y,x]=R*(1-np.sqrt(1-r2/R**2))

    z = z -np.amin(z)
    return z

'''
def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    obj_ind = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj_list = []
    
    for i in range(len(obj_ind)):
        z_single_obj = z.copy()
        z_single_obj[np.where(z_labeled!=i+1)] = 0
        obj_list.append(z_single_obj[obj_ind[i]])
    print('identObj found ' + str(len(obj_list)) + ' objects')
    return obj_list, z_labeled, obj_ind


def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    obj_ind = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj_list = []
    
    # prevent overlap of single structures
    for i in range(len(obj_ind)):
        z_single_obj = z.copy()
        z_single_obj[np.where(z_labeled!=i+1)] = 0
        
        # remove zero rows and columns
        z_single_obj = z_single_obj[~np.all(z_single_obj == 0, axis=1)]
        z_single_obj = z_single_obj[:, ~np.all(z_single_obj == 0, axis=0)]
        
        obj_list.append(z_single_obj)

    return obj_list, z_labeled, obj_ind
'''

def identObj(z, thres, Npx_min):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    obj_ind = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
        
    obj_list = []
    keep_obj_nr = []
    # prevent overlap of single structures
    for i in range(len(obj_ind)):
        z_single_obj = z.copy()
        z_single_obj[np.where(z_labeled!=i+1)] = 0
        
        # remove zero rows and columns
        z_single_obj = z_single_obj[~np.all(z_single_obj == 0, axis=1)]
        z_single_obj = z_single_obj[:, ~np.all(z_single_obj == 0, axis=0)]
        
        # keep only objects with a minimum pixel size
        if sum(z_single_obj.shape) >= Npx_min:
            obj_list.append(z_single_obj)
            keep_obj_nr.append(i)
    
    # relabel objects considering only kept objects
    keep_obj_nr = np.array(keep_obj_nr)
    i = 0
    while i < np.max(z_labeled):
        if i in keep_obj_nr:
            i += 1
            continue
        z_labeled[np.where(z_labeled==i+1)] = 0
        z_labeled[np.where(z_labeled>i+1)] += -1
        keep_obj_nr += -1
        i += 1
        
    obj_ind = [obj_ind[i] for i in keep_obj_nr]
        
    print('identObj found ' + str(len(obj_list)) + ' objects')   

    return obj_list, z_labeled, obj_ind

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
    plt.show()
    
def plotfalsecol(z,pxlen, figx=0, figy=0):
    if figx!=0: plt.figure(figsize=(figx,figy))
    else: plt.figure()
    plt.title(r'$ L_{px}= '+str(pxlen)+'nm $')
    plt.axis('equal')
    
    X = np.linspace(0,pxlen*np.shape(z)[1],num=np.shape(z)[1])
    Y = np.linspace(0,pxlen*np.shape(z)[0],num=np.shape(z)[0])
    X, Y = np.meshgrid(X, Y)
    
    plt.pcolormesh(X,Y,z)
    clb = plt.colorbar()    
    
    clb.set_label('Z (nm)')
    plt.xlabel('X (nm)')
    plt.ylabel('Y (nm)')
    plt.show()
    
    
def plotThres(z, z_labeled, pxlen):
    plotfalsecol(z, pxlen)
    for i in range(1, int(np.max(z_labeled)) + 1):
        obj_i_edge = (z_labeled==i) & mph.binary_dilation(z_labeled!=i, structure=np.ones([3,3])) # edge is part of structure
        index_x = np.where(obj_i_edge==1)[1] + 0.5
        index_y = np.where(obj_i_edge==1)[0] + 0.5
        plt.plot(index_x*pxlen, index_y*pxlen, color='r', marker='s', linestyle='None', markersize=2)
    plt.tight_layout()