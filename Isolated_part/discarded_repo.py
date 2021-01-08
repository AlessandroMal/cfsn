# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:09:37 2021

@author: alego
"""

def paramvsTip(surf, pxlen, thres, tipFunc, h, rtipmin, rtipmax, rtipstep):
    R = np.linspace(rtipmin, rtipmax, rtipstep)
    surfParams = calcParams(surf, pxlen, thres)
    imgParams = []
 
    tippos=str(tipFunc).find('Tip') #costuisco titolo file
    ind=tippos-1
    tipname=''
    while str(tipFunc)[ind]!=' ':
        ind-=1
    ind+=4
    while ind<tippos+3:
        tipname+=str(tipFunc)[ind]
        ind+=1
        
    out=open('paramvs'+tipname+'.dat', 'w')
    for rtip in R:
        out.write(str(rtip)+' ')
        tip = tipFunc(pxlen, h, r=rtip)
        img = mph.grey_dilation(surf, structure=-tip)
        imgParams.append(calcParams(img,pxlen,thres))
#        for i in imgParams[0]:
#            if surfParams[i]==0: imgParams[-1][i]=0
#            imgParams[-1][i]=imgParams[-1][i]/surfParams[i]
        print('progress: ' + str((np.where(R==rtip)[0][0] + 1)/len(R) * 100) + '%')
    
    out.write('Rtip (Surface par on last position)\n')
    for i in imgParams[0]:
        for el in imgParams:
            out.write(str(el[i])+' ')
        out.write(str(surfParams[i])+' ')
        out.write(i+'\n')
    out.close()
#    imgParams.append('image val/surface val')
#    Rlist = []
#    Rlist.append(R)
#    Rlist.append(r'$R_{tip}$')
#    return Rlist, imgParams
    return 'paramvs'+tipname+'.dat'


def paramvsNpart(Npx, pxlen, rmin, rmax, N_part_min, N_part_max, N_partstep,
                 tipFunc, h, parname, tipmin, tipmax, tipstep,
                 N_sample, paramFunc, y_label, err, relative=False):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    N_part = np.linspace(N_part_min, N_part_max, N_partstep)
    tip = np.linspace(tipmin, tipmax, tipstep)
    
    plt.figure()
    np.random.seed(123)
    plt_colors = [np.random.random(3) for _ in range(len(tip) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        z_param = []
        img_param = []
        xyr = []
        z=mf.genFlat(Npx)

        for N in N_part:
            print( 'N_part = ' + str(int(N)) )
            z, xyr = mf.genUnifIsolSph(z,pxlen,int(N),rmin,rmax, xyr, True)
      #      mf.plotfalsecol(z,pxlen)
            z_param.append(paramFunc(z))
        #    print('max height surface=',h_max(z,10))

            for m in tip:
                if parname=='angle': tip_ar = tipFunc(pxlen,h,angle=m)
                if parname=='r': tip_ar = tipFunc(pxlen,h,r=m)
                img_ar = mph.grey_dilation(z, structure=-tip_ar)
                
                if relative:
                    if z_param[-1]==0: img_param.append(0)
                    else: img_param.append(paramFunc(img_ar)/z_param[-1]) 
                else: img_param.append(paramFunc(img_ar))
         #       print('max height image=',h_max(z,10))
                
        if not(relative): #plotto anche valori di superficie
            plt_label = 'surface' if i==0 else '' # visualizza label solo una volta
            plt.errorbar(N_part/mf.Ncp(pxlen,Npx,(rmin+rmax)/2), z_param, yerr=err*np.array(z_param), color=plt_colors[-1], label=plt_label)
        for j in range(len(tip)):
            if parname=='r':
                plt_label = r'$R_{tip}/ \langle R_{part} \rangle=$'+str( round(2*tip[j]/(rmin+rmax), 2) )  if i==0 else '' # visualizza label solo una volta
            if parname=='angle':
                plt_label = r'$\alpha=$'+str( round(tip[j]*180/3.14159, 1) )+'°'  if i==0 else '' # visualizza label solo una volta
            
            plt.errorbar(N_part/mf.Ncp(pxlen,Npx,(rmin+rmax)/2), img_param[j::len(tip)], yerr=err*np.array(img_param[j::len(tip)]), color=plt_colors[j], label = plt_label)

    plt.title(r'res: $L_{map}=$'+str(Npx)+',  ' + r'$\langle R_{part} \rangle=$'+str(round( (rmin+rmax)/2/pxlen ,3))+',  '+r'$h_{tip}=$'+str(round(h/pxlen, 3)) )
    plt.xlabel(r'$N_{part}/N_{cp}$')
    plt.ylabel(y_label)
    plt.grid()
    plt.legend()
    plt.tight_layout()

def paramvsRpart(Npx, pxlen, N_part, R_part_min, R_part_max, R_part_step,
                 tipFunc, h, parname, tipmin, tipmax, tipstep,
                 N_sample, paramFunc, y_label, err, relative=False):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    
    R_part = np.linspace(R_part_min, R_part_max, R_part_step)
    tip = np.linspace(tipmin, tipmax, tipstep)
 #   zflat=mf.genFlat(Npx)
    
    plt.figure()
    np.random.seed(123)
    plt_colors = [np.random.random(3) for _ in range(len(tip) + 1)] # +1 per la superficie stessa
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        z_param = []
        img_param = []

        for rad in R_part:
            print( 'R_part = ' + str(rad), '   Ncp = '+str(mf.Ncp(pxlen,Npx,rad)) )
            z = mf.genUnifIsolSph(mf.genFlat(Npx),pxlen,N_part,rad,rad)
       #     mf.plotfalsecol(z,pxlen)
            z_param.append(paramFunc(z))
        #    print('max height surface=',h_max(z,10))

            for m in tip:
                if parname=='angle': tip_ar = tipFunc(pxlen,h,angle=m)
                if parname=='r': tip_ar = tipFunc(pxlen,h,r=m)
                img_ar = mph.grey_dilation(z, structure=-tip_ar)
                
                if relative:
                    if z_param[-1]==0: img_param.append(0)
                    else: img_param.append(paramFunc(img_ar)/z_param[-1]) 
                else: img_param.append(paramFunc(img_ar))
         #       print('max height image=',h_max(z,10))
                
        if not(relative): #plotto anche valori di superficie
            plt_label = 'surface' if i==0 else '' # visualizza label solo una volta
            plt.errorbar(R_part/pxlen, z_param, yerr=err*np.array(z_param), color=plt_colors[-1], label=plt_label)
        for j in range(len(tip)):
            if parname=='r':
                plt_label = r'$R_{tip}/L_{px}=$'+str( round(tip[j]/pxlen, 2) )  if i==0 else '' # visualizza label solo una volta
            if parname=='angle':
                plt_label = r'$\alpha=$'+str( round(tip[j]*180/3.14159, 1) )+'°'  if i==0 else '' # visualizza label solo una volta
            
            plt.errorbar(R_part/pxlen, img_param[j::len(tip)], yerr=err*np.array(img_param[j::len(tip)]), color=plt_colors[j], label = plt_label)

    plt.title(r'$N_{part}=$'+str(N_part)+',  '+ r'res: $L_{map}=$'+str(Npx)+',  '+ r'$h_{tip}=$'+str(round(h/pxlen, 3)) )
    plt.xlabel(r'$R_{part}/L_{px}$')
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
    
    
    
#MAP-------------------------------------------------

#z=mf.genFlat(Npx)
#z=mf.genHexSpikes(z,pxlen,20,2,20)
#z=mf.genNormNoise(z,pxlen,100,50)
#z=mf.genFlat(Npx)
#z=gaussian_filter(z,20) #non funziona?
#z=mf.genSphere(z,pxlen,np.array([100,100,20,70]),np.array([10,30]))
#z=mf.genUnifSph(z,pxlen,10,5,10) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genNormSph(z,pxlen,8,15,5) #, xmin=22, xmax=58, ymin=62, ymax=78)
#z=mf.genHexSpikes(z,pxlen,40,1,40)
#z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax) #, xmin=22, xmax=58, ymin=62, ymax=78)
#Tip------------------------------------------------
#occhio che h/a>>pxlen

#tip=mf.genParabolicTip(pxlen,50,r=8) 
#tip=mf.genPyramidTip(pxlen,50,angle=np.pi/4)
#tip=mf.genConeTip(pxlen,50,angle=np.pi/2)
#tip=mf.genSemisphTip(pxlen,50,r=20)

#IMG------------------------------------------------

#img = mph.grey_dilation(z, structure=-tip)

#PLOT-----------------------------------------------

#mf.plotview(z,pxlen,30,30)
#mf.plotview(z,pxlen,90,-90)
#mf.plotview(z,pxlen,0,90)
#mf.plotfalsecol(z,pxlen)
#mf.plotview(img,pxlen,30,30)
#mf.plotview(img,pxlen,90,-90)
#mf.plotview(img,pxlen,0,90)
#mf.plotfalsecol(img,pxlen)
#mf.plotview(tip,pxlen,30,30)
#mf.plotview(tip,pxlen,90,-90)
#mf.plotview(tip,pxlen,0,90)
#mf.plotfalsecol(tip,pxlen)

#PARAMS---------------------------------------------

#print('hmax ',par.h_max(z,10))
#print(mf.Ncp(pxlen,Npx,rmax))
#zParams = par.calcParams(z,pxlen,thres)
#print('\nSURFACE:\n', zParams)
#print(par.specArea(z,pxlen))
#print(time.perf_counter()-start)
#imgParams = par.calcParams(img,pxlen,thres)
#print('\nIMAGE:\n', imgParams)


#z=mf.genFlat(Npx)
#z=mf.genUnifIsolSph(z,pxlen,Npart,rmin,rmax)

'''
x,y=par.paramvsTip(z, pxlen, thres, mf.genParabolicTip, 50, 8, 30, 2)
x[0]=x[0]/rmax
x[1]=x[1]+r'$/ \langle R_{part} \rangle$'
par.plotParams(x,y)
'''

'''
file=par.paramvsTip(z, pxlen, thres, mf.genParabolicTip, 50, 8, 30, 2)
output=open(file, 'a')
output.write(str((rmax+rmin)/2)+' Rtip')
output.close()
'''

'''
par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 14,
                 mf.genPyramidTip, 50, 'angle', np.pi/8, np.pi/2, 6,
                 1, lambda surf: par.V(surf, pxlen), r'$rel. V$',0.04, True)
'''

'''
par.paramvsNpart(Npx, pxlen, rmin, rmax, 1, 40, 2,
                 mf.genPyramidTip, 50, 'angle', np.pi/8, np.pi/2, 2,
                 1, lambda surf: par.V(surf, pxlen), 'V.dat')
'''

'''
par.paramvsRpart(Npx, pxlen, 15, 10, 30, 2,
                 mf.genPyramidTip, 65, 'angle', np.pi/8, np.pi/2, 2,
                 1, lambda surf: par.V(surf, pxlen), 'V.dat')
'''

#obj = mf.identObj(img,thres)
#for i in obj: mf.plotfalsecol(i,pxlen)


def genHexSpikes(z, pxlen, hmin, hmax, dist, parmin, parmax, pbroken, **kwargs):
    par=kwargs.get('par','r')
    if par=='r': R=parmax
    if par=='angle': R=hmax*np.tan(parmax/2)
    xmin=kwargs.get('xmin',R)
    ymin=kwargs.get('ymin',R)
    xmax=kwargs.get('xmax',pxlen*len(z)-R)
    ymax=kwargs.get('ymax',pxlen*len(z)-R)

    bol=True
    while xmin<xmax:
        if bol: y0=ymin
        else: y0=ymin+dist/2
        while y0<ymax:
            
            h=uniform(hmin, hmax)
            if par=='r':
                R=uniform(parmin, parmax)
                m=h/R
            if par=='angle':
                R=h*np.tan(uniform(parmin, parmax)/2)
                m=h/R
            e=uniform(0,emax)
            r=R*np.sqrt(1-e**2)
            
            
            
            lwrx=int((xmin-R)/pxlen)-1  #ottimizza l'algoritmo, non ciclo su
            if lwrx<0: lwrx=0            #tutta la mappa, importante se ho
            uprx=int((xmin+R)/pxlen)+1  #alta risoluzione
            if uprx>len(z): uprx=len(z)
            lwry=int((y0-R)/pxlen)-1
            if lwry<0: lwry=0
            upry=int((y0+R)/pxlen)+1
            if upry>len(z): upry=len(z)
            
            if uniform(0,1)>pbroken: #spike non rotta
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                        if R**2 - (x*pxlen - xmin)**2 - (y*pxlen - y0)**2 > 0:
                            z[y,x]+=m*(R-np.sqrt( (x*pxlen - xmin)**2 + (y*pxlen - y0)**2))
            
            else: #spike rotta
                p=[]
                for i in range(3): #prendo 3 punti (centro in 0)
                    xrnd=xmin+R #rejection method for point in the circle
                    yrnd=y0+R
                    while (xrnd-xmin)**2+(yrnd-y0)**2>R**2:
                        xrnd=uniform(xmin-R,xmin+R)
                        yrnd=uniform(y0-R  ,y0+R)

                    p.append(xrnd)
                    p.append(yrnd)
                    p.append(m*(R-np.sqrt( (xrnd-xmin)**2 + (yrnd-y0)**2)))
                a=(p[4]-p[1])*(p[8]-p[2])-(p[5]-p[2])*(p[7]-p[1]) #build plane
                b=-(p[3]-p[0])*(p[8]-p[2])+(p[5]-p[2])*(p[6]-p[0])
                c=(p[3]-p[0])*(p[7]-p[1])-(p[4]-p[1])*(p[6]-p[0])
                d=-(a*p[0]+b*p[1]+c*p[2])
                p.clear()
                for x in range(lwrx,uprx):
                    for y in range(lwry,upry):
                       # if R**2 - (x*pxlen - xmin)**2 - (y*pxlen - y0)**2 > 0:
                        if ((x*pxlen-xmin)/
                          
                            if -(a*x*pxlen+b*y*pxlen+d)/c>0:
                                z[y,x]+=min(m*(R-np.sqrt( (x*pxlen - xmin)**2 + (y*pxlen - y0)**2)), -(a*x*pxlen+b*y*pxlen+d)/c)

            y0+=dist
        xmin+=dist*np.sqrt(3)/2
        bol=not(bol)
    return z