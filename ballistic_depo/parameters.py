import numpy as np
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
import mapsfunction as mf
import skimage.measure as sk
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
#from numba.typed import List
#from numba import guvectorize, float64

def skew(z):
    std=np.std(z)
    if std!=0: z_skew = np.mean((z - np.mean(z))**3)/ std**3
    else: z_skew=0
    return z_skew

def V(z,pxlen): return pxlen**len(np.shape(z)) *np.sum(z)

def specArea(z):
#first order specific area
    A=0
    for x in range(np.shape(z)[1]):
        for y in range(np.shape(z)[0]):
            A+=np.linalg.norm(np.array([-z[y,x+1]+z[y,x], -z[y+1,x]+z[y,x], 1]))/2
            A+=np.linalg.norm(np.array([z[y+1,x]-z[y+1,x+1], z[y,x+1]-z[y+1,x+1], 1]))/2
    return A/( (np.shape(z)[1]-1)*(np.shape(z)[0]-1) )

def coverage(z, thres):
    N_tot = np.shape(z)[1]*np.shape(z)[0]
    N = np.sum(z>thres)
    cov = N/N_tot
    return cov

def h_max(z,n):
    top=[]
    for x in range(np.shape(z)[1]):
        for y in range(np.shape(z)[0]):
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

#@njit
#@guvectorize([(float64[:,:], float64, float64[:], float64[:])], '(n,n),()->(n),(n)')
def C_fast(z, indz, bin_size):
    print('computing height correlation function')
    r=np.arange(0, 2/3* min(np.shape(z)[1], np.shape(z)[0])+bin_size, bin_size)
    counter=np.zeros(len(r)-1)
    hh=np.zeros(len(r)-1)
    
    wh=np.where(z>=np.amin(z))
#    np.savetxt('wh.dat', wh[0][0::np.shape(z)[1]])

    for x0 in range(np.shape(z)[1]):
        if (x0+1)%int(np.shape(z)[1]/10)==0: print(x0/np.shape(z)[1]*100,'%',end=' ')
        for y0 in range(np.shape(z)[0]):
            for i in range(len(r)-1):
                hh[i]=np.sum(z[indz])*z[y0,x0]
                counter[i]+=z[indz].size
    
#    rout=np.delete(r,-1)+bin_size
#    hout/=counter
#    return np.delete(r, -1) + bin_size, hh/counter

    
def G_fast(z, bin_size):
    print('computing height-height correlation function')
    r=np.arange(0, 2/3* min(np.shape(z)[1], np.shape(z)[0])+bin_size, bin_size)
    counter=np.zeros(len(r)-1)
    hh=np.zeros(len(r)-1)
    
    wh=np.where(z>=np.amin(z))
#    np.savetxt('wh.dat', wh[0][0::np.shape(z)[1]])

    for x0 in range(np.shape(z)[1]):
        if (x0+1)%int(np.shape(z)[1]/10)==0 and x0+1!=np.shape(z)[1]: print(round((x0+1)/np.shape(z)[1]*100, 0),'%',end=' ')
        for y0 in range(np.shape(z)[0]):
            for i in range(len(r)-1):
#                index=(wh[max(0,y0-int(r[i+1])-1):min(np.shape(z)[0],y0+int(r[i+1])+2), max(0,x0-int(r[i+1])-1):min(np.shape(z)[1],x0+int(r[i+1])+2)] <r[i+1])
                index= (np.sqrt((wh[0][0:] - x0)**2 + (wh[1][0:] - y0)**2) >= r[i]) & (np.sqrt((wh[0][0:] - y0)**2 + (wh[1][0:] - x0)**2) < r[i+1])
                index= index.reshape([np.shape(z)[0],np.shape(z)[1]])
                
                hh[i]=np.sum((z[y0,x0]-z[index])**2)
                counter[i]+=z[index].size
    
#    rout=np.delete(r,-1)+bin_size
#    hout/=counter
    return np.delete(r, -1) + bin_size, hh/counter

def C(z, bin_size):
    print('computing height correlation function')
    r=np.arange(0, 2/3* min(np.shape(z)[1], np.shape(z)[0])+bin_size, bin_size)
    counter=np.zeros(len(r)-1)
    hh=np.zeros(len(r)-1)
    
    wh=np.where(z>=np.amin(z))
#    np.savetxt('wh.dat', wh[0][0::np.shape(z)[1]])

    for x0 in range(np.shape(z)[1]):
        if (x0+1)%int(np.shape(z)[1]/10)==0: print(x0/np.shape(z)[1]*100,'%',end=' ')
        for y0 in range(np.shape(z)[0]):
            for i in range(len(r)-1):
#                index=(wh[max(0,y0-int(r[i+1])-1):min(np.shape(z)[0],y0+int(r[i+1])+2), max(0,x0-int(r[i+1])-1):min(np.shape(z)[1],x0+int(r[i+1])+2)] <r[i+1])
                index= (np.sqrt((wh[0][0:] - x0)**2 + (wh[1][0:] - y0)**2) >= r[i]) & (np.sqrt((wh[0][0:] - y0)**2 + (wh[1][0:] - x0)**2) < r[i+1])
                index= index.reshape([np.shape(z)[0],np.shape(z)[1]])
                
                hh[i]+=np.sum(z[y0,x0]*z[index])
                counter[i]+=z[index].size
    
#    rout=np.delete(r,-1)+bin_size
#    hout/=counter
    return np.delete(r, -1) + bin_size, hh/counter

    
def G(z, bin_size):
    print('computing height-height correlation function')
    r=np.arange(0, 2/3* min(np.shape(z)[1], np.shape(z)[0])+bin_size, bin_size)
    counter=np.zeros(len(r)-1)
    hh=np.zeros(len(r)-1)
    
    wh=np.where(z>=np.amin(z))
#    np.savetxt('wh.dat', wh[0][0::np.shape(z)[1]])

    for x0 in range(np.shape(z)[1]):
        if (x0+1)%int(np.shape(z)[1]/10)==0 and x0+1!=np.shape(z)[1]: print(round((x0+1)/np.shape(z)[1]*100, 0),'%',end=' ')
        for y0 in range(np.shape(z)[0]):
            for i in range(len(r)-1):
#                index=(wh[max(0,y0-int(r[i+1])-1):min(np.shape(z)[0],y0+int(r[i+1])+2), max(0,x0-int(r[i+1])-1):min(np.shape(z)[1],x0+int(r[i+1])+2)] <r[i+1])
                index= (np.sqrt((wh[0][0:] - x0)**2 + (wh[1][0:] - y0)**2) >= r[i]) & (np.sqrt((wh[0][0:] - y0)**2 + (wh[1][0:] - x0)**2) < r[i+1])
                index= index.reshape([np.shape(z)[0],np.shape(z)[1]])
                
                hh[i]+=np.sum((z[y0,x0]-z[index])**2)
                counter[i]+=z[index].size
    
#    rout=np.delete(r,-1)+bin_size
#    hout/=counter
    return np.delete(r, -1) + bin_size, hh/counter

def C_profile(z, pxlen, bin_size):
    L_x=np.shape(z)[1]
    r=np.arange(bin_size/2, L_x*pxlen, bin_size)
    counter=np.zeros(len(r))
    hh=np.zeros(len(r))
    av=np.mean(z)
    
    for x0 in range(0, L_x):
        for x in range(x0,L_x):
            bin_hist=int(abs(x0-x)*pxlen/bin_size)
            if bin_hist<len(r):
                counter[bin_hist]+=np.shape(z)[0]
                hh[bin_hist]+=np.sum( (z[0:,x0]-av)*(z[0:,x]-av) )
                
    for y0 in range(0, L_x):
        for y in range(y0,L_x):
            bin_hist=int(abs(y0-y)*pxlen/bin_size)
            if bin_hist<len(r):
                counter[bin_hist]+=np.shape(z)[1]
                hh[bin_hist]+=np.sum( (z[y0,0:]-av)*(z[y,0:]-av) )
    
    return r-0.5, hh/counter

def G_profile(z, pxlen, bin_size):
    L_x=np.shape(z)[1]
    r=np.arange(bin_size/2, L_x*pxlen, bin_size)
    counter=np.zeros(len(r))
    hh=np.zeros(len(r))
    
    for x0 in range(0, L_x):
        for x in range(x0,L_x):
            bin_hist=int(abs(x0-x)*pxlen/bin_size)
            if bin_hist<len(r):
                counter[bin_hist]+=np.shape(z)[0]
                hh[bin_hist]+=np.sum((z[0:,x0]-z[0:,x])**2)
                
    for y0 in range(0, L_x):
        for y in range(y0,L_x):
            bin_hist=int(abs(y0-y)*pxlen/bin_size)
            if bin_hist<len(r):
                counter[bin_hist]+=np.shape(z)[1]
                hh[bin_hist]+=np.sum((z[y0,0:]-z[y,0:])**2)
    
    return r-0.5, hh/counter

def C_1d(z):
    Corr=np.array([])
    av=np.mean(z)
    
    for x0 in range(len(z)):
        s=0
        for x in range(len(z)-x0):
            s+=(z[x]-av)*(z[x+x0]-av)
        Corr=np.append(Corr, s/(len(z)-x0))

    return Corr

def G_1d(z):
    Corr=np.array([])
    
    for x0 in range(len(z)):
        s=0
        for x in range(len(z)-x0):
            s+=(z[x]-z[x+x0])**2
        Corr=np.append(Corr, s/(len(z)-x0))

    return Corr

def wavelength(z,pxlen,direction):
    if direction=='x': 
        av=np.zeros(np.shape(z)[0])
        for y in range(len(av)):
            av[y]=(z[y,0]-z[y,-1])/(np.shape(z)[1]*pxlen)
        
        phi_av=np.mean(av)
        sm=0
        for y in range(np.shape(z)[0]):
            for x in range(np.shape(z)[1]-1):
                sm+=( (z[y,x+1]-z[y,x])/pxlen - phi_av )**2
                
    if direction=='y': 
        av=np.zeros(np.shape(z)[1])
        for x in range(len(av)):
            av[x]=(z[0,x]-z[-1,x])/(np.shape(z)[0]*pxlen)
        
        phi_av=np.mean(av)
        sm=0
        for x in range(np.shape(z)[1]):
            for y in range(np.shape(z)[0]-1):
                sm+=( (z[y+1,x]-z[y,x])/pxlen - phi_av )**2
        
    deltaq=np.sqrt(sm/(np.shape(z)[1]-1)/(np.shape(z)[0]-1))
    return 2*np.pi*np.std(z)/deltaq

def calcParams(z,pxlen,thres):
    param_list = [h_max(z,10),
                  np.mean(z),
                  np.std(z),
                  skew(z),
                  V(z,pxlen),
                  specArea(z,pxlen),
                  coverage(z,thres)]
    return param_list

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
        print('progress: ' + str((np.where(R==rtip)[0][0] + 1)/len(R) * 100) + '%')
    
    out.write('Rtip (Surface par on last position)\n')
    for i in imgParams[0]:
        for el in imgParams:
            out.write(str(el[i])+' ')
        out.write(str(surfParams[i])+' ')
        out.write(i+'\n')
    out.close()
    return 'paramvs'+tipname+'.dat'

def paramvsNpart(Npx, pxlen, rmin, rmax, N_part_min, N_part_max, N_partstep,
                 tipFunc, h, tipparname, tipmin, tipmax, tipstep,
                 N_sample, paramFunc, filename):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    N_part = np.linspace(N_part_min, N_part_max, N_partstep)
    tip = np.linspace(tipmin, tipmax, tipstep)
    
    out=open(filename, 'w')
    out.write(str(Npx)+' '+str(pxlen)+' '+str((rmax+rmin)/2)+' '+str(h)+' '+str(mf.Ncp(pxlen,Npx,(rmin+rmax)/2))+' ' )
    out.write('#Npx, pxlen, avR_part, h_tip, Ncp(avR_part)\n')
    for m in tip:
        out.write(str(m)+' ')
    out.write('#tip_'+tipparname+' (rows)\n')
        
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        z_param = []
        img_param = []
        xyr = []
        z=mf.genFlat(Npx)

        for N in N_part:
            out.write(str(int(N))+' ')
            print( 'N_part = ' + str(int(N)) )
            z, xyr = mf.genUnifIsolSph(z,pxlen,int(N),rmin,rmax, xyr, True)
      #      mf.plotfalsecol(z,pxlen)
            z_param.append(paramFunc(z))
        #    print('max height surface=',h_max(z,10))

            for m in tip:
                if tipparname=='angle': tip_ar = tipFunc(pxlen,h,angle=m)
                if tipparname=='r': tip_ar = tipFunc(pxlen,h,r=m)
                img_ar = mph.grey_dilation(z, structure=-tip_ar)
                img_param.append(paramFunc(img_ar))
         #       print('max height image=',h_max(z,10))
         
        out.write('#Npart\n')
        for j in range(len(N_part)):
            out.write(str(z_param[j])+' ')
        out.write('#Surface par\n')
        
        for i in range(len(tip)):
            for j in range(len(N_part)):
                out.write(str(img_param[i+ j*len(tip)])+' ')
            out.write('\n')
        out.write('\n')
    out.close()
    print('datas printed in '+filename)

def paramvsRpart(Npx, pxlen, N_part, R_part_min, R_part_max, R_part_step,
                 tipFunc, h, tipparname, tipmin, tipmax, tipstep,
                 N_sample, paramFunc, filename):
    '''tipFunc deve essere una funzione.
    paramFunc deve essere una funzione sola della superficie/ dell'immagine.
    '''
    R_part = np.linspace(R_part_min, R_part_max, R_part_step)
    tip = np.linspace(tipmin, tipmax, tipstep)
    
    out=open(filename, 'w')
    out.write(str(Npx)+' '+str(pxlen)+' '+str(N_part)+' '+str(h)+' ' )
    out.write('#Npx, pxlen, N_part, h_tip\n')
    for m in tip:
        out.write(str(m)+' ')
    out.write('#tip_'+tipparname+' (rows)\n')
    
    for i in range(N_sample):
        print('N_sample = ' + str(i + 1))
        z_param = []
        img_param = []

        for rad in R_part:
            out.write(str(rad)+' ')
            print( 'R_part = ' + str(rad), '   Ncp = '+str(mf.Ncp(pxlen,Npx,rad)) )
            z = mf.genUnifIsolSph(mf.genFlat(Npx),pxlen,N_part,rad,rad)
       #     mf.plotfalsecol(z,pxlen)
            z_param.append(paramFunc(z))
        #    print('max height surface=',h_max(z,10))

            for m in tip:
                if tipparname=='angle': tip_ar = tipFunc(pxlen,h,angle=m)
                if tipparname=='r': tip_ar = tipFunc(pxlen,h,r=m)
                img_ar = mph.grey_dilation(z, structure=-tip_ar)
                img_param.append(paramFunc(img_ar))
         #       print('max height image=',h_max(z,10))
        
        out.write('#Rpart\n')
        for j in range(len(R_part)):
            out.write(str(z_param[j])+' ')
        out.write('#Surface par\n')
        
        for i in range(len(tip)):
            for j in range(len(R_part)):
                out.write(str(img_param[i+ j*len(tip)])+' ')
            out.write('\n')
        out.write('\n')
    out.close()
    print('datas printed in '+filename)

def capPar(z,pxlen,thres):
    h=np.amax(z) #già in unità fisiche l'asse z, come anche V
    props=sk.regionprops( sk.label(z>thres) )
    e=props[0].eccentricity #eccentricità
    a=props[0].equivalent_diameter /2 #raggio equivalente cap
    return h, a*pxlen, props[0].area *pxlen**2, V(z,pxlen), e, sk.label(z>thres)

def revimg_filter(obj,pxlen,thres, etol, msqertol, relVtol):
    filtered=[]
    listvol=[]
    h=[]
    for i in range(len(obj)): #filtro con eccentricità
        listvol.append(capPar(obj[i],pxlen,thres)[3])
        h.append(capPar(obj[i],pxlen,thres)[0])
        if capPar(obj[i],pxlen,thres)[4] < etol: filtered.append(i)
    
    print('eccentricity_filter ha tenuto',len(filtered),'particelle')
    logv= np.log10(np.array(listvol))
    logv= logv.reshape((-1, 1))
    h = np.array(h)
    
    model  = LinearRegression().fit(logv, h)
    h_pred = model.predict(logv)
    msqer  = mean_squared_error(h, h_pred)
    
    regrfiltered=[]
    for i in filtered: #filtro con regressione h(logV) lineare
        if abs(h_pred[i] - h[i]) < msqertol*msqer: regrfiltered.append(i)
        
    print('linear_regression_filter ha tenuto',len(regrfiltered),'particelle')
    filtered.clear()
    for i in regrfiltered: #filtro per broken spike: errore relativo Vcap
        if 1 - 6*listvol[i]/(np.pi*h[i]**3 + 3*h[i]*capPar(obj[i],pxlen,thres)[2]) < relVtol: filtered.append(i)
        
    filtered_obj=[]
    for i in filtered:
        filtered_obj.append(obj[i])
    print('relativeVerr_filter ha tenuto',len(filtered_obj),'particelle')
    return filtered_obj
    
def singlePartAnalysis(surf, img, pxlen, thres, h_tip, ar_tip, N_part, N_pxl, 
                        param_string_list, binwidth_list, bin_min_list, bin_max_list):
    
    # get single particles 
    surf_obj_list, surf_labeled = mf.identObj(surf, thres)
    img_obj_list, img_labeled = mf.identObj(img, thres)
    mf.plotThres(surf, surf_labeled, pxlen, 'surface (found ' + str(len(surf_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm)')
    mf.plotThres(img, img_labeled, pxlen, 'image (found ' + str(len(img_obj_list)) + ' of ' + str(N_part) +' particles, thres=' + str(thres)+'nm)')

    surf_param_list = [calcParams(obj_i, pxlen, thres) for obj_i in surf_obj_list]
    img_param_list = [calcParams(obj_i, pxlen, thres) for obj_i in img_obj_list]
    
    # plot histograms of parameters
    N_param = len(surf_param_list[0])
    for i in range(N_param): # iterate over all parameters (mean, std, etc.)
        surf_dist = [surf_param_list[j][i] for j in range(len(surf_param_list))]
        img_dist = [img_param_list[j][i] for j in range(len(img_param_list))]
        
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True, sharey=True)
        bins = np.arange(min(surf_dist + img_dist), max(surf_dist + img_dist) + binwidth_list[i], binwidth_list[i])
        ax1.hist(surf_dist, color='b', bins=bins, edgecolor='black', linewidth=2)
        ax2.hist(img_dist, color='r', bins=bins, edgecolor='black', linewidth=2)
        ax1.set_xlim(bin_min_list[i], bin_max_list[i])
        ax2.set_xlim(bin_min_list[i], bin_max_list[i])
        ax2.set_xlabel(param_string_list[i])
        ax1.set_ylabel('frequency')
        ax2.set_ylabel('frequency')
        ax1.set_title('surface (found ' + str(len(surf_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm)')
        ax2.set_title('image (found ' + str(len(img_obj_list)) + ' of ' + str(N_part) +' particles, thres=' + str(thres)+'nm)')
        fig.suptitle(r'$h_{tip}=$' + str(h_tip) + r'$nm,\ a.r.=$' + str(ar_tip) + 
                     r'$,\ N_{pxl}=$' + str(N_pxl) + r'$,\ l_{px}=$' + str(pxlen) + 
                     r'$nm,\ thres=$' + str(thres) + r'$nm$')
    
    