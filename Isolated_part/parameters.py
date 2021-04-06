import numpy as np
import scipy.ndimage.morphology as mph
import skimage.measure as sk
import mapsfunction as mf
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

def skew(z,pxlen):
    std=h_std(z,pxlen)
    if std!=0: z_skew = (h_av(z**3,pxlen)-3*h_av(z,pxlen)*std**2-h_av(z,pxlen)**3)/ std**3
    else: z_skew=0
    return z_skew

def V(z,pxlen):
    somma=0
    
    for x in range(np.shape(z)[1]):
        for y in range(np.shape(z)[0]):
            somma+=z[y,x]
    #z ha gia unità fisiche
    return pxlen*pxlen*somma

def specArea(z):
#first order specific area
    A=0
    for x in range(np.shape(z)[1]-1):
        for y in range(np.shape(z)[0]-1):
    #        z_sorted=np.reshape(z[x:x+2, y:y+2], 4)
     #       z_sorted=np.sort(z_sorted)
            A+=np.linalg.norm(np.array([-z[y,x+1]+z[y,x], -z[y+1,x]+z[y,x], 1]))/2
            A+=np.linalg.norm(np.array([z[y+1,x]-z[y+1,x+1], z[y,x+1]-z[y+1,x+1], 1]))/2
    return A /((np.shape(z)[1]-1)*(np.shape(z)[0]-1))

def coverage(z, thres):
    N_tot = np.shape(z)[1]*np.shape(z)[0]
    N = np.sum(z>thres)
    cov = N/N_tot
    return cov

def h_av(z,pxlen):
    summ=0
    for x in range(np.shape(z)[1]):
        for y in range(np.shape(z)[0]):
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
        #      'skew': skew(z,pxlen),
              'V': V(z,pxlen),
          #    'specArea': specArea(z),
              'coverage': coverage(z,thres)}
    return params
    
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
    h=np.amax(z)-thres #già in unità fisiche l'asse z, come anche V
    ind= z>thres
    props=sk.regionprops( sk.label(ind) )
    e=props[0].eccentricity #eccentricità
    a=props[0].equivalent_diameter /2 #raggio equivalente cap
    return h, a*pxlen, props[0].area *pxlen**2, np.sum(z[ind]-thres)* pxlen**2, e

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
    
    model = LinearRegression().fit(logv, h)
    h_pred = model.predict(logv)
    msqer=mean_squared_error(h, h_pred)
    
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