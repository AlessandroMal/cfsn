# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 14:09:37 2021

@author: alego
"""

def partDep(Npx, pxlen, step_sim, N_part_min, N_part_max, N_part_step, R_mu, R_sigma, firstmap='', usefile=False, savefile=False):
    N_part = np.linspace(np.log10(N_part_min), np.log10(N_part_max), N_part_step)
    N_part = np.round(10**N_part)
    N_part.astype(int, copy=False)
    
    N_est = np.array([])
    V_est = np.array([])
    rms_est = np.array([])
    h_est = np.array([])
    h_top = np.array([])
    
#    L_corr_est = np.array([])
#    L_corr_err = np.array([])
#    alfa_est=np.array([])
#    alfa_err=np.array([])
    
    C_true=[]
    G_true=[]
    C2_true=[]
    G2_true=[]
    C_list=[]
    G_list=[]
    for i in range(step_sim):
        if firstmap=='': #initialize map
            z = mf.genFlat(Npx)
            N_prec=0
            V_real=0
        else:
            z = np.loadtxt(firstmap)
            dotpos=firstmap.find('.')
            start=dotpos-1
            for i in range(dotpos-1):
                if firstmap[start]=='_': break
                start-=1
            N_prec=int(firstmap[start+1:dotpos])
            
            out=open(firstmap)
            head=out.readline()
            out.close()
            start=head.find('V=')+2
            V_real=float(head[start:head.find(';', start)])

        for N in N_part:
            print('Sim.',i+1,'; N=',str(N)[:len(str(N))-2], end=' ')
            if usefile and isfile('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat'):
                print('map from file ...', end=' ')
                z= np.loadtxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')
                out=open('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')
                head=out.readline()
                out.close()
                start=head.find('V=')+2
                V_real=float(head[start:head.find(';', start)])
            else:
                print('generating map ...', end=' ')
                z, R_part_real = mf.genLogNormSolidSph(z,pxlen,int(N-N_prec),R_mu,R_sigma)
                V_real += np.sum(4/3 * np.pi * R_part_real**3)
                if savefile: np.savetxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat', z, header='V='+str(V_real)+'; Npx,pxlen,Npart in filename')
            N_prec=N
            
            print('computing parameters ...')
            N_est=np.append(N_est, partNum(z,pxlen,R_mu,R_sigma)[0]/N)
            V_est=np.append(V_est, par.V(z,pxlen)/V_real)
            rms_est=np.append(rms_est, np.std(z))
            h_est=np.append(h_est, np.mean(z))
            h_top=np.append(h_top, np.amax(z))
            
            #print('computing correlations ...')
            l,C=par.C_profile(z, 2, 800)
#            l*=pxlen
            
            C_list.append(C)
            
#            slope, intrc, r_val, p_val, errlinregr= linregress(l, C)
#            L_corr_est=np.append( L_corr_est, - slope**-1 )
#            L_corr_err=np.append( L_corr_err, errlinregr/slope**2)
#            L_corr_est=np.append(L_corr_est, 1)
            
            l,C=par.G_profile(z, 2, 800)
#            l*=pxlen

            G_list.append(C)
            
#            Lforfit=np.ones(len(l))*L_corr_est[-1]
#            rmsforfit=np.ones(len(l))*rms_est[-1]
#            opt,cov= curve_fit(G_gauss_prof, (l, Lforfit, rmsforfit), C)
#            alfa_est=np.append(alfa_est, *opt)
#            alfa_err=np.append(alfa_err, *cov[0])
        if i==0: 
            for el in C_list:
                C_true.append(el)
                C2_true.append(el**2)
            for el in G_list:
                G_true.append(el)
                G2_true.append(el**2)
        else:
            for i in range(len(C_true)):
                C_true[i]=C_true[i]+C_list[i]
                C2_true[i]=C2_true[i]+C_list[i]**2
            for i in range(len(G_true)):
                G_true[i]=G_true[i]+G_list[i]
                G2_true[i]=G2_true[i]+G_list[i]**2
        C_list.clear()
        G_list.clear()

    

    filename=['N_relVsN.dat', 'V_relVsN.dat', 'rmsVsN.dat', 'hVsN.dat', 'maxhVsN.dat']
    est=[N_est, V_est, rms_est, h_est, h_top]
    
    for j in range(len(est)):
      #  if j<4:
        err=np.array([])
        for i in range(len(N_part)):
            if step_sim==1: err=np.append(err, 0)
            else: err=np.append(err, np.std(est[j][i::len(N_part)]))
            est[j][i]=np.mean(est[j][i::len(N_part)])
#        else:
 #           if j==4: err=L_corr_err
  #          if j==5: err=alfa_err
   #         for i in range(len(N_part)):
    #            est[j][i]=np.mean(est[j][i::len(N_part)])
     #           err[i]=   np.mean(err[i::len(N_part)]) /np.sqrt(step_sim -1)

    
        np.savetxt(filename[j], np.array([ est[j][:len(N_part)], N_part, err ]),
                   header=str(R_mu) + ' ' + str(R_sigma) + ' ' + str(Npx) + ' ' + str(pxlen) + '\n' +
                   r'$\mu _R$ $\sigma _R$ $N_{px}$ $L_{px}$')
        print('data saved in '+ filename[j] +' and in folder maps/')
    
    C_true=np.array(C_true)/step_sim
    G_true=np.array(G_true)/step_sim
    C2_true=np.array(C2_true)/step_sim - C_true**2
    G2_true=np.array(G2_true)/step_sim - G_true**2
    np.savetxt('correlations.dat', np.array([ l*pxlen, C_true, G_true, C2_true, G2_true]), fmt='%s')
   


def identObj(z, thres):
    z_binary = (z>thres) #vero se elemento e piu grande della soglia
    z_labeled = scipy.ndimage.label(z_binary)[0] #numera le singole particelle
    objInd = scipy.ndimage.find_objects(z_labeled) #trova gli indici delle particelle
    
    obj = [z[i] for i in objInd]
    print('identObj ha trovato ' + str(len(obj)) + ' particelle')

    return obj


def Cm(z, bin_size, px_cut):
    r=np.arange(0, min(px_cut, np.shape(z)[1]/2, np.shape(z)[0]/2), bin_size)
    hh=np.zeros(len(r)-1)

    for x0 in range(np.shape(z)[1]):
        for y0 in range(np.shape(z)[0]):
            for i in range(len(r)-1):
                m=np.where( (z[max(0,int(y0-r[i+1])-1):min(np.shape(z)[0],int(y0+r[i+1])-1), max(0,int(x0-r[i+1])-1):min(np.shape(z)[1],int(x0+r[i+1])-1)] >=r[i]) and
                            (z[max(0,int(y0-r[i+1])-1):min(np.shape(z)[0],int(y0+r[i+1])-1), max(0,int(x0-r[i+1])-1):min(np.shape(z)[1],int(x0+r[i+1])-1)] <r[i+1]) )
                m=np.ma.masked_array(z[max(0,int(y0-r[i+1])-1):min(np.shape(z)[0],int(y0+r[i+1])-1), max(0,int(x0-r[i+1])-1):min(np.shape(z)[1],int(x0+r[i+1])-1)],
                    mask=[~(np.sqrt( (np.ndenumerate(a))[0]-y0)**2 + (list(np.ndenumerate(a))[0]-x0)**2)>=r[i] and  np.sqrt( (list(np.ndenumerate(a))[0]-y0)**2 + ( list(np.ndenumerate(a))[1]-x0)**2)<r[i+1]) for a in z])

                hh[i]+=np.mean(z[y0,x0]*m)
        if (x0+1)%int(np.shape(z)[1]/20) ==0: print(round((x0+1)/np.shape(z)[1]*100, 0), '%')

    return r+bin_size/2, hh/ (np.shape(z)[1]*np.shape(z)[0])



inp=open('sphOnsph.dat') #leggo i dati
file=inp.read()
inp.close()

start=file.find('[') +1
end  =file.find(']', start)
R_tip=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
N_obj=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
h=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
r_part=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
A=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
V=np.fromstring(file[start:end].replace('\n',' '), sep=' ')
start=file.find('[', end+1) +1
end  =file.find(']', start)
dV=np.fromstring(file[start:end].replace('\n',' '), sep=' ')

R_rel=R_tip / h[0]

x_prof= []
z_prof= []
for i in range(len(R_rel)+1):
    start=file.find('[', end+1) +1
    end  =file.find(']', start)
    x_prof.append(np.fromstring(file[start:end].replace('\n',' '), sep=' '))
    start=file.find('[', end+1) +1
    end  =file.find(']', start)
    z_prof.append(np.fromstring(file[start:end].replace('\n',' '), sep=' '))

plt.figure()
plt.title(r'$V_{cap}$ error funciton')
plt.plot(R_rel,dV)
plt.xlabel(r'$ R_{tip} / R_{part} $')
plt.ylabel(r'$ \Delta V $')
plt.grid()
plt.show()

plt.figure()
plt.title('Deformation')
plt.plot(R_rel,r_part/h)
plt.xlabel(r'$ R_{tip} / R_{part} $')
plt.ylabel(r'$ R_{dil} / R_{true} $')
plt.grid()
plt.show()

'''
plt.figure()
plt.title(r'particle profile at $d/L_{px}=$'+str(round(2*R/pxlen[-1], 3)))
for i in range(0,len(z_prof)):
    if i==0: plt.plot(x_prof[i], z_prof[i],  label='digital map')
    else: plt.plot(x_prof[i], z_prof[i], label=r'$ R_{tip} / R_{part} = $'+str(round(R_rel[i-1], 2)))

plt.xlabel('r [npx]')
plt.xlabel('z [npx]')
plt.legend()
plt.grid()
plt.show()
'''




    maxh=0
    posmax=np.append(posmax, 0)
    R0=par.capPar(obj[0],pxlen,thres)[0]
    for x in range(np.shape(obj[0])[0]):
        height.append(np.amax(obj[0][0:,x])/R0)
        if maxh<height[-1]:
            maxh=height[-1]
            posmax[-1]=x
    profiles.append(np.array(height))
    height.clear()


np.savetxt('sphOnsph.dat', np.array([R_tip, N_obj , h_arr, r_arr, A_arr, V_arr, dV_arr]), fmt='%s',
           header='on rows: R_tip, N_obj, h, r, A, V, dV')

out=open('profiles.dat', 'w')
for i in range(len(R_tip)):
    out.write(str(pxlen * (np.arange(0,len(profiles[i])) - posmax[i]) )+'\n')
    out.write(str(profiles[i])+'\n'+'\n')
out.close()

def VfracDep(Npx, pxlen, step_sim, N_part_min, N_part_max, N_part_step, R_mean, R_std, V_real=0, loadmap=''): 
    if loadmap=='':
        z = mf.genFlat(Npx)
        N_prec=0
    else:
        z = np.loadtxt(loadmap)
        dotpos=loadmap.find('.')
        start=dotpos-1
        for i in range(dotpos-1):
            if loadmap[start]=='_': break
            start-=1
        N_prec=int(loadmap[start+1:dotpos])
    
    
    N_part = np.linspace(np.log10(N_part_min), np.log10(N_part_max), N_part_step)
    N_part = np.round(10**N_part)
    N_part.astype(int, copy=False)
    V_rel_list = []
    
    R_mu    = np.log(R_mean / np.sqrt(1 + R_std **2 / R_mean**2))  # recalculated gaussian
    R_sigma = np.sqrt(np.log(1 + (R_std/R_mean)**2)) # reculaculated gaussian
    
    for i in range(step_sim):
        print('simulation',i+1)
        for N in N_part:
            print('N=',str(N)[:len(str(N))-2], '...')
            z, R_part_real = mf.genLogNormSolidSph(z,pxlen,int(N-N_prec),R_mean,R_std)
            if( not(isfile('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')) ): np.savetxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat', z, header='Npx,pxlen,Npart in filename')
            N_prec=N
            V_real += np.sum(4/3 * np.pi * R_part_real**3)
            V_rel_list.append(par.V(z,pxlen)/V_real)
    
    V_std=[]
    for i in range(len(N_part)):
        if step_sim==1: V_std.append(0)
        else: V_std.append(np.std(np.array(V_rel_list[i::len(N_part)])))
        V_rel_list[i]=np.mean(np.array(V_rel_list[i::len(N_part)]))
    
    np.savetxt( 'Vap_onVtrueVsN.dat', np.array([np.array(V_rel_list[:len(N_part)]),N_part,np.array(V_std)]),
              header=str(R_mu) + ' ' + str(R_sigma) + ' ' + str(R_mean) + ' ' + str(R_std) + ' ' + str(Npx) + ' ' + str(pxlen) + '\n' +
              r'$\mu_R$ $\sigma_R$ $R_{mean}$ $R_{std}$ $N_{px}$ $L_{px}$' )
    print('data saved in file Vap_onVtrueVsN.dat and in folder maps/')
    
def GrowthDep(Npx, pxlen, step_sim, N_part_min, N_part_max, N_part_step, R_mean, R_std, loadmap=''):
    if loadmap=='':
        z = mf.genFlat(Npx)
        N_prec=0
    else:
        z = np.loadtxt(loadmap)
        dotpos=loadmap.find('.')
        start=dotpos-1
        for i in range(dotpos-1):
            if loadmap[start]=='_': break
            start-=1
        N_prec=int(loadmap[start+1:dotpos])
    
    N_part = np.linspace(np.log10(N_part_min), np.log10(N_part_max), N_part_step)
    N_part = np.round(10**N_part)
    N_part.astype(int, copy=False)
    RMS_list = []
    
    R_mu    = np.log(R_mean / np.sqrt(1 + R_std **2 / R_mean**2))  # recalculated gaussian
    R_sigma = np.sqrt(np.log(1 + (R_std/R_mean)**2)) # reculaculated gaussian

    for i in range(step_sim):
        print('simulation',i+1)
        for N in N_part:
            print('N=',str(N)[:len(str(N))-2], '...')
            z, R_part_real = mf.genLogNormSolidSph(z,pxlen,int(N-N_prec),R_mean,R_std)
            if( not(isfile('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')) ): np.savetxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat', z, header='Npx,pxlen,Npart in filename')
            N_prec=N
            RMS_list.append(np.std(z))

    RMS_std=[]
    for i in range(len(N_part)):
        if step_sim==1: RMS_std.append(0)
        else: RMS_std.append(np.std(np.array(RMS_list[i::len(N_part)])))
        RMS_list[i]=np.mean(np.array(RMS_list[i::len(N_part)]))

    np.savetxt( 'rmsVsN.dat', np.array([np.array(RMS_list[:len(N_part)]),N_part,np.array(RMS_std)]),
              header=str(R_mu) + ' ' + str(R_sigma) + ' ' + str(R_mean) + ' ' + str(R_std) + ' ' + str(Npx) + ' ' + str(pxlen) + '\n' +
              r'$\mu_R$ $\sigma_R$ $R_{mean}$ $R_{std}$ $N_{px}$ $L_{px}$' )
    print('data saved in file rmsVsN.dat and in folder maps/')

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


import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import numpy as np
import matplotlib.pyplot as plt
#from scipy.ndimage import gaussian_filter
#import time

plt.close('all') #chiude tutte le figure aperte
#probabilmente devo avere Npxtip<280 circa (16gb ram con 25% circa occupato)
#278 funziona, 296 no

def volsphere(R): return 5/3*np.pi*R**3
def volsemisph(R): return 2/3*np.pi*R**3

filename='conv_sphere.dat'
L=100 #lato mappa
q=1.3/10 #rapporto raggio/lato
Npx=np.linspace(60,400,8)

pxlen=L/Npx
thres=0 #soglia
#-----------------------------------------------
R=q*L
htip=2*R*1.02
R_tip=np.linspace(R/5,3*R, 8)

V_red_dil = []
height = []
profiles = []
posmax = []

z=mf.genFlat(int(Npx[-1]))
z=mf.genSphere(z,pxlen[-1],np.array([L/2, L/2]),np.array([R]))
obj = mf.identObj(z,thres)[0]
V_red_calc= par.V(obj, pxlen[-1]) / volsphere(R)

maxh=0
posmax.append(0)
for x in range(np.shape(obj)[0]):
    height.append(np.amax(obj[0:,x])/R)
    if maxh<height[-1]:
        maxh=height[-1]
        posmax[-1]=x
profiles.append(np.array(height))
height.clear()

for rtip in R_tip:
    for i in range(len(Npx)): #iterazione su risoluzione
        print('R_tip=', rtip , ' Npx=', int(Npx[i]))
        z=mf.genFlat(int(Npx[i]))
        z=mf.genSphere(z,pxlen[i],np.array([L/2, L/2]),np.array([R]))
        
        tip=mf.genSemisphTip(pxlen[i],htip,r=rtip)
        z = mph.grey_dilation(z, structure=-tip)
        if rtip==R_tip[-1] and i==0: mf.plotfalsecol(z,pxlen[i])
        obj = mf.identObj(z,thres)[0]
        V_red_dil.append(par.V(obj, pxlen[i]) / volsphere(R))
        
        if i==len(Npx)-1:
            maxh=0
            posmax.append(0)
            for x in range(np.shape(obj)[0]):
                height.append(np.amax(obj[0:,x])/R)
                if maxh<height[-1]:
                    maxh=height[-1]
                    posmax[-1]=x
            profiles.append(np.array(height))
            height.clear()

R_tip=np.append(R_tip, R)
print('printing data in',filename)
out=open(filename, 'w')
out.write('r_tip (R_part in last pos): '+str(R_tip)+'\n')
out.write('Lpx: '+str(pxlen)+'\n')
out.write('digital map volume: ['+str(np.array(V_red_calc))+']\n')
out.write('dilation map volume: '+str(np.array(V_red_dil))+'\n')
out.write('x_surf, surface profile, x dilated, map profile\n')
out.write(str(pxlen[-1] * (np.arange(0,len(profiles[0])) - posmax[0]) )+'\n')
out.write(str(profiles[0])+'\n'+'\n')
R_tip=np.delete(R_tip, -1)
for i in range(1,len(R_tip)+1):
    out.write(str(pxlen[-1] * (np.arange(0,len(profiles[i])) - posmax[i]) )+'\n')
    out.write(str(profiles[i])+'\n'+'\n')
out.close()