import numpy as np
import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
from os.path import isfile
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error

def reconstLogNorm(z, pxlen, thres, N_part, R_mu_real, R_sigma_real):
    z_obj_list, z_labeled = mf.identObj(z, thres)
    mf.plotThres(z, z_labeled, pxlen, 'found ' + str(len(z_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm')
    
    R_list = [(np.max(np.shape(obj_i)))*pxlen/2 for obj_i in z_obj_list]
    
    R_mean = np.mean(R_list)
    R_std = np.std(R_list)
    
    R_mu = np.log(R_mean / np.sqrt(1 + R_std**2 / R_mean**2)) # recalculated gaussian
    R_sigma = np.sqrt(np.log(1 + (R_std/R_mean)**2)) # reculaculated gaussian
    
    x = np.linspace(R_mean - 3*R_std, R_mean + 3*R_std, 1000)
    pdf = 1 / (x * R_sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - R_mu)**2 / (2*R_sigma**2))
    pdf_real = 1 / (x * R_sigma_real * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - R_mu_real)**2 / (2*R_sigma_real**2))
    
    plt.figure()
    plt.hist(R_list, bins=12, density=True, edgecolor='black', linewidth=2, color='grey', alpha=0.5)
    plt.plot(x, pdf, color='r', linewidth=3.5, label='empiric distribution (R_mu = {0:.3f}, R_sigma = {1:.3f})'.format(R_mu, R_sigma))
    plt.plot(x, pdf_real, color='green', linewidth=3.5, label='real distribution (R_mu = {0:.3f}, R_sigma = {1:.3f})'.format(R_mu_real, R_sigma_real))
    plt.xlabel(r'$R_{part} [nm]$')
    plt.ylabel('frequency')
    plt.title('found ' + str(len(z_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm')
    plt.legend(loc=1)
    plt.tight_layout()

    return R_mean, R_std, R_mu, R_sigma
    
def partNum(z, pxlen, R_mu, R_sigma):
    A = len(z)**2 * pxlen**2
    V = par.V(z, pxlen)
    
    V_mean = 4/3 * np.pi * np.exp(3*R_mu + 3**2*R_sigma**2/2)
    
    N_part = V / V_mean
    eff_cov = N_part / A
    
    return N_part, eff_cov

def partDep(Npx, pxlen, step_sim, N_part_min, N_part_max, N_part_step, R_mean, R_std, par, firstmap='', usefile=False):
    N_part = np.linspace(np.log10(N_part_min), np.log10(N_part_max), N_part_step)
    N_part = np.round(10**N_part)
    N_part.astype(int, copy=False)
    
    R_mu    = np.log(R_mean / np.sqrt(1 + R_std **2 / R_mean**2))  # recalculated gaussian
    R_sigma = np.sqrt(np.log(1 + (R_std/R_mean)**2)) # reculaculated gaussian

    est_list = []
    wl_list = []
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
            print('Sim.',i+1,'N=',str(N)[:len(str(N))-2], end=' ')
            if usefile and isfile('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat'):
                print('map from file ...')
                z= np.loadtxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')
                out=open('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')
                head=out.readline()
                out.close()
                start=head.find('V=')+2
                V_real=float(head[start:head.find(';', start)])
            else:
                print('generating map ...')
                z, R_part_real = mf.genLogNormSolidSph(z,pxlen,int(N-N_prec),R_mean,R_std)
                V_real += np.sum(4/3 * np.pi * R_part_real**3)
                if( not(isfile('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat')) ): np.savetxt('maps/lognorm_'+str(Npx)+'_'+str(pxlen)+'_'+str(N)[:len(str(N))-2]+'.dat', z, header='V='+str(V_real)+'; Npx,pxlen,Npart in filename')
            N_prec=N
            
            if par=='N' or par=='N_part' or par=='Npart': est_list.append(partNum(z,pxlen,R_mu,R_sigma)[0])
            if par=='V' or par=='V_rel' or par=='V_frac': est_list.append(par.V(z,pxlen)/V_real)
            if par=='rms' or par=='RMS' or par=='std': est_list.append(np.std(z))
            wl_list.append(par.wavelength(z,pxlen,'x'))

    if par=='N' or par=='N_part' or par=='Npart': filename='N_relVsN_real.dat'
    if par=='V' or par=='V_rel' or par=='V_frac': filename='V_relVsV_real.dat'
    if par=='rms' or par=='RMS' or par=='std': filename='rmsVsN_real.dat'
    err=[]
    for i in range(len(N_part)):
        if step_sim==1: err.append(0)
        else: err.append(np.std(np.array(est_list[i::len(N_part)])))
        est_list[i]=np.mean(np.array(est_list[i::len(N_part)]))
        wl_list[i] =np.mean(np.array( wl_list[i::len(N_part)]))
        
    np.savetxt( filename, np.array([np.array(est_list[:len(N_part)]),N_part,np.array(err)],np.array(wl_list[:len(N_part)])),
              header=str(R_mu) + ' ' + str(R_sigma) + ' ' + str(R_mean) + ' ' + str(R_std) + ' ' + str(Npx) + ' ' + str(pxlen) + '\n' +
              r'$\mu_R$ $\sigma_R$ $R_{mean}$ $R_{std}$ $N_{px}$ $L_{px}$'+
              '\n wavelength on last line')
    print('data saved in '+ filename +'and in folder maps/')

def plotTipDep(z, pxlen, h, R_mu, R_sigma, N_part_real):
    R_tip = np.linspace(0.01, 20, 10)
    N_part_est = []
    for R in R_tip:
        tip = mf.genParabolicTip(pxlen,h,r=R) 
        img = mph.grey_dilation(z, structure=-tip)
        N_part_est.append(partNum(img, pxlen, R_mu, R_sigma)[0])
        print(R)
    
    plt.figure()
    plt.plot(R_tip, np.array(N_part_est) / N_part_real, color='r', marker='o')
    plt.xlabel(r'$R_{tip} [nm]$')
    plt.ylabel(r'$N_{part,est} / N_{part,real}$')
    plt.title(r'$N_{part,real} = $' + str(N_part_real) + r'$, \mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$')
    plt.grid()

def plotpar_fromfile(filename, xlab, ylab, xscale='linear', yscale='linear', a=0, b=0):
    est, Ntrue, err, wl = np.loadtxt(filename)
    
    out=open(filename)
    line1=out.readline()
    line2=out.readline()
    param=np.fromstring(line1[2:], sep=' ')
    name=line2[2:].split()
    out.close()
    
    plt.figure()
    if filename.find('N_rel')!=-1: est=est/Ntrue
    plt.errorbar(Ntrue, est, yerr=err)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(name[0]+'='+str(round(param[0], 2))+', '+name[1]+'='+str(round(param[1], 2))+', '+
              name[2]+'='+str(param[2])+', '+name[3]+'='+str(param[3])+', '+name[4]+'='+str(param[4])+', '+
              name[5]+'='+str(param[5]) )
    if a<b:
        passedA=False
        for i in range(len(Ntrue)):
            if Ntrue[i]>a and not(passedA):
                a=i
                passedA=True
                
            if Ntrue[i]>b: break
            else: b=i
        
        model  = LinearRegression().fit(np.log10(Ntrue[a:b]), np.log10(est[a:b]) )
        plt.plot(np.log10(Ntrue[a:b]),model.predict(Ntrue[a:b]), color='red')
        print('coefficient:',model.coef_)
        print('mean squared error:', mean_squared_error(np.log10(est[a:b]), model.predict(np.log10(Ntrue[a:b]))) )
    plt.show()