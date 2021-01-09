import numpy as np
import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt

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
    
    V_mean = 4/3 * np.pi * np.exp(3*(2*R_mu + 3*R_sigma**2)/2)
    
    N_part = V / V_mean
    eff_cov = N_part / A
    
    return N_part, eff_cov

def plotTipDep(z, pxlen, h, R_mu, R_sigma, N_part_real):
    R_tip = np.linspace(0.01, 20, 10)
    N_part_est = []
    for R in R_tip:
        tip = mf.genParabolicTip(pxlen,h,r=R) 
        img = mph.grey_dilation(z, structure=-tip)
        N_part_est.append(partNum(img, pxlen, R_mu, R_sigma)[0])
    
    plt.figure()
    plt.plot(R_tip, np.array(N_part_est) / N_part_real, color='r', marker='o')
    plt.xlabel(r'$R_{tip} [nm]$')
    plt.ylabel(r'$N_{part,est} / N_{part,real}$')
    plt.title(r'$N_{part,real} = $' + str(N_part_real) + r'$, \mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$')
    plt.grid()

def plotNpartDep(Npx, pxlen, R_mean, R_std, R_mu, R_sigma):
    Npart = np.linspace(1000, 5000, 10)
    N_part_est = []
    for N in Npart:
        z=mf.genFlat(Npx)
        z = mf.genLogNormSph(z, pxlen,int(N),R_mean,R_std)
        N_part_est.append(partNum(z, pxlen, R_mu, R_sigma)[0])
    
    plt.figure()
    plt.plot(Npart, np.array(N_part_est) / Npart, color='r', marker='o')
    plt.xlabel(r'$N_{part,real}$')
    plt.ylabel(r'$N_{part,est} / N_{part,real}´$')
    plt.title(r'$\mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$' + r'$, R_{mean} = $' + str(R_mean) + r'$nm, R_{std} = $' + str(R_std) + r'$nm$')
    plt.grid()