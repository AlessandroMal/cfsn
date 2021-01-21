import numpy as np
import mapsfunction as mf
import parameters as par
import scipy.ndimage.morphology as mph
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from skimage.feature import peak_local_max

def reconstLogNorm(z, pxlen, R_tip, sigma_noise, thres, N_part, R_median_real, R_mu_real, R_sigma_real):
    z_obj_list, z_labeled = mf.identObj(z, thres)
    mf.plotThres(z, z_labeled, pxlen, 'found ' + str(len(z_obj_list)) + ' of ' + str(N_part) + ' particles, thres=' + str(thres)+'nm')

    R_list = np.array([np.max(obj_i)/2 for obj_i in z_obj_list])
    
    # R_list_index = peak_local_max(z, min_distance=1)
    # print(len(R_list_index))
    # mf.plotfalsecol(z, pxlen)
    # plt.scatter(pxlen*R_list_index[:,1] + 1/2*pxlen, pxlen*R_list_index[:,0] + 1/2*pxlen, color='red', marker='x')
    # R_list = z[R_list_index[:,0], R_list_index[:,1]] / 2
    
    logR_list = np.log(R_list)
    
    R_mu_fit, R_sigma_fit = norm.fit(logR_list)
    R_median_fit = np.exp(R_mu_fit)
    
    x = np.linspace(R_mu_fit - 3*R_sigma_fit, R_mu_fit + 3*R_sigma_fit, 1000)
    pdf_fit = norm.pdf(x, R_mu_fit, R_sigma_fit)
    pdf_real = norm.pdf(x, R_mu_real, R_sigma_real)
    
    plt.figure()
    plt.hist(logR_list, bins=12, density=True, edgecolor='black', linewidth=2, color='grey', alpha=0.5)
    plt.plot(x, pdf_fit, color='r', linewidth=3.5, label='empiric distribution (R_mu = {0:.3f}, R_sigma = {1:.3f})'.format(R_mu_fit, R_sigma_fit))
    plt.plot(x, pdf_real, color='green', linewidth=3.5, label='real distribution (R_mu = {0:.3f}, R_sigma = {1:.3f})'.format(R_mu_real, R_sigma_real))
    plt.xlabel(r'$\ln(R_{part}/nm)$')
    plt.ylabel('frequency')
    plt.title('gaussian fit to log data, ' + r'$\sigma_{noise} = $' + str(sigma_noise) + ' nm, ' + r'$R_{tip} = $' + str(R_tip) + ' nm')
    plt.legend(loc=1)
    plt.tight_layout()
    
    R_std_fit = np.sqrt((np.exp(R_sigma_fit**2)-1)*np.exp(2*R_mu_fit)*np.exp(R_sigma_fit**2))
    x_log = np.linspace(R_median_fit - 3*R_std_fit, R_median_fit + 3*R_std_fit, 1000)
    pdf_log_fit = lognorm.pdf(x_log, scale=R_median_fit, s=R_sigma_fit)
    pdf_log_real = lognorm.pdf(x_log, scale=R_median_real, s=R_sigma_real)
    
    plt.figure()
    plt.hist(R_list, bins=12, density=True, edgecolor='black', linewidth=2, color='grey', alpha=0.5)
    plt.plot(x_log, pdf_log_fit, color='r', linewidth=3.5, label='empiric distribution (R_median = {0:.3f}, R_sigma = {1:.3f})'.format(R_median_fit, R_sigma_fit))
    plt.plot(x_log, pdf_log_real, color='green', linewidth=3.5, label='real distribution (R_median = {0:.3f}, R_sigma = {1:.3f})'.format(R_median_real, R_sigma_real))
    plt.xlabel(r'$R_{part} [nm]$')
    plt.ylabel('frequency')
    plt.title('from fit calculated lognorm, ' + r'$\sigma_{noise} = $' + str(sigma_noise) + ' nm, ' + r'$R_{tip} = $' + str(R_tip) + ' nm')
    plt.legend(loc=1)
    plt.tight_layout()

    return R_median_fit, R_mu_fit, R_sigma_fit
    
def partNum(z, pxlen, R_mu, R_sigma):
    A = len(z)**2 * pxlen**2
    V = par.V(z, pxlen)
    
    V_mean = 4/3 * np.pi * np.exp(3*R_mu + 3**2*R_sigma**2/2)
    
    N_part = V / V_mean
    eff_cov = N_part / A
    
    return N_part, eff_cov

def plotNpartDep(Npx, pxlen, N_part_max, N_part_step, R_mean, R_std, R_mu, R_sigma):
    z = mf.genFlat(Npx)
    N_part = np.arange(0, N_part_max + 1, N_part_step)
    N_est_list = []
    for N in N_part:
        z, R_part_real = mf.genLogNormSolidSph(z, pxlen,N_part_step,R_mean,R_std)
        N_est_list.append(partNum(z,pxlen,R_mu,R_sigma)[0])
        print(N)
        mf.plotfalsecol(z,pxlen)
    
    plt.figure()
    plt.plot(N_part, np.array(N_est_list)/N_part, color='r', marker='o')
    plt.xlabel(r'$N_{part,real}$')
    plt.ylabel(r'$N_{est} / N_{real}$')
    plt.title(r'$\mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$' + r'$, R_{mean} = $' + str(R_mean) + r'$nm, R_{std} = $' + str(R_std) + r'$nm$')
    plt.grid()

def plotVfracDep(Npx, pxlen, N_part_max, N_part_step, R_mean, R_std, R_mu, R_sigma):
    z = mf.genFlat(Npx)
    V_real = 0
    V_frac_list = []
    N_part = np.arange(0, N_part_max + 1, N_part_step)
    for N in N_part:
        z, R_part_real = mf.genLogNormSolidSph(z, pxlen,N_part_step,R_mean,R_std)
        V_real += np.sum(4/3 * np.pi * R_part_real**3)
        V_frac_list.append(par.V(z,pxlen)/V_real)
        print(N)
    
    plt.figure()
    plt.plot(N_part, V_frac_list, color='r', marker='o')
    plt.xlabel(r'$N_{part,real}$')
    plt.ylabel(r'$V_{tot,est} / V_{tot,real}$')
    plt.title(r'$\mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$' + r'$, R_{mean} = $' + str(R_mean) + r'$nm, R_{std} = $' + str(R_std) + r'$nm$')
    plt.grid()
    
def plotGrowthDep(Npx, pxlen, N_part_max, N_part_step, R_mean, R_std, R_mu, R_sigma):
    z = mf.genFlat(Npx)
    RMS_list = []
    N_part = np.arange(0, N_part_max + 1, N_part_step)
    for N in N_part:
        z, R_part_real = mf.genLogNormSolidSph(z, pxlen,N_part_step,R_mean,R_std)
        RMS_list.append(np.std(z))
        print(N)
        
    plt.figure()
    plt.plot(N_part, RMS_list, color='r', marker='o')
    plt.xlabel(r'$N_{part,real}$')
    plt.ylabel(r'$RMS [nm]$')
    plt.title(r'$\mu_R = $' + str(R_mu) + r'$nm, \sigma_R = $' + str(R_sigma) + r'$nm$' + r'$, R_{mean} = $' + str(R_mean) + r'$nm, R_{std} = $' + str(R_std) + r'$nm$')
    plt.grid()

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

