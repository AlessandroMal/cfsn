import numpy as np
import matplotlib.pyplot as plt

plt.close('all')

mu = 1 # gaussian
sigma = 0.5 # gaussian

mean = np.exp(mu + sigma**2 / 2) # lognorm
std = np.sqrt((np.exp(sigma**2) - 1) * np.exp(2*mu + sigma**2)) # lognorm

inverse_mu = np.log(mean / np.sqrt(1 + std**2 / mean**2)) # recalculated gaussian
inverse_sigma = np.sqrt(np.log(1 + (std/mean)**2)) # reculaculated gaussian

#STARTING FROM GAUSSIAN PARAMETERS
plt.figure()
rand_values = np.random.lognormal(mu, sigma, size=10000)

x = np.linspace(mean - 3*std, mean + 3*std, 100)
pdf = 1 / (x * sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - mu)**2 / (2*sigma**2))
plt.hist(rand_values, bins=x, density=True)
plt.plot(x, pdf, color='r')
plt.xlabel('x')
plt.title('lognorm distribution from gaussian parameters: mu=' + str(mu) + ' , sigma=' + str(sigma))

#STARTING FROM LOGNORM PARAMETERS
plt.figure()
rand_values = np.random.lognormal(inverse_mu, inverse_sigma, size=10000)

x = np.linspace(mean - 3*std, mean + 3*std, 100)
pdf = 1 / (x * sigma * np.sqrt(2*np.pi)) * np.exp(-(np.log(x) - mu)**2 / (2*sigma**2))
plt.hist(rand_values, bins=x, density=True)
plt.plot(x, pdf, color='r')
plt.xlabel('x')
plt.title('lognorm distribution from lognormal parameters: mean = ' + str(mean) + ', std = ' + str(std))