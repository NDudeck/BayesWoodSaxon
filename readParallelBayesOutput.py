import numpy as np
import corner
import matplotlib.pyplot as plt

samples = np.genfromtxt(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\samples_bayes_out.csv')
chainV = np.genfromtxt(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\chainV_bayes_out.csv')
chainA = np.genfromtxt(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\chainA_bayes_out.csv')
chainR = np.genfromtxt(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\chainR_bayes_out.csv')
get_log_prob = np.genfromtxt(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\log_prob_bayes_out.csv')

chain = np.array([chainV, chainA, chainR])

(ndim, nwalkers, nsteps) = chain.shape

# make a corner plot with the posterior distribution
fig1, axs = plt.subplots(ndim, ndim)
corner.corner(samples, labels=["V", "a", 'r'], truths=[1, 1.25, 0.65], show_titles=True, fig=fig1)

axs[0, ndim-2].text(0.5, 0.8, 'A,Z = ' + str(6) + ', ' + str(3), horizontalalignment='center')
axs[0, ndim-2].text(0.5, 0.6, 'a,z = ' + str(1) + ', ' + str(0), horizontalalignment='center')
axs[0, ndim-2].text(0.5, 0.4, 'n,l,2j = ' + str(0) + ', ' + str(0) + ', ' + str(1), horizontalalignment='center')

axs[0, ndim-1].text(0.5, 0.2, 'ndim = ' + str(ndim), horizontalalignment='center')
axs[0, ndim-1].text(0.5, 0.4, 'nwalkers = ' + str(nwalkers), horizontalalignment='center')
axs[0, ndim-1].text(0.5, 0.6, 'nsteps = ' + str(nsteps), horizontalalignment='center')
axs[0, ndim-1].text(0.5, 0.8, 'sigma = ' + str(0.2), horizontalalignment='center')

axs[1, ndim-1].text(0.5, 1, 'Truths', horizontalalignment='center')
axs[1, ndim-1].text(0.5, 0.8, 'V0 = ' + str(1), horizontalalignment='center')
axs[1, ndim-1].text(0.5, 0.6, 'a0 = ' + str(1.25), horizontalalignment='center')
axs[1, ndim-1].text(0.5, 0.4, 'r0 = ' + str(0.65), horizontalalignment='center')

fig2, axs = plt.subplots(3)
for j in range(ndim):
    for i in range(nwalkers):
        axs[j].plot(range(nsteps), chain[j, i, :], linewidth=0.5)

axs[0].plot(range(nsteps), np.full((nsteps, 1), 1), 'r.')
axs[1].plot(range(nsteps), np.full((nsteps, 1), 1.25), 'r.')
axs[2].plot(range(nsteps), np.full((nsteps, 1), 0.65), 'r.')

fig3 = plt.figure()
plt.plot(range(nsteps), get_log_prob)

plt.show()
