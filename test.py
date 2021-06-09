import math

import get_wspot_schro as gws
from scipy import integrate
import matplotlib.pyplot as plt


# Initialize physical parameters
import numpy as np

A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 1, 1]  # Desired quantum numbers
theta = [-39.11, 3.30, 0.75, 4.55, 1.12, 0.88]
emin = -90
emax = -1

wf, e, rms = gws.gwf(A, Z, a, z, nl2j[0], nl2j[1], nl2j[2],
                           theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], emin, emax)

# Calculate root-mean-squared
wf_exp_in = np.genfromtxt(r"/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS "
                          r"4399/6Li_7Li_overlaps_NNLOopt_Nmax12.csv", skip_header=1,
                          delimiter=",")
sf = integrate.simps(wf_exp_in[:, 3] * np.square(np.linspace(0, 10, 100)) * wf_exp_in[:, 3], np.linspace(0, 10, 100))
# make a corner plot with the posterior distribution
fig1 = plt.figure(1)
plt.plot(np.linspace(0, 20, 200), np.log(wf * math.sqrt(sf)))
plt.errorbar(np.linspace(0, 10, 100), np.log(wf_exp_in[:, 3]), yerr=wf_exp_in[:, 4], label='SA-NCSM', ms=5,
             c='gray', fmt='.', zorder=1)
plt.savefig(r'/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS 4399/testplot.jpg')

print(e)
print(rms)
print(sf)

