import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import math
import get_wspot_schro as gws


# Initialize physical parameters
A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 1, 1]  # Desired quantum numbers
theta_brida_p1 = [-41.80, 3.18, 0.85, 2 * 2.45, 1.35, 0.17]
theta_brida_p3 = [-69.55, 1.89, 1.17, 2 * 2.13, 2.36, 0.21]
emin = -90
emax = -1


wf_exp_in = np.genfromtxt(r"/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS "
                          r"4399/6Li_7Li_overlaps_NNLOopt_Nmax12.csv", skip_header=1,
                          delimiter=",")
wf_exp_p1 = wf_exp_in[:, 3]
wf_exp_p1_err = wf_exp_in[:, 4]

wf_exp_p3 = wf_exp_in[:, 1]
wf_exp_p3_err = wf_exp_in[:, 2]

sf_p1 = integrate.simps(np.square(wf_exp_p1) * np.square(wf_exp_in[:, 0]), wf_exp_in[:, 0])
sf_p3 = integrate.simps(np.square(wf_exp_p3) * np.square(wf_exp_in[:, 0]), wf_exp_in[:, 0])

wf_br_p1, e_1, _ = gws.gwf(A, Z, a, z,
                        nl2j[0], nl2j[1], 1,
                        theta_brida_p1[0], theta_brida_p1[1], theta_brida_p1[2], theta_brida_p1[3], theta_brida_p1[4], theta_brida_p1[5], emin, emax)
wf_br_p3, e_3, _ = gws.gwf(A, Z, a, z,
                        nl2j[0], nl2j[1], 3,
                        theta_brida_p3[0], theta_brida_p3[1], theta_brida_p3[2], theta_brida_p3[3], theta_brida_p3[4], theta_brida_p3[5], emin, emax)

wf_br_p1 = wf_br_p1 * math.sqrt(sf_p1)
wf_br_p3 = wf_br_p3 * math.sqrt(sf_p3)


plt.plot(np.linspace(0, 10, 100), wf_br_p1, label='Brida p=1/2', c='tab:blue')
plt.plot(np.linspace(0, 10, 100), wf_br_p3, label='Brida p=3/2', c='tab:orange')
plt.errorbar(np.linspace(0, 10, 100), wf_exp_p1, yerr=wf_exp_p1_err, label='SA-NCSM p=1/2', c='tab:blue', ecolor='gray', elinewidth=1, marker='o', markersize=3, linestyle='none')
plt.errorbar(np.linspace(0, 10, 100), wf_exp_p3, yerr=wf_exp_p3_err, label='SA-NCSM p=3/2', c='tab:orange', ecolor='gray', elinewidth=1, marker='o', markersize=3, linestyle='none')
plt.legend()
plt.xlabel('r')
plt.ylabel('Relative Motion Wavefunction $\Psi^{-3/2}$')
plt.title('Brida vs SA-NCSM')
plt.savefig(r'/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS 4399/wf.jpg')

print(e_1)
print(e_3)