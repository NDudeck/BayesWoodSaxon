from SALib.sample import saltelli
from SALib.analyze import sobol
from SALib.test_functions import Ishigami
import numpy as np
import get_wspot_schro as gws
import math
import scipy.integrate as integrate

# Define the model inputs
problem = {
    'num_vars': 6,
    'names': ["Vws", "Rws", "aws", "Vso", "Rso", "aso"],
    'bounds': [[-60, -40],
               [1.0, 4.0],
               [1.0, 3.0],
               [3.0, 10.0],
               [0.0, 1.5],
               [0.05, 2.0]]
}

# Generate samples
param_values = saltelli.sample(problem, 1000)
Y = np.zeros([param_values.shape[0]])


def get_wf_ws(theta):
    wf, e = gws.gwf(6, 3, 1, 0, 0, 1, 1, theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], -10, -1)

    wf = wf * math.sqrt(sf)

    ll_wf = -0.5 * np.sum(((wf_exp - wf) / sigma_wf) ** 2) + 0.5 * 100 * np.log(2 * np.pi * sigma_wf ** 2)
    ll_e = -0.5 * np.sum(((e_exp - e) / sigma_e) ** 2) + 0.5 * 100 * np.log(2 * np.pi * sigma_e ** 2)

    if np.isnan(e):
        e = -1000
        # print("NaN e")
    if np.isnan(wf).any():
        # print("NaN in rad.dat")
        return np.full((100, 1), -1000), e
    return ll_wf + ll_e


# Load Experimental data
wf_exp_in = np.genfromtxt(r"/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS "
                          r"4399/6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12.csv", skip_header=1,
                          delimiter=",")
wf_exp = wf_exp_in[:, 2]  # p = 1/2
# wf_exp = wf_exp_in[:, 1]  # p = 3/2
e_exp = -7.2499

sigma_e = 0.1 * e_exp
sigma_wf = 0.01 * max(wf_exp)

# sigma_wf = 0.04
# sigma_e = 0.3

sf = integrate.simps(np.square(wf_exp) * np.square(wf_exp_in[:, 0]), wf_exp_in[:, 0])

# Run model (example)
for i, X in enumerate(param_values):
    Y[i] = get_wf_ws(X)

# Perform analysis
Si = sobol.analyze(problem, Y, print_to_console=True)

# Print the first-order sensitivity indices
print(Si['S1'])
