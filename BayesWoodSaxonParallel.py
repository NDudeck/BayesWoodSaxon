import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import emcee
import scipy.integrate as integrate
import get_wspot_schro as gws
import math


# Method for running wood-saxon executable and grabbing output wavefunction
def get_wf_ws(theta):
    wf, e = gws.gwf(A, Z, a, z,
                    nl2j[0], nl2j[1], nl2j[2],
                    theta[0], theta[1], theta[2], theta[3], theta[4], theta[5], emin, emax)

    wf = wf * math.sqrt(sf)

    # norm = integrate.simps(np.square(wf_out[:, 1]) * np.square(wf_out[:, 0]), wf_out[:, 0])
    # print(norm)

    if np.isnan(e):
        e = -500
        # print("NaN e")
    if np.isnan(wf).any():
        # print("NaN in rad.dat")
        return np.full((100, 1), -500), e
    return wf, e


# Start run timer
start_time = time.time()

# Initialize physical parameters
A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0  # Proton is 1,1
nl2j = [0, 1, 1]  # Desired quantum numbers
# TODO both spins at same-ish time

# Load in sampling parameters
ndim = 6  # number of parameters in the model
nwalkers = 32  # number of MCMC walkers
nburn = 10000  # "burn-in" period to let chains stabilize
nsteps = 10000  # number of MCMC steps to take
N = 100  # Resolution of wavefunction. ws_schro needs to be recompiled to change
wfr = np.linspace(0, 10, 100)
np.random.seed(1)

# Build parameter space theta [V_ws, a_ws, r_ws, V_so, a_so, r_so,]
min_theta = np.array([-85, 0.5, 0.5, 0.5, 0.5, 0.05])
max_theta = np.array([-35, 3.5, 1.5, 10, 2.5, 0.5])
emin = -90
emax = - 1
# min_theta = np.array([-42.00, 0.80, 3.15, 0.05, 0.05, 1.00])
# max_theta = np.array([-40.00, 0.90, 3.20, 20, 0.5, 3.00])
volume_theta = np.prod(max_theta - min_theta)
mu_prior = 0.5 * (min_theta + max_theta)

# Load Experimental data
wf_exp_in = np.genfromtxt(r"/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS "
                          r"4399/6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12 2021-02-05 21_00_39.csv", skip_header=1,
                          delimiter=",")
if nl2j[2] == 3:
    wf_exp = wf_exp_in[:, 1]
elif nl2j[2] == 1:
    wf_exp = wf_exp_in[:, 2]
else:
    print("Need valid 2j")
    raise ValueError
e_exp = -7.2499

sigma_e = 0.10 * e_exp
sigma_wf = 0.10 * max(wf_exp)
# sigma_wf = 0.04
# sigma_e = 0.3


# sf = integrate.simps(np.square(wf_exp) * np.square(wf_exp_in[:, 0]), wf_exp_in[:, 0])


# For synthetic data
sf = 1
if nl2j[2] == 1:
    theta_truths = [-41.80, 3.18, 0.85, 2*2.45, 1.35, 0.17]
else:
    theta_truths = [-69.55, 1.89, 1.17, 2*2.13, 2.36, 0.21]

wf_exp_in, e_exp = get_wf_ws(theta_truths)
wf_exp = wf_exp_in + np.random.normal(0, 0.0, 100)

def log_prior_uniform(theta):
    # Flat prior
    if np.logical_and(min_theta <= theta, theta <= max_theta).all():
        return np.log(1 / volume_theta)
    else:
        return -np.inf


def log_prior_gaussian(theta):
    # Gaussian prior
    if np.logical_and(min_theta <= theta, theta <= max_theta).all():
        return -0.5 * np.sum(((theta - mu_prior) / sigma_wf) ** 2) - 0.5 * N * np.log(2 * np.pi * sigma_wf ** 2)
    else:
        return -np.inf


def log_likelihood(theta):
    # print("Params " + str(theta))
    weight = 1
    wf_ws, e_ws = get_wf_ws(theta)  # Grab wavefunction with provided theta
    # wf_ws = wf_ws[:, 1] * wf_ws[:, 0]  # turn phi(r)/r to psi(r) then square

    # Update dynamic plot
    # line1.set_xdata(np.linspace(0, 10, 100))
    # line1.set_ydata(wf_ws)
    # figure.canvas.draw()
    # figure.canvas.flush_events()
    # time.sleep(0.1)

    try:
        # Gaussian Log Likelihood
        # Sigma has to be < (2pi)^-0.5 (or about 0.39) to keep ll < 0
        ll_wf = -0.5 * np.sum(((wf_exp - wf_ws) / sigma_wf) ** 2) + 0.5 * N * np.log(2 * np.pi * sigma_wf ** 2)
        ll_e = weight * -0.5 * np.sum(((e_exp - e_ws) / sigma_e) ** 2) + 0.5 * N * np.log(2 * np.pi * sigma_e ** 2)
        # print(ll_wf + ll_e)
        return ll_wf + ll_e
    except ValueError:
        print("value Error")
        return -np.inf


def log_posterior(theta):
    return log_prior_uniform(theta) + log_likelihood(theta)


def main():
    with Pool() as pool:
        # Start at random locations within the prior volume
        starting_guesses = np.random.uniform(min_theta, max_theta, (nwalkers, ndim))

        # Initialize sampler using parameters
        print("MCMC sampling using emcee (affine-invariant ensemble sampler) with {0} walkers".format(nwalkers))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, pool=pool,  moves=emcee.moves.DESnookerMove())
        #
        # "burn-in" period; save final positions and then reset
        pos, prob, state = sampler.run_mcmc(starting_guesses, nburn, progress=True)
        sampler.reset()

        # sampling period
        sampler.run_mcmc(pos, nsteps, progress=True)

    print("Mean acceptance fraction: {0:.3f} (in total {1} steps)"
          .format(np.mean(sampler.acceptance_fraction), nwalkers * nsteps))

    print("time elapsed: {:.2f}s".format(time.time() - start_time))

    # discard burn-in points and flatten the walkers; the shape of samples is (nwalkers*nsteps, ndim)
    samples = sampler.chain.reshape((-1, ndim))

    # save the output in a compressed and easy to read-in format for analysis
    np.savez_compressed(r'/mnt/c/Users/noeld/eclipse-workspace/get_wspot_schro/src/parallelout/bayes_out.npz',
                        samples=samples, chain=sampler.chain, log_prob=sampler.get_log_prob(),
                        azaz=[A, Z, a, z], nl2j=nl2j, sigma_wf=sigma_wf, sigma_e=sigma_e,
                        wfr=wfr, wf_exp=wf_exp, e_exp=e_exp, sf=sf, e=[emin, emax])

    tau = sampler.get_autocorr_time()
    print(tau)


if __name__ == '__main__':
    main()
    plt.show()

# To compile code
#  python3 -m numpy.f2py -c get_wspot_schro.f -m get_wspot_schro
# To run code
# python3 BayesWoodSaxonParallel.py
# python3 "/mnt/c/Users/noeld/OneDrive - Louisiana State University/PHYS 4399/Code/BayesWoodSaxon/BayesWoodSaxonParallel.py"
