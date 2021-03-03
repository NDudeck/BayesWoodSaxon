import subprocess
import time
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
import get_wspot_schro as gws


# Method for running wood-saxon executable and grabbing output wavefunction
def get_wf_ws(theta):
    wfr,wf,e = gws.get_wspot_wf(A, Z, a, z, nl2j[0], nl2j[1], nl2j[2], theta[0], theta[1], theta[2])

    wf_out = np.array((wfr,wf)).T
    wf_out = wf_out[:100,:]
    # print(wf_out)
    # print(e)

    if np.isnan(wf).any():
        print("NaN in rad.dat")
        return np.full((100, 2, -np.inf))
    return wf_out,e


# Start run timer
start_time = time.time()

# Initialize physical parameters
A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 0, 1]  # Desired quantum numbers
# TODO both spins at same-ish time

# Load in sampling parameters
ndim = 3  # number of parameters in the model
nwalkers = 36  # number of MCMC walkers
nburn = 200  # "burn-in" period to let chains stabilize
nsteps = 10000  # number of MCMC steps to take
sigma = 0.2
N = 100  # Resolution of wavefunction. ws_schro needs to be recompiled to change

# Load Experimental data
# wf_exp_in = np.genfromtxt("6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12.csv", skip_header=1, delimiter=",")
wf_exp_in,e_exp = get_wf_ws([1, 1.25, 0.65])
wf_exp = np.square(wf_exp_in[:, 1] * wf_exp_in[:, 0] + np.random.normal(0, 0.02, 100))

# Plot experimental data
# plt.ion()
# figure, ax = plt.subplots(figsize=(8, 6))
# line1, = ax.plot(wf_exp_in[:, 0], wf_exp)
# line2 = ax.scatter(wf_exp_in[:, 0], wf_exp)
# plt.title("WS Psi vs Exp Psi", fontsize=25)
# plt.xlabel("r", fontsize=18)
# plt.ylabel("psi^2(r)", fontsize=18)

# Build parameter space theta [V_N, a_0, r_0]
min_theta = np.array([0.85, 1.1, 0.50])
max_theta = np.array([1.15, 1.4, 0.80])
# min_theta = np.array([0.5, 1.0, 0.45])
# max_theta = np.array([1.9, 1.9, 0.95])
mu_prior = 0.5 * (min_theta + max_theta)
volume_theta = np.prod(max_theta - min_theta)


def log_prior_uniform(theta):
    # Flat prior
    if np.logical_and(min_theta <= theta, theta <= max_theta).all():
        return np.log(1 / volume_theta)
    else:
        return -np.inf


def log_prior_gaussian(theta):
    # Flat prior
    if np.logical_and(min_theta <= theta, theta <= max_theta).all():
        return -0.5 * np.sum(((theta - mu_prior) / sigma) ** 2) - 0.5 * N * np.log(2 * np.pi * sigma ** 2)
    else:
        return -np.inf


def log_likelihood(theta):
    # print("Params " + str(theta))

    wf_ws,e_ws = get_wf_ws(theta)  # Grab wavefunction with provided theta
    wf_ws = np.square(wf_ws[:, 1] * wf_ws[:, 0])  # turn phi(r)/r to psi(r) then square

    # Update dynamic plot
    # line1.set_xdata(np.linspace(0, 10, 100))
    # line1.set_ydata(wf_ws)
    # figure.canvas.draw()
    # figure.canvas.flush_events()
    # time.sleep(0.1)

    try:
        # Gaussian Log Likelyhood
        # Sigma has to be < (2pi)^-0.5 (or about 0.39) to keep ll < 0
        ll_wf = -0.5 * np.sum(((wf_exp - wf_ws) / sigma) ** 2) + 0.5 * N * np.log(2 * np.pi * sigma ** 2)
        ll_e = -0.5 * np.sum(((e_exp - e_ws) / sigma) ** 2) + 0.5 * N * np.log(2 * np.pi * sigma ** 2)
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

        # Initialize sampler using paramters
        print("MCMC sampling using emcee (affine-invariant ensamble sampler) with {0} walkers".format(nwalkers))
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, pool=pool, moves=emcee.moves.DESnookerMove())

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

    tau = sampler.get_autocorr_time()
    print(tau)

    # make a corner plot with the posterior distribution
    fig1 = corner.corner(samples, labels=["V", "a", 'r'], truths=[1, 1.25, 0.65], show_titles=True)

    fig2, axs = plt.subplots(3)
    for j in range(ndim):
        for i in range(nwalkers):
            axs[j].plot(range(nsteps), sampler.chain[i, :, j], linewidth=0.5)

    axs[0].plot(range(nsteps), np.full((nsteps, 1), 1), 'r.')
    axs[1].plot(range(nsteps), np.full((nsteps, 1), 1.25), 'r.')
    axs[2].plot(range(nsteps), np.full((nsteps, 1), 0.65), 'r.')

    fig3 = plt.figure()
    plt.plot(range(nsteps), sampler.get_log_prob())

    plt.show()

    np.savetxt('samples_bayes_out.csv', samples)
    np.savetxt('chainV_bayes_out.csv', sampler.chain[:, :, 0])
    np.savetxt('chainA_bayes_out.csv', sampler.chain[:, :, 1])
    np.savetxt('chainR_bayes_out.csv', sampler.chain[:, :, 2])
    np.savetxt('log_prob_bayes_out.csv', sampler.get_log_prob())

if __name__ == '__main__':
    main()
    plt.show()



# To compile code
#  python3 -m numpy.f2py -c get_wspot_schro.f -m get_wspot_schro
# To run code
# python3 BayesWoodSaxonParallel.py
