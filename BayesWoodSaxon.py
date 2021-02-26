import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner


# Method for running wood-saxon executable and grabbing output wavefunction
def get_wf_ws(theta):
    # Write file input with given parameters
    l7_dia = open(r"./io/l7.dai", 'w')
    l7_dia.write(str(A) + " " + str(Z) + " " + str(a) + " " + str(z) + "\n\n")
    l7_dia.write(' '.join(map(str, nl2j)) + ' ' + ' '.join(map(str, theta)))
    l7_dia.close()

    # Call executable without outputting into console
    subprocess.run(r"./io/wspot_schro.exe l7", cwd=r'./io', stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    try:
        wf = np.genfromtxt(r"./io/rad.dat").astype(float)  # Grab output
        wf = wf[:100, :]  # Trim to match input
        if np.isnan(wf).any():
            print("NaN in rad.dat")
            return np.full((100, 2, -np.inf))
        return wf
    except ValueError:
        print("Error Reading in rad.dat")
        # If parameters aren't good, return flat wavefunction
        return np.full((100, 2, -np.inf))


# Start run timer
start_time = time.time()

# Initialize physical parameters
A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 1, 1]  # Desired quantum numbers
# TODO both spins at same-ish time

# Load in sampling parameters
ndim = 3  # number of parameters in the model
nwalkers = 12  # number of MCMC walkers
nburn = 10  # "burn-in" period to let chains stabilize
nsteps = 200  # number of MCMC steps to take
sigma = 0.45
N = 100  # Resolution of wavefunction. ws_schro needs to be recompiled to change

# Load Experimental data
# wf_exp_in = np.genfromtxt("6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12.csv", skip_header=1, delimiter=",")
wf_exp_in = get_wf_ws([1, 1.25, 0.65])
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
max_theta = np.array([1.35, 1.4, 0.80])
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


def log_likelihood(theta, wf_exp):
    # print("Params " + str(theta))

    wf_ws = get_wf_ws(theta)  # Grab wavefunction with provided theta
    wf_ws = np.square(wf_ws[:, 1] * wf_ws[:, 0])  # turn phi(r)/r to psi(r) then square

    # Update dynamic plot
    # line1.set_xdata(np.linspace(0, 10, 100))
    # line1.set_ydata(wf_ws)
    # figure.canvas.draw()
    # figure.canvas.flush_events()
    # time.sleep(0.1)

    try:
        # Gaussian Log Likelyhood
        # Sigma has to be > (2pi)^-0.5 (or about 0.39) to keep ll < 0
        ll = -0.5 * np.sum(((wf_exp - wf_ws) / sigma) ** 2) - 0.5 * N * np.log(2 * np.pi * sigma ** 2)
        # print(ll)
        return ll
    except ValueError:
        print("value Error")
        return -np.inf


def log_posterior(theta, wf_exp):
    return log_prior_uniform(theta) + log_likelihood(theta, wf_exp)


# Start at random locations within the prior volume
starting_guesses = np.random.uniform(min_theta, max_theta, (nwalkers, ndim))

# Initialize sampler using paramters
print("MCMC sampling using emcee (affine-invariant ensamble sampler) with {0} walkers".format(nwalkers))
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[wf_exp])

# "burn-in" period; save final positions and then reset
pos, prob, state = sampler.run_mcmc(starting_guesses, nburn, progress=True)
sampler.reset()

# sampling period
sampler.run_mcmc(pos, nsteps, progress=True)
print("Mean acceptance fraction: {0:.3f} (in total {1} steps)"
      .format(np.mean(sampler.acceptance_fraction), nwalkers * nsteps))

# discard burn-in points and flatten the walkers; the shape of samples is (nwalkers*nsteps, ndim)
samples = sampler.chain.reshape((-1, ndim))

print("time elapsed: {:.2f}s".format(time.time() - start_time))

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
