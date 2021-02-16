import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner


# Method for running wood-saxon executable and grabbing output wavefunction
def get_wf_ws(theta):
    # Write file input with given parameters
    l7_dia = open("l7.dai", 'w')
    l7_dia.write(str(A) + " " + str(Z) + " " + str(a) + " " + str(z) + "\n\n")
    l7_dia.write(' '.join(map(str, nl2j)) + ' ' + ' '.join(map(str, theta)))
    l7_dia.close()

    # Call executable without outputting into console
    subprocess.check_output("wspot_schro.exe l7")

    try:
        wf = np.genfromtxt("rad.dat").astype(float)  # Grab output
        wf = wf[:100, :]  # Trim to match input
        return wf
    except ValueError:
        # If parameters aren't good, return flat wavefunction
        return np.zeros((100, 2))


# Start run timer
start_time = time.time()

# Initialize physical parameters
A = 6  # Lithium 6
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 1, 1]  # Desired quantum numbers

# Load in sampling parameters
ndim = 3  # number of parameters in the model
nwalkers = 24  # number of MCMC walkers
nburn = 50  # "burn-in" period to let chains stabilize
nsteps = 400  # number of MCMC steps to take
sigma = 0.02

# Load Experimental data
# wf_exp_in = np.genfromtxt("6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12.csv", skip_header=1, delimiter=",")
wf_exp_in = get_wf_ws([1, 1.25, 0.65])
wf_exp = np.square(wf_exp_in[:, 1] * wf_exp_in[:, 0] + np.random.normal(0, 0.02, 100))

# Plot experimental data
plt.ion()
figure, ax = plt.subplots(figsize=(8, 6))
line1, = ax.plot(wf_exp_in[:, 0], wf_exp)
line2 = ax.scatter(wf_exp_in[:, 0], wf_exp)
plt.title("WS Psi vs Exp Psi", fontsize=25)
plt.xlabel("r", fontsize=18)
plt.ylabel("psi^2(r)", fontsize=18)

# Build parameter space theta [V_N, a_0, r_0]
min_theta = np.array([0.5, 0.5, 0.25])
max_theta = np.array([1.5, 2, 1.25])
volume_theta = np.prod(max_theta - min_theta)


def log_prior(theta):
    # Flat prior
    if np.logical_and(min_theta <= theta, theta <= max_theta).all():
        return np.log(1 / volume_theta)
    else:
        return -np.inf


def log_likelihood(theta, wf_exp):
    # print("Params " + str(theta))

    wf_ws = get_wf_ws(theta)  # Grab wavefunction with provided theta

    wf_ws = np.square(wf_ws[:, 1] * wf_ws[:, 0])  # turn phi(r)/r to psi(r) then square
    # Update dynamic plot
    line1.set_xdata(np.linspace(0, 10, 100))
    line1.set_ydata(wf_ws)
    figure.canvas.draw()
    figure.canvas.flush_events()
    # time.sleep(0.1)

    try:
        # Chi squared log likelyhood
        (N,) = np.shape(wf_ws)
        diff = wf_exp - wf_ws
        ll = -0.5 * np.sum(((diff) / sigma) ** 2) + 0.5 * N * np.log(2 * np.pi * sigma ** 2)
        print(ll)
        return ll
    except ValueError:
        print("Value Error")
        return -np.inf


def log_posterior(theta, wf_exp):
    return log_prior(theta) + log_likelihood(theta, wf_exp)


# Start at random locations within the prior volume
starting_guesses = np.random.uniform(min_theta, max_theta, (nwalkers, ndim))

# Initialize sampler using paramters
print("MCMC sampling using emcee (affine-invariant ensamble sampler) with {0} walkers".format(nwalkers))
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior, args=[wf_exp])

# "burn-in" period; save final positions and then reset
pos, prob, state = sampler.run_mcmc(starting_guesses, nburn)
sampler.reset()

# sampling period
sampler.run_mcmc(pos, nsteps)
print("Mean acceptance fraction: {0:.3f} (in total {1} steps)"
      .format(np.mean(sampler.acceptance_fraction), nwalkers * nsteps))

# discard burn-in points and flatten the walkers; the shape of samples is (nwalkers*nsteps, ndim)
samples = sampler.chain.reshape((-1, ndim))

print("time elapsed: {:.2f}s".format(time.time() - start_time))

# make a corner plot with the posterior distribution
fig1 = corner.corner(samples, labels=["V", "a", 'r'], truths=[1, 1.25, 0.65], show_titles=True)

# fig2, axs = plt.subplots(3)
# for j in range(ndim):
#     for i in range(nwalkers):
#         axs[j].plot(range(nsteps), sampler.chain[i, :, j], linewidth=0.5)
