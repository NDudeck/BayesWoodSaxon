import time
import numpy as np
from multiprocessing import Pool
import emcee


def log_prob(theta):
    t = time.time() + np.random.uniform(0.005, 0.008)
    while True:
        if time.time() >= t:
            break
    return -0.5 * np.sum(theta ** 2)


def main():
    np.random.seed(42)
    initial = np.random.randn(32, 5)
    nwalkers, ndim = initial.shape
    nsteps = 100
    serial_time = 21.1

    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
        start = time.time()
        sampler.run_mcmc(initial, nsteps, progress=True)
        end = time.time()
        multi_time = end - start
        print("Multiprocessing took {0:.1f} seconds".format(multi_time))
        print("{0:.1f} times faster than serial".format(serial_time / multi_time))


if __name__ == '__main__':
    main()
