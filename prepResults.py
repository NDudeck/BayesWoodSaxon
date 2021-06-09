import get_wspot_schro as gws
import numpy as np
import corner
import math

bayes_out = np.load(r'/mnt/c/Users/noeld/eclipse-workspace/get_wspot_schro/src/parallelout/bayes_out.npz')
azaz = bayes_out["azaz"]
nl2j = bayes_out["nl2j"]
samples = bayes_out["samples"]
sf = bayes_out["sf"]
emin, emax = bayes_out["e"]


# Method for running wood-saxon executable and grabbing output wavefunction
def get_wf_ws(theta):
    wf, e, rms = gws.gwf(azaz[0], azaz[1], azaz[2], azaz[3], nl2j[0], nl2j[1], nl2j[2], theta[0], theta[1], theta[2],
                         theta[3], theta[4], theta[5], emin, emax)

    wf = wf * math.sqrt(sf)
    print(rms)

    return wf, e


theta_th = np.mean(samples, axis=0)
theta_brida_p1 = [-41.80, 3.18, 0.85, 2 * 2.45, 1.35, 0.17]
theta_brida_p3 = [-69.55, 1.89, 1.17, 2 * 2.13, 2.36, 0.21]

theta_q1 = corner.quantile(samples[:, 0], 0.25)
theta_q3 = corner.quantile(samples[:, 0], 0.75)

wf_th, e_th = get_wf_ws(theta_th)
theta_th[0] = theta_q1
wf_q1, e_q1 = get_wf_ws(theta_th)

theta_th[0] = theta_q3
wf_q3, e_q3 = get_wf_ws(theta_th)

theta_th[0] = np.mean(samples[:, 0])

wf_br_p1, e_br_p1 = get_wf_ws(theta_brida_p1)
wf_br_p3, e_br_p3 = get_wf_ws(theta_brida_p3)


np.savez_compressed(r'/mnt/c/Users/noeld/eclipse-workspace/get_wspot_schro/src/parallelout/bayes_out_prep.npz',
                    wf=wf_th, wf_q1=wf_q1, wf_q3=wf_q3, truths=theta_th, e=e_th, wf_br_p1=wf_br_p1, e_br_p1=e_br_p1,
                    wf_br_p3=wf_br_p3, e_br_p3=e_br_p3)
