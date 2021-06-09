import numpy as np
import corner
import matplotlib.pyplot as plt
from math import log

bayes_out = np.load(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\parallelout\bayes_out.npz')
bayes_out_prep = np.load(r'C:\Users\noeld\eclipse-workspace\get_wspot_schro\src\parallelout\bayes_out_prep.npz')
samples = bayes_out["samples"]
chain = bayes_out["chain"]
get_log_prob = bayes_out["log_prob"]
azaz = bayes_out["azaz"]
nl2j = bayes_out["nl2j"]
sigma_wf = bayes_out["sigma_wf"]
sigma_e = bayes_out["sigma_e"]
truths = bayes_out_prep["truths"]
index = bayes_out["index"]
autocorr = bayes_out["autocorr"]
#print(bayes_out_prep["e_br"])


N = 100

(nwalkers, nsteps, ndim) = chain.shape

if nl2j[2] == 1:
    truths = [-41.80, 3.18, 0.85, 2 * 2.45, 1.35, 0.17]
elif nl2j[2] == 3:
    truths = [-69.55, 1.89, 1.17, 2 * 2.13, 2.36, 0.21]
else:
    truths = bayes_out["truths"]

print(bayes_out_prep["truths"] - truths)

# make a corner plot with the posterior distribution
fig1, axs = plt.subplots(ndim, ndim)
# fig1.canvas.set_window_title('Correlation Plots')
# corner.corner(samples, labels=["Vws", "r0ws", "a0ws", "Vso", "r0so", "a0so"], show_titles=True, fig=fig1, levels=(0.68,))
# axs[0, ndim - 1].text(0.5, 0.8, 'A,Z,a,z = ' + str(azaz[0]) + ', ' + str(azaz[1]) + ',' + str(azaz[2]) + ', '
#                       + str(azaz[3]), horizontalalignment='center')
# axs[0, ndim - 1].text(0.5, 0.7, 'n,l,2j = ' + str(nl2j[0]) + ', ' + str(nl2j[1]) + ', ' + str(nl2j[2]),
#                       horizontalalignment='center')
# axs[0, ndim-1].text(0.5, 0.6, r'$V_0,a_0,r_0 = $' + str(truths[0]) + ', ' + str(truths[1]) + ', ' + str(truths[2]),
# horizontalalignment='center')
# axs[0, ndim - 1].text(0.5, 0.4, 'nsteps = ' + str(nsteps), horizontalalignment='center')
# axs[0, ndim - 1].text(0.5, 0.3, 'nwalkers = ' + str(nwalkers), horizontalalignment='center')
# axs[0, ndim - 1].text(0.5, 0.2, 'sigma = ' + str(sigma_wf), horizontalalignment='center')
# plt.savefig("corner.png")

# Plot chains for each varable
# TODO update for 6 variables
# fig2, axs = plt.subplots(ndim, sharex=True)
# for j in range(ndim):
#     for i in range(nwalkers):
#         axs[j].plot(range(nsteps), chain[i, :, j], linewidth=0.5)
#
# # axs[0].plot(range(nsteps), np.full((nsteps, 1), truths[0]), 'r.')
# axs[0].set(ylabel="Vws")
# # axs[1].plot(range(nsteps), np.full((nsteps, 1), truths[1]), 'r.')
# axs[1].set(ylabel="a0ws")
# # axs[2].plot(range(nsteps), np.full((nsteps, 1), truths[2]), 'r.')
# axs[2].set(ylabel="r0ws")
# axs[3].set(ylabel="Vso")
# axs[4].set(ylabel="a0so")
# axs[5].set(ylabel="r0so")
# plt.xlabel("Steps")

# Plot likelyhood
fig3 = plt.figure()
fig3.canvas.set_window_title('Likelihood Progression')
#plt.plot(range(nsteps), get_log_prob)
#plt.plot(range(nsteps), np.full((nsteps, 1), (-0.5 * N * np.log(2 * np.pi * sigma_wf ** 2) + (- 0.5 * N * np.log(2 * np.pi * sigma_e ** 2)))), 'k')
plt.xlabel("Steps")
plt.ylabel("Likelyhood")

# Compare to data
fig4 = plt.figure()
fig4.canvas.set_window_title('Wood-Saxon Fit')
plt.plot(bayes_out["wfr"], bayes_out_prep["wf"][:100], label='WS Fit', c='red', zorder=2)
if nl2j[2] == 1:
    plt.errorbar(bayes_out["wfr"], bayes_out["wf_exp"][:, 3], yerr=bayes_out["wf_exp"][:, 4], label='SA-NCSM', ms=5, c='gray', fmt='.', zorder=1)
    plt.plot(bayes_out["wfr"], bayes_out_prep["wf_br_p1"], label='Brita')
if nl2j[2] == 3:
    plt.errorbar(bayes_out["wfr"], bayes_out["wf_exp"][:, 1], yerr=bayes_out["wf_exp"][:, 2], label='SA-NCSM', ms=5, c='gray', fmt='.', zorder=1)
    plt.plot(bayes_out["wfr"], bayes_out_prep["wf_br_p3"][:100], label='Brita')
# plt.fill_between(bayes_out["wfr"], bayes_out_prep["wf_q1"], bayes_out_prep["wf_q3"], color='gray', alpha=0.2)
plt.xlabel("r (fm)")
plt.ylabel("Overlap ($fm^{-3/2}$)")
plt.legend()
plt.title("WS e = " + str(format(bayes_out_prep["e"], 'f')) + "\n" + "Exp e = " + str(bayes_out["e_exp"]))
# if nl2j[2] == 1:
#     plt.ylim((0, 0.15))
# else:
#     plt.ylim((0, 0.2))
# plt.xlim((0, 6))
# plt.savefig("wf.png")

# sf = integrate.simps(np.square(bayes_out["wf_exp"]) * np.square(bayes_out["wfr"]), bayes_out["wfr"])
# print(sf)
#
# norm = integrate.simps(np.square(bayes_out_prep["wf"]) * np.square(bayes_out["wfr"]), bayes_out["wfr"])
# print(norm)

# Plot autocorrelation time
fig5 = plt.figure()
fig5.canvas.set_window_title("Autocorrelation")
n = 100 * np.arange(1, index + 1)
y = autocorr[:index]
plt.plot(n, n / 100.0, "--k")
plt.plot(n, y)
# /plt.xlim(0, n.max())
plt.xlabel("number of steps")
plt.ylabel(r"mean $\hat{\tau}$")

# Plot Differentials between fit and exp data
fig6 = plt.figure()
fig6.canvas.set_window_title("Differential")
diff = bayes_out["wfr"]
if nl2j[2] == 1:
    diff = (bayes_out_prep["wf"] - bayes_out["wf_exp"][:, 4]) / bayes_out["wf_exp"][:, 3]
if nl2j[2] == 3:
    diff = (bayes_out_prep["wf"][:100] - bayes_out["wf_exp"][:100, 1]) / bayes_out["wf_exp"][:, 2]
plt.plot(bayes_out["wfr"], diff, label='Diff')
plt.title('Fit Differential')
plt.ylabel('Deviations from fit')
plt.xlabel('r (fm)')
plt.figtext(0.5, 0.01, '*Negative values indicate underfitting', ha='center', fontsize=8)

# Compare to data
fig7 = plt.figure()
fig7.canvas.set_window_title('Wood-Saxon Fit: Log Plot')
plt.plot(bayes_out["wfr"], np.log(bayes_out_prep["wf"]), label='WS Fit', c='red', zorder=2)
if nl2j[2] == 1:
    plt.errorbar(bayes_out["wfr"], np.log(bayes_out["wf_exp"][:, 3]), yerr=bayes_out["wf_exp"][:, 4], label='SA-NCSM', ms=5, c='gray', fmt='.', zorder=1)
    # plt.plot(bayes_out["wfr"], np.log(bayes_out_prep["wf_br_p1"]), label='Brita')
if nl2j[2] == 3:
    plt.errorbar(bayes_out["wfr"], np.log(bayes_out["wf_exp"][:, 1]), yerr=bayes_out["wf_exp"][:, 2], label='SA-NCSM', ms=5, c='gray', fmt='.', zorder=1)
    # plt.plot(bayes_out["wfr"], np.log(bayes_out_prep["wf_br_p3"]), label='Brita')
plt.fill_between(bayes_out["wfr"], np.log(bayes_out_prep["wf_q1"][:100]), np.log(bayes_out_prep["wf_q3"][:100]), color='gray', alpha=0.2)
plt.xlabel("r (fm)")
plt.ylabel("Log Overlap ($fm^{-3/2}$)")
plt.legend()
plt.title("WS e = " + str(format(bayes_out_prep["e"], 'f')) + "\n" + "Exp e = " + str(bayes_out["e_exp"]))

plt.show()
