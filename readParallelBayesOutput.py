import numpy as np
import corner
import matplotlib.pyplot as plt
import scipy.integrate as integrate

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
# truths = [-41.80, 0.85, 3.18, 2.45, 0.17, 1.35]
print(bayes_out_prep["e_br"])
print(bayes_out_prep["truths"] - [-41.80, 3.18, 0.85, 2*2.45, 1.35, 0.17])
N = 100

(nwalkers, nsteps, ndim) = chain.shape
# samples = samples[int(3*nsteps/4):]

# make a corner plot with the posterior distribution
fig1, axs = plt.subplots(ndim, ndim)
corner.corner(samples, labels=["Vws", "r0ws", "a0ws", "Vso", "r0so", "a0so"], show_titles=True, fig=fig1, levels=(0.68,), truths=[-41.80, 3.18, 0.85, 2*2.45, 1.35, 0.17])
#
# axs[1, 0].set(visible='off')
# axs[2, 0].set(visible='off')
# axs[3, 0].set(visible='off')
# axs[4, 0].set(visible='off')
#
# axs[2, 1].set(visible='off')
# axs[3, 1].set(visible='off')
# axs[4, 1].set(visible='off')
#
# axs[3, 2].set(visible='off')
# axs[4, 2].set(visible='off')
#
# axs[4, 3].set(visible='off')

axs[0, ndim - 1].text(0.5, 0.8, 'A,Z,a,z = ' + str(azaz[0]) + ', ' + str(azaz[1]) + ',' + str(azaz[2]) + ', '
                      + str(azaz[3]), horizontalalignment='center')
axs[0, ndim - 1].text(0.5, 0.7, 'n,l,2j = ' + str(nl2j[0]) + ', ' + str(nl2j[1]) + ', ' + str(nl2j[2]),
                      horizontalalignment='center')
# axs[0, ndim-1].text(0.5, 0.6, r'$V_0,a_0,r_0 = $' + str(truths[0]) + ', ' + str(truths[1]) + ', ' + str(truths[2]), horizontalalignment='center')

axs[0, ndim - 1].text(0.5, 0.4, 'nsteps = ' + str(nsteps), horizontalalignment='center')
axs[0, ndim - 1].text(0.5, 0.3, 'nwalkers = ' + str(nwalkers), horizontalalignment='center')
axs[0, ndim - 1].text(0.5, 0.2, 'sigma = ' + str(sigma_wf), horizontalalignment='center')
# plt.savefig("corner.png")


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
#
# fig3 = plt.figure()
# plt.plot(range(nsteps), get_log_prob)
# plt.plot(range(nsteps), np.full((nsteps, 1), (-0.5 * N * np.log(2 * np.pi * sigma_wf ** 2) + (- 0.5 * N * np.log(2 * np.pi * sigma_e ** 2)))), 'k')
# plt.xlabel("Steps")
# plt.ylabel("Likelyhood")


# Compare to data
fig4 = plt.figure()
plt.plot(bayes_out["wfr"], bayes_out_prep["wf"])
plt.plot(bayes_out["wfr"], bayes_out["wf_exp"])
plt.plot(bayes_out["wfr"], bayes_out_prep["wf_br"])
plt.fill_between(bayes_out["wfr"], bayes_out_prep["wf_q1"], bayes_out_prep["wf_q3"], color='gray', alpha=0.2)
plt.legend(["WS", "NCNSM", "Brida"])
plt.title("th e = " + str(bayes_out_prep["e"]) + "\n" + "exp e = " + str(bayes_out["e_exp"]))
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

plt.show()
