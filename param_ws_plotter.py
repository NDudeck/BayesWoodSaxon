import numpy as np
import subprocess
import matplotlib.pyplot as plt
import time

# Load in wf_ws
A = 6  # Lithium 7
Z = 3
a = 1  # Neutron is 1,0
z = 0
nl2j = [0, 1, 1]

# Load Experimental data
wf_exp_in = np.genfromtxt("6Li_7Li_overlaps_NNLOopt_hw=10MeV_Nmax12.csv", skip_header=1, delimiter=",")
wf_exp = wf_exp_in[:, 1]

plt.ion()
figure, ax = plt.subplots(figsize=(8, 6))
line1, = ax.plot(wf_exp_in[:, 0], wf_exp)
line2, = ax.plot(wf_exp_in[:, 0], wf_exp)
plt.title("Dynamic Plot of psi", fontsize=25)
plt.xlabel("r", fontsize=18)
plt.ylabel("psi(r)", fontsize=18)

min_theta = np.array([0, 0, 0])
max_theta = np.array([2, 5, 5])


def get_wf_ws(theta):
    l7_dia = open("io/l7.dai", 'w')
    l7_dia.write(str(A) + " " + str(Z) + " " + str(a) + " " + str(z) + "\n\n")
    l7_dia.write(' '.join(map(str, nl2j)) + ' ' + ' '.join(map(str, theta)))
    l7_dia.close()

    subprocess.check_output("wspot_schro.exe l7")

    try:
        wf = np.genfromtxt("rad.dat")
        # wf[:, 1] = wf[:, 1] * wf[:, 0]
        wf = wf[:100, :]
        return wf
    except ValueError:
        return np.zeros((100, 2))


steps = 100
v = 1.01
a0 = 1.41
r = 1.01

outp = np.zeros((steps, 2))
index = 0

for r in np.linspace(1, 10, steps):
    wf_ws = get_wf_ws([v, a0, r])
    line1.set_xdata(wf_ws[:, 0])
    line1.set_ydata(wf_ws[:, 1])
    figure.canvas.draw()
    figure.canvas.flush_events()
    outp[index, 0] = r
    outp[index, 1] = np.mean(np.abs(wf_ws[:, 1] - wf_exp))
    index = index + 1

print('Min r = ' + str(outp[outp[:, 1] == min(outp[:, 1]), 0]))
