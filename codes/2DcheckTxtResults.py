#!/usr/bin/env python
# -*- coding: utf-8 -*-
# from mpl_toolkits import mplot3d

from fipy import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection
from matplotlib import colors as mcolors
import scipy.sparse as sp;
import scipy.sparse.linalg as splin;
from scipy.optimize import newton
from scipy.optimize import fsolve
import os
import time

tested_param = 'pf'
typeOfS = 'homo'
a = 0.25
pf = 10.
delta = 10.
kr = 100.
kc = 100.
dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/with_growth_promotion/V_gompertz/with_Sr/parameter' + tested_param
# dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/without_growth_promotion/V_constant/with_Sr/parameter' + tested_param
# dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/deltaVariable/V_constant/with_Sr/parameter' + tested_param
filename = '2D' + typeOfS + 'SCParams' + tested_param + str(pf)\
           # + 'b1'
           # + 'g_mu1'
dataPath = dirpath + '/' + filename + ".txt"

data = np.loadtxt(open(dataPath,'rt').readlines()[:-1], skiprows=0)
# data = np.loadtxt(open(dataPath,'rt').readlines()[0:7010], skiprows=0)
data[:, 0]
Nt = data.shape[0]
# 0: time
# 1: dt time step
# 2: dz space step
# 3: mu_0
# 4: mu_1
# 5: mu1E
# 6: mu1B
# 7: Strength
# 8: I

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams.update({'font.size': 18})

# plt.plot(data[:, 0], data[:, 4] / np.max(data[:, 4]), label="$\mu_1$", linewidth=2)
# plt.plot(data[:, 0], data[:, 6] / np.max(data[:, 6]), label="$\mu_{c_r}$", linewidth=2)
# plt.plot(data[:, 0], data[:, 7] / np.max(data[:, 7]), '--', label="$\mu_c$", linewidth=1)
# plt.plot(data[:, 0], a * np.ones(Nt) / np.max(data[:, 7]), '--', label="$a$, div. rate", linewidth=1)
# plt.plot(data[:, 0], data[:, 8] / np.max(data[:, 8]), '--', label="$I$", linewidth=1)
#
# plt.legend()

truncate =30
n=2
size_of_fig = 8.
fig, ax1 = plt.subplots()
fig.set_size_inches(size_of_fig+2, size_of_fig)

# plt.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 5]/data[range(Nt-truncate), 6])

color1 = 'tab:red'
ax1.set_xlabel('t')
ax1.set_ylabel('$\mu_1$', color=color1)
ax1.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 4], '-', color=color1, linewidth=4, label='$t \mapsto \mu_1(t)$')
# ax1.legend()
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()

color2 = 'tab:blue'
color3 = 'tab:olive'
ax2.set_ylabel('$\overline{\mu_{c}}$', color=color2)
# ax2.plot(data[np.arange(0, Nt-truncate, 20), 0], data[np.arange(0, Nt-truncate, 20), 8],'-', marker='o', markersize=3, color=color2, linewidth=4, label='$t \mapsto I(t)$')
# ax2.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 6], '--', marker='o', markersize=2, color='blue', linewidth=4, label='$t \mapsto \mu_{c_r}$')
ax2.plot(data[np.arange(0, Nt-truncate, 1), 0], data[np.arange(0, Nt-truncate, 1), 7], '--', marker='o', markersize=2, color=color2, linewidth=4, label='$t \mapsto \overline{\mu_{c}}$')

ax2.plot(data[np.arange(0, Nt-truncate, 1), 0], data[-1, 7]* np.ones(np.arange(0, Nt-truncate, 1).shape[0]), '--', marker='o', markersize=0, color='k', linewidth=4, label='$\lambda$')

ax2.tick_params(labelcolor=color2)
ax1.legend(loc="upper left")
ax2.legend()
fig.tight_layout()


# truncate = 0
n=2
size_of_fig = 8.
fig, ax1 = plt.subplots()
fig.set_size_inches(size_of_fig+2, size_of_fig)
color1 = 'tab:red'
ax1.set_xlabel('t')
ax1.set_ylabel('$\mu_{c_r}$', color=color1)
ax1.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 6], '-', color=color1, linewidth=4, label='$t \mapsto \mu_{c_r}(t)$')
# ax1.legend()
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()

color2 = 'tab:blue'
color3 = 'tab:olive'
ax2.set_ylabel('$\mu_c$', color=color2)
# ax2.plot(data[:, 0], data[:, 8],'-', marker='o', markersize=3, color=color2, linewidth=4, label='$t \mapsto I(t)$')
ax2.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 7], '--', marker='o', markersize=2, color=color2, linewidth=4, label='$t \mapsto \mu_{c}$')

# ax2.plot(data[np.arange(0, Nt-truncate, 20), 0], a * np.ones(np.arange(0, Nt-truncate, 20).shape[0]), '--', marker='o', markersize=0, color='k', linewidth=4, label='$a$')

ax2.tick_params(labelcolor=color2)
ax1.legend(loc="upper right")
ax2.legend(loc="upper left")
fig.tight_layout()


# fig = plt.figure()
# plt.plot(data[:, 0], data[:, 5],  linewidth=4, label='$c$')
# plt.plot(data[:, 0], data[:, 6],  linewidth=4, label='$c_r$')
# plt.xlabel('$t$')
# plt.ylabel('$c_r, c$')

# fig = plt.figure()
# ax = plt.axes(projection="3d")
# ax.set_xlabel('$\mu_1$')
# ax.set_ylabel('$C$')
# ax.set_zlabel('$C_r$')
# ax.plot3D(data[:, 4], data[:, 7], data[:, 6])