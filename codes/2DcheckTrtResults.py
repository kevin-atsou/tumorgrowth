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

# T1 = 1.
# T2 = 7.
# t0 = 0.
# da = 0.002
# di = 0.

das = [0.002, 0.5, 1., 2., 3., 5., 7., 10.]
dis = [0.002, 0.04, 0.08, 0.12, 0.16, 0.2, 0.24, 0.28, 0.32, 0.36, 0.4]
t0s = [0., 5., 10., 15.]
Tf = 45.
# typeOfTreat = 'Anergic'
typeOfTreat = 'anergic'
# typeOfTreat = 'cytokine'
typeOfS = 'homo'
# typeOfTreat = 'Cytokine'

# Take the value without treatment

tested_param = 'a'
a = 4.
dirpath_w = 'TestResults_Protumor/' + typeOfS + 'geneousS/Treatments/V_gompertz/without_treatment/'
filename_w = '2D' + typeOfS + 'SCParams' + tested_param + str(a)
dataPath_w = dirpath_w + '/' + filename_w + ".txt"

data_w = np.loadtxt(open(dataPath_w, 'rt').readlines()[:-1], skiprows=0)
values_w = np.where(np.logical_and(data_w[:, 0] >= Tf, data_w[:, 0] <= Tf + 1e-2))
mu1_w = data_w[values_w[0][0], 4]

del data_w, values_w
# data = np.loadtxt(open(dataPath,'rt').readlines()[0:7010], skiprows=0)

# Nt = data.shape[0]
Effs = []
if typeOfTreat == "cytokine":
    for di in dis:
        Effs_tmp = []
        for t0 in t0s:
            dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/Treatments/V_gompertz/' + typeOfTreat + '_treatment/treatment_t0_' + str(
                t0)
            filename = '2D' + typeOfTreat + '_treat_' + 'di_' + str(di)
            dataPath = dirpath + '/' + filename + ".txt"
            data = np.loadtxt(open(dataPath, 'rt').readlines()[:-1], skiprows=0)

            values = np.where(np.logical_and(data[:, 0] >= Tf, data[:, 0] <= Tf + 1e-2))

            Eff = (mu1_w - data[values[0][0], 4])/mu1_w
            # Eff = data[values[0][0], 4]
            # Eff = data[values[0][0], 4]/mu1_w

            Effs_tmp.append(Eff)

        Effs.append(Effs_tmp)
elif typeOfTreat == "anergic":
    for da in das:
        Effs_tmp = []
        for t0 in t0s:
            dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/Treatments/V_gompertz/' + typeOfTreat + '_treatment/treatment_t0_' + str(
                t0)
            filename = '2D' + typeOfTreat + '_treat_' + 'da_' + str(da)
            dataPath = dirpath + '/' + filename + ".txt"
            data = np.loadtxt(open(dataPath, 'rt').readlines()[:-1], skiprows=0)

            values = np.where(np.logical_and(data[:, 0] >= Tf, data[:, 0] <= Tf + 1e-2))

            Eff = (mu1_w - data[values[0][0], 4]) / mu1_w
            # Eff = data[values[0][0], 4]
            # Eff = data[values[0][0], 4]/mu1_w

            Effs_tmp.append(Eff)

        Effs.append(Effs_tmp)

del data, values
Effs = np.asarray(Effs)
a, b, c, d = 0., 15., 0.002, 10. 
# plt.imshow(Effs, origin='lower', aspect='auto', extent=[a, b, c, d], cmap='Spectral_r')
figure, axes = plt.subplots()
caxes = axes.matshow(Effs, cmap='Spectral_r', interpolation ='nearest')
figure.colorbar(caxes)

if typeOfTreat == "cytokine":
    axes.set_xticklabels(['']+t0s)
    axes.set_yticklabels(['']+dis)

    size_of_fig = 8.

    fig, ax1 = plt.subplots()
    fig.set_size_inches(size_of_fig+2, size_of_fig)

    ax1.plot(t0s, Effs[1], 'r--', linewidth=4, label='$q ='+str(dis[1]) +'$');
    ax1.legend(); ax1.set_xlabel('$t_0$'); ax1.set_ylabel('$E$')
    plt.plot(t0s, Effs[2], 'k-', linewidth=4, label='$q ='+str(dis[2]) +'$');
    ax1.legend(); ax1.set_xlabel('$t_0$'); ax1.set_ylabel('$E$')
    ax1.plot(t0s, Effs[3], '-', color='tab:blue', linewidth=4, label='$q ='+str(dis[3]) +'$');
    ax1.legend(); ax1.set_xlabel('$t_0$'); ax1.set_ylabel('$E$')
    plt.axvline(x=t0s[1], color='gray', linestyle='--')
    plt.axvline(x=t0s[2], color='gray', linestyle='--')
    plt.axvline(x=t0s[3], color='gray', linestyle='--')

    # ax2 = ax1.twinx()
    #
    # ax2.plot(t0s, Effs[0], 'b:', linewidth=4, label='$q ='+str(dis[0]) +'$');
    # ax2.legend(); ax2.set_xlabel('$t_0$'); ax2.set_ylabel('$E$')
    # ax2.tick_params(labelcolor='blue')
elif typeOfTreat == "anergic":
    axes.set_xticklabels([''] + t0s)
    axes.set_yticklabels([''] + das)

    size_of_fig = 8.

    fig, ax1 = plt.subplots()
    fig.set_size_inches(size_of_fig + 2, size_of_fig)

    ax1.plot(t0s, Effs[1], 'r--', linewidth=4, label='$q =' + str(das[1]) + '$');
    ax1.legend();
    ax1.set_xlabel('$t_0$');
    ax1.set_ylabel('$E$')
    plt.plot(t0s, Effs[2], 'k-', linewidth=4, label='$q =' + str(das[2]) + '$');
    ax1.legend();
    ax1.set_xlabel('$t_0$');
    ax1.set_ylabel('$E$')
    ax1.plot(t0s, Effs[3], '-', color='tab:blue', linewidth=4, label='$q =' + str(das[3]) + '$');
    ax1.legend();
    ax1.set_xlabel('$t_0$');
    ax1.set_ylabel('$E$')
    plt.axvline(x=t0s[1], color='gray', linestyle='--')
    plt.axvline(x=t0s[2], color='gray', linestyle='--')
    plt.axvline(x=t0s[3], color='gray', linestyle='--')

    # ax2 = ax1.twinx()
    #
    # ax2.plot(t0s, Effs[0], 'b:', linewidth=4, label='$q =' + str(das[0]) + '$');
    # ax2.legend();
    # ax2.set_xlabel('$t_0$');
    # ax2.set_ylabel('$E$')
    # ax2.tick_params(labelcolor='blue')

dec_a = 0.05 # decay of the treatment in the microenvironment
dec_i = 0.0105# decay of the treatment in the microenvironment

T1 = 1.
T2 = 7.
t0 = 5.
da = 0.5
di = 0.28

N_opt = (np.log(1. - (dec_i*(1. - np.exp(dec_i*T2))/(di*np.exp(dec_i*t0)*(np.exp(dec_i*T1) - 1.))))/(dec_i*T2)) - 1.

Tf = 45.

dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/Treatments/V_gompertz/'+typeOfTreat+'_treatment/treatment_t0_'+str(t0)

if typeOfTreat == "cytokine":
    filename = '2D' + typeOfTreat + '_treat_'+'di_'+str(di)
elif typeOfTreat == "anergic":
    filename = '2D' + typeOfTreat + '_treat_' + 'da_' + str(da)
elif typeOfTreat == "anergiccytokine":
    filename = '2D' + typeOfTreat + '_treat_' + 'da_' + str(da)+'_'+'di_'+str(di)

dataPath = dirpath + '/' + filename + ".txt"

data = np.loadtxt(open(dataPath, 'rt').readlines()[:-1], skiprows=0)
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
# 9: Trt_a
# 10: mu1_a

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
truncate = 0
n=2
fig, ax1 = plt.subplots()
fig.set_size_inches(size_of_fig+2, size_of_fig)
color1 = 'tab:red'
ax1.set_xlabel('t')
ax1.set_ylabel('$\mu_1$', color=color1)
ax1.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 4], '-', color=color1, linewidth=4, label='$t \mapsto \mu_1(t)$')

ax1.axvline(x=t0, color='k', linestyle='-.', linewidth=4)
# ax1.legend()
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()

color2 = 'tab:blue'
color3 = 'tab:olive'
ax2.set_ylabel('$\mu_{c}$', color=color2)
# ax2.plot(data[np.arange(0, Nt-truncate, 20), 0], data[np.arange(0, Nt-truncate, 20), 8],'-', marker='o', markersize=3, color=color2, linewidth=4, label='$t \mapsto I(t)$')
# ax2.plot(data[range(Nt-truncate), 0], data[range(Nt-truncate), 6], '--', marker='o', markersize=2, color='blue', linewidth=4, label='$t \mapsto \mu_{c_r}$')
# ax2.plot(data[np.arange(0, Nt-truncate, 1), 0], data[np.arange(0, Nt-truncate, 1), 10], '--', marker='o', markersize=2, color=color2, linewidth=4, label='$t \mapsto \overline{\mu_{c}}$')
ax2.plot(data[np.arange(0, Nt-truncate, 1), 0], data[np.arange(0, Nt-truncate, 1), 7], '--', marker='o', markersize=2, color=color2, linewidth=4, label='$t \mapsto \overline{\mu_{c}}$')

# ax2.plot(data[np.arange(0, Nt-truncate, 1), 0], a * np.ones(np.arange(0, Nt-truncate, 1).shape[0]), '--', marker='o', markersize=0, color='k', linewidth=4, label='$a$')

# ns = np.where(np.logical_and(data[:, 0] >= t0, data[:, 0] <= t0 + 1e-2))
# if t0 == 0.:
#     ax2.annotate('treatment', xy=(t0, data[ns[0][0], 7]), xytext=(t0, data[ns[0][0], 7] + 0.25),
#                  arrowprops=dict(facecolor='black', shrink=0.01))
# elif t0> 0.:
#     ax2.annotate('treatment', xy=(t0, data[ns[0][0], 7]), xytext=(t0 - 2., data[ns[0][0], 7] + 0.25),
#                 arrowprops=dict(facecolor='black', shrink=0.05))
ax2.tick_params(labelcolor=color2)
# ax1.legend(loc="upper left")
# ax2.legend()
fig.tight_layout()
