#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from scipy import optimize as opt
from scipy import stats as st
from matplotlib import pyplot as plt
import pandas as pd
from sklearn.linear_model import LinearRegression as linreg

# ------------------------------------------------------------------------------------------------------------------
# ------------------------------ Estimation of the normal influx rate of Im. cells R  ------------------------------
# ------------------------------------------------------------------------------------------------------------------

data_Tcells = np.genfromtxt('tumor_volume_vs_Im_cells_rate.csv', delimiter=',')
x_Tcells_tmp = data_Tcells[:, 0]*1e9
y_Tcells_tmp = data_Tcells[:, 1]*1e6

# Detection of the outliers using zscore
idx_outliers_x = np.where(np.abs(st.zscore(x_Tcells_tmp)) <= 1.8)[0]
idx_outliers_y = np.where(np.abs(st.zscore(y_Tcells_tmp)) <= 2.)[0]

idxs = []
for y in idx_outliers_y:
    if y in idx_outliers_x and y in idx_outliers_y:
        idxs.append(y)

x_Tcells = x_Tcells_tmp[idxs]
y_Tcells = y_Tcells_tmp[idxs]


def f_affine(x, p, q):
    return p*x + q

def f_lin(x, p):
    return p*x


def f_poly(x, b, c):
     return b*x**2 + c*x

# --------------------------------- Linear regression without intercept --------------------------------
popt, pcov = opt.curve_fit(f_lin, x_Tcells, y_Tcells, method='lm', absolute_sigma=True)
slope = np.mean(x_Tcells*y_Tcells)/np.mean(x_Tcells**2)

# --------------------------------- Computation of the Confidence interval ------------------------------------
n = y_Tcells.shape[0]
sig_chap = np.sum((y_Tcells - f_lin(x_Tcells, slope))**2)/(n-1)

SCE = np.sum((f_lin(x_Tcells, slope))**2)
SCT = np.sum(( y_Tcells)**2)
R_square  = SCE/SCT

# 95% confidence interval
t_alpha = st.t.ppf(0.975, n-1)
Sxx = np.sum((x_Tcells)**2)
slope_int_size = t_alpha * np.sqrt(sig_chap/Sxx)


#  The confidence interval
conf_inter = [slope - slope_int_size, slope + slope_int_size]
print('confidence interval: ',conf_inter)

# std err
std_err = np.sqrt(sig_chap/Sxx)
print('the std err: ',std_err)

#  The pvalue
test_stat = slope/(np.sqrt(sig_chap/Sxx))
pvalue = 2*(1 - st.t.cdf(test_stat, n-1))


# -----------------------Visualization of the regression result --------------------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
# plt.rcParams["text.usetex"] = True
plt.rcParams.update({'font.size': 20})

size_of_fig=10
fig = plt.figure(2)
fig.set_size_inches(size_of_fig+5, size_of_fig)
plt.scatter(x_Tcells/1e11, y_Tcells, s=np.asarray([400.]), c='k')
plt.plot(x_Tcells/1e11, f_lin(x_Tcells, slope), 'k-', label='slope=%5.9f'% tuple(np.array([popt[0]])), linewidth=4)
# plt.plot(x_Tcells, f_poly(x_Tcells, *popt))
plt.xlabel('$\mu_1$ in $\mu m^3$ ( $\\times$ $10^{11}$)')
plt.ylabel('$Y$ in $cell_c \cdot day^{-1}$')
plt.legend()


# ------------------------------------------------------------------------------------------------------------------
# ---------------- Estimation of the tumor growth parameters (the division rate a and growth rate V)----------------
# ------------------------------------------------------------------------------------------------------------------

data = pd.read_csv('tumor_time_to_event_data.csv')
# data = pd.read_csv('Vero_data_tumor_volume_Allmice_means.csv')

# function for retreiving each tumor evolution data from each cinetic
def get_tumor_by_cinetique_id(data, cin, tum_id):
    tumors = []
    for i in range(data.shape[0]):
        if data.get('A')[i] == cin and data.get('B')[i]==tum_id:
            tumors.append([data.get('C')[i], data.get('D')[i], data.get('E')[i]])

    return np.asarray(tumors)

# Different kinetics: C2, C3, C4
# Tumors in each kinetic: C2 (T1, T2, T3, T4, T5) -- C3 (T1, T2, T3) -- C4 (T1, T2, T3, T4, T5)


tumor_data = get_tumor_by_cinetique_id(data, 'C4', 'T4')

time_mu1s = tumor_data[0:3, 0]
mu1s = tumor_data[0:3, 1]*1e9

time_mu0s = tumor_data[0:3, 0]
mu0s = tumor_data[0:3, 2]*1e6

scale = 1e8
size_of_fig_2 = 6

# ----------------------------------------------------------------------

mu1s_mean = mu1s.mean()
mu1s_std = mu1s.std()

mu0s_mean = mu0s.mean()
mu0s_std = mu0s.std()
# mu0s = np.concatenate([tumor_data[0:2, 2], [tumor_data[3, 2]]])


def f_mu0s(t, a):
    return mu0s[0]*np.exp(a*(t - time_mu0s[0]))


def f_mu1s(t, V):
    return mu1s[0] + (mu0s[0]*V/a)*(np.exp(a * (t - time_mu1s[0])) - 1.)


# def obj_func(x):
#     f = 0.
#     for i in range(mu0s.shape[0]):
#         f = f + ((f_mu0s(xdata_mu0s[i], x[0]) - mu0s[i])**2) + (f_mu1s(xdata_mu0s[i], x[1]) - mu1s[i])**2
#     return f
#
#
# es = cma.CMAEvolutionStrategy(2*[1.], 0.5)
# es.optimize(obj_func)
method = 'lm'
popt_mu0, pcov_mu0 = opt.curve_fit(f_mu0s, time_mu0s, mu0s ,method=method, absolute_sigma=True)
a = popt_mu0[0]

scale_mu0 = 1e8
scale_mu1 = 1e11

fig = plt.figure(1)
fig.set_size_inches(size_of_fig_2 + 3, size_of_fig_2)
plt.plot(time_mu0s, mu0s/scale_mu0, 'k-', label='data')
plt.scatter(time_mu0s, mu0s/scale_mu0, s=np.asarray([100.]),  c='k')
plt.plot(time_mu0s, f_mu0s(time_mu0s, a)/scale_mu0, 'k--', label='fit LS: a=%5.3f $day^{-1}$' % tuple(np.array([a])), linewidth=3)
plt.xlabel('$t$ in $day$')
plt.ylabel('$\mu_0$ and $Y^{(0)}$ ( $\\times$ $10^{8}$)')
plt.legend()


popt_mu1, pcov_mu1 = opt.curve_fit(f_mu1s, time_mu1s, mu1s, method=method, absolute_sigma=True)
V = popt_mu1[0]

fig = plt.figure(2)
fig.set_size_inches(size_of_fig_2 + 3, size_of_fig_2 )
plt.plot(time_mu1s, mu1s / scale_mu1, 'k-', label='data')
plt.scatter(time_mu1s, mu1s / scale_mu1, s=np.asarray([100.]), c='k')
plt.plot(time_mu1s, f_mu1s(time_mu1s, V) / scale_mu1, 'k--', label='fit LS: V=%5.3f $\mu m^3/day$'
                                                                   % tuple(np.array([V])), linewidth=3)

plt.xlabel('$t$ in $day$')
plt.ylabel('$\mu_1$ and $Y^{(1)}$ ( $\\times$ $10^{11}$)')
plt.legend()

