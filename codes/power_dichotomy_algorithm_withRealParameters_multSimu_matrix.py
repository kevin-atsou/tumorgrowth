##!/usr/bin/env python
# -*- coding: utf-8 -*-

from fipy import *
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
import scipy.linalg as lin
from scipy import signal as sign
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import sys
import os
import time
import parameterFunctions.divisionRate as ad
import parameterFunctions.immuneResponse as delt
import parameterFunctions.gammaF as gammaF
import parameterFunctions.pF as pF
import parameterFunctions.sigmaF as sigmaF
import parameterFunctions.divisionOperatorKernel as K
import parameterFunctions.TestFeature as TF
# import pint as pt

# u = pt.UnitRegistry()

def plus(z):
    return 0.5*(z+np.abs(z))


def minus(z):
    return 0.5*(-z+np.abs(z))

# Stationary profile


def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)


def dichotomy(mu_a, mu_b, eps, mesK, Q, aDiffE, aConVE, Id, gF, t_s, E, delta_tild):
    while mu_b - mu_a > eps:
        mu = (mu_b + mu_a) / 2.

        bE = (mesK * Q * mu).T
        E_new = la.spsolve((- aDiffE.tocsc() - mu * aConVE + Id.multiply(mesK * gF * t_s)), bE)
        E.setValue(E_new)
        E.updateOld()

        F_mu_m = numerix.sum(mesK * delta_tild.value * E.value) - 1.

        bE = (mesK * Q * mu_a).T
        E_new = la.spsolve((- aDiffE.tocsc() - mu_a * aConVE + Id.multiply(mesK * gF * t_s)), bE)
        E.setValue(E_new)
        E.updateOld()

        F_mu_a = numerix.sum(mesK * delta_tild.value * E.value) - 1.

        # print ('{0} x {1}'.format(F_mu_m, F_mu_a))

        if F_mu_m * F_mu_a <= 0:
            mu_b = mu
        else:
            mu_a = mu
    return mu


def not_converge(x, y):
    if  (np.abs(x-y)/y)<= 1e-4:
        return True
    else:
        return False


# xt, yt = mesh.cellCenters
# nVol = mesh.numberOfCells

# ----------------------Configuration the test cases and the results saving --------------------

secs = 1.
# a = 4.
a = 0.1/ secs # / u.day (the tumor cells division rate)
# V = 713.608 / secs  # * (u.um**3)/u.day
V = np.float64(308.526) / secs  # * (u.um**3)/u.day (the tumor cells growth rate)
delta = np.float64(2.) / secs  # / u.day

pf = 0.25 # This parameter is no longer used in the case with real parameters

# S = np.float64(7.573e-8)/secs # /(u.day * (u.mm**3) * u.micrometer**3)
S = np.float64(6.11e-7)/secs # /(u.day * (u.mm**3) * u.micrometer**3) IC4 (the influx rate R)
# S = np.float64(1.231e-6)/secs # /(u.day * (u.mm**3) * u.micrometer**3)
chi = np.float64(8.64e1)/secs  # * (u.mm**2)/(u.mmol * u.day) (the chemotactic coefficient \chi)
D = np.float64(8.64e-5) / secs # * (u.mm**2)/(u.day) (Space diffusion of Im. Cells D)
# D = np.float64(8.64e-6) / secs # * (u.mm**2)/(u.day)
# D = np.float64(2.5e-2) / secs # * (u.mm**2)/(u.day)
sF = np.float64(5e-17) /secs # * u.mmol/(u.day * u.micrometer**3) (the strength of the chemical signal A_{\sigma})
gF = np.float64(2e-2)/secs  # /u.day (the death rate of Im. cells \gamma)
d = np.float64(1e-2)/secs  # *(u.mm**2)/(u.day) (the diffusion rate of the attractive potential K)
init_value = 5.
alpha = 0.
beta = 0.
tested_param_name = 'a'

a_max =  np.float64(0.351)/ secs # / u.day
V_max = np.float64(2521.975) / secs  # * (u.um**3)/u.day
delta_max = np.float64(60.) / secs  # / u.day
pf_max = 0.25
# S_max = np.float64(1.520e-6)/secs # /(u.day * (u.mm**3) * u.micrometer**3)
S_max = np.float64(9.74e-7)/secs # /(u.day * (u.mm**3) * u.micrometer**3)
chi_max = np.float64(8.64e6)/secs  # * (u.mm**2)/(u.mmol * u.day)
D_max = np.float64(1e-3) / secs # * (u.mm**2)/(u.day)
sF_max = np.float64(0.625e-6) /secs # * u.mmol/(u.day * u.micrometer**3)
gF_max = np.float64(1.)/secs  # /u.day
d_max = np.float64(1.)/secs  # *(u.mm**2)/(u.day)

n_step = 100.
a_s = np.arange(a, a_max, (a_max - a)/n_step)
# V_s = np.arange(308.526, 2521.975, step)
delta_s = np.arange(delta, delta_max, (delta_max-delta)/n_step)
S_s = np.arange(S, S_max, (S_max - S)/n_step)
chi_s = np.arange(chi, chi_max, (chi_max - chi)/n_step)
D_s = np.arange(D, D_max, (D_max - D)/n_step)
sF_s = np.arange(sF, sF_max, (sF_max - sF)/n_step)
gF_s = np.arange(gF, gF_max, (gF_max - gF)/n_step)
K_s = np.arange(d, d_max, (d_max-d)/n_step)

t_r = 0.


params = TF.TestFeature(a, V, delta, pf, S, chi, D, sF, gF, d, tested_param_name, init_value, alpha, beta)

mus_mat = []
n = 1

# tested_vars = a_s

for a in a_s:
    mus = []
    # for S in S_s[[0, 1, 5, 10, 20]]:
    for gF in gF_s[[0, 24, 48, 72, 96]]:
        t_s = 1./a
        x_s = np.sqrt(D*t_s)
        c_s = 1./(t_s*(x_s**2)*delta)
        # nu = (D/a)*np.sqrt(D/a)*(S*d*delta)/(chi*sF)
        # mu1_s = mu1_tild
        mu1_s = c_s/(S*t_s)
        mu0_s = a*mu1_s/V

        phi_s = (x_s**2)/(mu1_s*t_s*chi)

        # Q = S*mu1_s*t_s/c_s
        Q = 1.
        U = (sF*x_s**2)/(d*phi_s)

        radius = 1./x_s
        cellSize = radius/10.

        mesh = Gmsh2D('''
                       cellSize = %(cellSize)g;
                       radius = %(radius)g;
                       Point(1) = {0, 0, 0, cellSize};
                       Point(2) = {-radius, 0, 0, cellSize};
                       Point(3) = {0, radius, 0, cellSize};
                       Point(4) = {radius, 0, 0, cellSize};
                       Point(5) = {0, -radius, 0, cellSize};
                       Circle(6) = {2, 1, 3};
                       Circle(7) = {3, 1, 4};
                       Circle(8) = {4, 1, 5};
                       Circle(9) = {5, 1, 2};
                       Line Loop(10) = {6, 7, 8, 9};
                       Plane Surface(11) = {10};
                       ''' % locals())

        # c_s = t_s*S*mu1_s
        # Diffusion Matrix construction
        x = mesh.cellCenters
        xt, yt = mesh.cellCenters

        nVol = mesh.numberOfCells
        nFaces = mesh.numberOfFaces

        intF = mesh.interiorFaceIDs
        extF = np.arange(0, nFaces, 1)[np.array(mesh.exteriorFaces)]

        intFacesCells = mesh.faceCellIDs[:,intF]
        extFacesCells = mesh.faceCellIDs[:,extF]

        TKL = mesh._calcFaceAreas()/mesh._calcFaceToCellDistAndVec()[0].sum(axis=0)
        mes_edge = mesh._calcFaceAreas()
        mesK = mesh.cellVolumes

        # ------------------------------------------ The Chemical Potential ------------------------------
        aDiffP = np.zeros((nVol, nVol))
        aDiffP = sp.csc_matrix(aDiffP)
        aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
        aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
        aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
        aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

        # -----------------------------------Neumann Boundary condition------------------------------------------
        aDiffP = aDiffP + sp.coo_matrix((0.*TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

        e = np.ones((1, nVol))
        EaDiffP = sp.csc_matrix(np.concatenate((np.concatenate((d*aDiffP.T.todense(), (mesK*e).T), axis=1),
                                                np.array([np.append((mesK*e).T, 0.)])), axis=0))

        # -----------------------------------Dirichlet Boundary condition------------------------------------------
        test = CellVariable(mesh=mesh, value = 0.)
        phi = CellVariable(name="$\phi(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
        # sF = sigmaF.SigmaF2D(params.sF, xt, yt, Rs=0.05)
        F = sigmaF.SigmaF2D(1./x_s, xt, yt, Rs=0.05/(x_s**2))

        extendedF = np.append(mesK * U * F, 0.)
        phi_new = la.spsolve(EaDiffP,  extendedF)
        phi.setValue(phi_new[0:nVol])
        phi.updateOld()

        # ------------------------------------------ The Chemoattractant ------------------------------
        aDiffE = np.zeros((nVol, nVol))
        aDiffE = sp.csc_matrix(aDiffE)
        aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
        aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
        aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
        aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

        # -----------------------------------Dirichlet Boundary condition------------------------------------------
        aDiffE = aDiffE + sp.coo_matrix((-TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

        aConVE = np.zeros((nVol, nVol))
        aConVE = sp.csc_matrix(aConVE)

        dPhi_int = numerix.dot(phi.faceGrad.value, mesh.faceNormals)[intF]

        aConVE = aConVE + sp.coo_matrix((mes_edge[intF] * plus(dPhi_int), (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
        aConVE = aConVE + sp.coo_matrix((-mes_edge[intF] * minus(dPhi_int), (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
        aConVE = aConVE + sp.coo_matrix((-mes_edge[intF] * plus(dPhi_int), (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
        aConVE = aConVE + sp.coo_matrix((mes_edge[intF] * minus(dPhi_int), (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

        dPhi_ext = numerix.dot(phi.faceGrad.value, mesh.faceNormals)[extF]

        aConVE = aConVE + sp.coo_matrix((mes_edge[extF] * plus(dPhi_ext), (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

        Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

        # ---------------Variables and parameters for the Immune Cells Displacement equation---------
        E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0., hasOld=1)

        delta_tild = CellVariable(name="$\delta_t(x,y)$", mesh=mesh, value=0.)
        delta_tild.setValue(delt.GaussianImmuneResponse2D(1./x_s, xt, yt, Ra=0.02/x_s**2))

        pf = params.pf * np.ones(nVol)
        # gF = params.gf * np.ones(nVol)

        Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)


        # ---------------------------------------------- Dichotomie Method -----------------------------------------------

        mu_a = 0.
        mu_b = 1.

        # F_mu_m = 0.
        # F_mu_a = 0.

        eps = 1e-10
        Tol = 10000

        # Strength = 0.

        mu = dichotomy(mu_a, mu_b, eps, mesK, Q, aDiffE, aConVE, Id, gF, t_s, E, delta_tild)
        while not_converge(mu, mu_b):
            print(not_converge(mu, mu_b))
            mu_b = (mu_a + mu_b)/2.
            mu = dichotomy(mu_a, mu_b, eps, mesK, Q, aDiffE, aConVE, Id, gF, t_s, E, delta_tild)
            print('***({0},{1})'.format(mu, mu_b))
        print('step {0} -> mu = {1}'.format(n, mu))

        # print(mu*mu1_s*1e-9)
        mus.append(2. * np.power(3. * (mu * mu1_s * 1e-9) / (4. * np.pi), 1. / 3.))
        n=n+1
    mus_mat.append(mus)
# np.savetxt("data/mu1_vs_a.txt", mus)
# np.savetxt("data/martix_mu1_vs_a_R.txt", mus_mat)
# np.savetxt("data/martix_mu1_vs_a_A_new.txt", mus_mat)
# np.savetxt("data/martix_mu1_vs_A_R_IC1.txt", mus_mat)
np.savetxt("data/martix_mu1_vs_a_gamma_IC4.txt", mus_mat)
# np.savetxt("data/martix_mu1_vs_A_gamma_IC4.txt", mus_mat)

# Y = np.loadtxt("data/mu1_vs_a.txt")
# Y = np.loadtxt("data/martix_mu1_vs_a_R.txt")
# Y = np.loadtxt("data/martix_mu1_vs_a_A_new.txt")
Y = np.loadtxt("data/martix_mu1_vs_a_gamma_IC4.txt")


# -------------------------------------- Result visualization -----------------------------------

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams.update({'font.size': 22})

fig = plt.figure(2)
# color1 = 'tab:red'
color1 = 'k'
size_of_fig=8
fig.set_size_inches(size_of_fig+6, size_of_fig+2)
for i in range(5):
    peaks = sign.find_peaks(Y[:,i])
    for j in peaks[0]:
        Y[j, i] = (Y[j-1, i] + Y[j+1, i])/2.
# for i in np.arange(0,100, 1):

# tests = [0, 1, 5, 10, 20]
tests = [0, 24, 48, 72, 96]
for i in [0,1,2,3,4]:
    plt.plot(a_s,Y[:,i], color=color1, linewidth=4)
    # if i ==0:
    # plt.text(np.max(a_s) - 3 * (a_s[-1] - a_s[-2]), np.max(Y[:, i])+0.0008,
    #              '$\gamma={0}$'.format(round(A_s[tests[i]], 8)))

    plt.text(np.max(a_s) - 5. * (a_s[-1] - a_s[-2]), np.max(Y[:, i]) ,
                 '$\gamma={0}$'.format(round(gF_s[tests[i]], 8)))
    # else:
    #     plt.text(np.max(delta_s) - 3 * (delta_s[-1] - delta_s[-2]), np.min(Y[:,i])+0.0001, '$R={0}$'.format(round(S_s[tests[i]], 8)))
    # plt.text(np.max(a_s) - 5 * (a_s[-1] - a_s[-2]), np.max(Y[:,i]), '$A={0}$'.format(round(delta_s[tests[i]], 8)))
    plt.xlabel('division rate $a$ in $day^{-1}$')
    # plt.xlabel('Im. cells influx rate $R$ in $\\frac{cell_c \cdot mm^{-3}}{\mu m^3} \cdot day^{-1}$')
    # plt.xlabel('Im. cells death rate. $\gamma$ in $day^{-1}$')
    # plt.xlabel('Im. strength $A$ in $cell_c^{-1} \cdot day^{-1}$')
    plt.ylabel('tumor diameter in $mm$')




