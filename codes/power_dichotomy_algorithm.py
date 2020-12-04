#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fipy import *
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
from matplotlib import pyplot as plt
import parameterFunctions.immuneResponse as delt
import parameterFunctions.sigmaF as sigmaF
import parameterFunctions.divisionOperatorKernel as K
import parameterFunctions.TestFeature as TF

a = 4.
V = .616
delta = 1.
pf = 0.25
S = 20.
chi = 8.64e-1
# chi = 1000
D = 0.025
sF = 200.e-5
gF = 0.18
# d = 2.16
d = 1.
init_value = 5.
alpha = 0.
beta = 0.

# tested_param = ''
tested_param_name = 'a'
params = TF.TestFeature(a, V, delta, pf, S, chi, D, sF, gF, d, tested_param_name, init_value, alpha, beta)


def norm_n(V, dx, n):
    if n == 0:
        # maximum des valeurs absolues
        c_max = np.max(np.abs(V))
        yield c_max
    else:
        norme = np.sum(np.abs(dx*V)**n)
        # racine n-ieme
        yield norme**(1./n)


def mo(x):
    return numerix.L2norm(x)


def plus(z):
    return 0.5*(z+np.abs(z))


def minus(z):
    return 0.5*(-z+np.abs(z))

# Stationary profile


def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)

if params.a/params.V > 1.:
    Lz = 5.
    nz = 7000
elif params.a/params.V < 1.:
    Lz = 100.
    nz = 7000
dz = Lz / nz
meshz = Grid1D(dz, nz)
z, = meshz.cellCenters


cellSize = .1
radius = 1.

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
# aDiffP = aDiffP + sp.coo_matrix(-TKL[extF], (extFacesCells[0], extFacesCells[0]), shape=(nVol, nVol))

test = CellVariable(mesh=mesh, value = 0.)
phi = CellVariable(name="$\phi(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
sF = sigmaF.SigmaF2D(params.sF, xt, yt, Rs=0.05)

F = sF
extendedF = np.append(mesK *F, 0.)
phi_new = la.spsolve(d * EaDiffP,  extendedF)
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
Idz = sp.spdiags(numerix.ones(nz), [0], nz, nz)

# ---------------Variables and parameters for the Immune Cells Displacement equation---------
E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
E_p = CellVariable(name="$\partial_{\mu_1} E(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
# S = CellVariable(name="$S(x,y)$",mesh=mesh, value=0.0, hasOld=1)
S = CellVariable(mesh=mesh, value=0.0, hasOld=1)

# center_x1 = 0.207
# center_x2 = 0.907
# S.setValue(params.S * numerix.exp(-(((xt - 0.) ** 2)/0.5 +((yt - 1.) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt - center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt + center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125))
# S.setValue(S/numerix.sum(mesK*S))


S.setValue(params.S / np.sum(mesK))

# nS = numerix.sum(mesK*S).value
# ring_vol = numerix.sum(mesK[numerix.where(((xt**2 + yt**2) > 0.5) & ((xt**2 + yt**2) < 0.75))])
# nS_per_vol = nS/ring_vol
#
# S.setValue(nS_per_vol, where=((xt**2 + yt**2) > 0.5) & ((xt**2 + yt**2) < 0.75))
# S.setValue(0., where=(xt**2 + yt**2) <= 0.5)
# S.setValue(0., where=(xt**2 + yt**2) >= 0.75)

delta = CellVariable(name="$\delta(x,y)$", mesh=mesh, value=0.)
delta.setValue(delt.GaussianImmuneResponse2D(params.delta, xt, yt, Ra=0.02))

# pf = pF.GaussianConversionRate2D(amppf, xt, yt, mean=0, Rp=1.)
pf = params.pf * np.ones(nVol)
# gF = gammaF.gaussianGammaF2D(ampgF, xt, yt, Rs=0.3)
gF = params.gf * np.ones(nVol)

# ---------------------------------- For the tumor stationnary problem -------------------------------------------------
a = CellVariable(mesh=meshz, value=params.a)
# a.setValue(params.a*z)
V = CellVariable(mesh=meshz, value=params.V)
V1 = V.harmonicFaceValue.value
vPlus = plus(V1)
vMinus = minus(V1)
# ConvT = sp.lil_matrix(sp.spdiags([-vPlus[0:nz], (vPlus[1:nz+1] + vMinus[0:nz]), -vMinus[1:nz+1]], [-1, 0, 1], nz, nz))

GPT = vPlus
GMT = vMinus
GP2T = np.concatenate([[0], vPlus[1:nz]])
GM2T = np.concatenate([[0], vMinus[1:nz]])
ConvT = sp.lil_matrix(sp.spdiags([-GPT[1:nz + 1], GPT[1:nz + 1] + GMT[0:nz], -GM2T], [-1, 0, 1], nz, nz))

G = sp.csc_matrix(K.loopGaussianKernel(z, dz))
# G = sp.csc_matrix(K.UniformKernel(z))
div = sp.lil_matrix(sp.spdiags(a.value, [0], nz, nz))


Coef = 1. + 4.*dz*np.max(np.sum(G.multiply(a.value), axis=1))*np.max(V.value)/np.min(V.value) - np.min(a.value)
# Coef = 1. + dz*np.max(np.sum(G.multiply(a.value), axis=1))*np.max(V.value)/np.min(V.value) - np.min(a.value)


# gain = sp.lil_matrix(sp.spdiags(4.*np.sum(G.multiply(dz*a.value), axis=1), [0], nz, nz))
operatorMat = ConvT.multiply(1./dz).tocsc() + div.tocsc() - 4*G.multiply(dz*a.value)
# operatorMat = ConvT.multiply(1./dz).tocsc() + div.tocsc() - 2*G.multiply(dz*a.value)
# operatorMat = ConvT.multiply(1./dz).tocsc() + div.tocsc() - 2.*numerix.dot(G.multiply(dz), div)
A = operatorMat.tocsc() + Idz.multiply(Coef)

Tol = 10000
# Vect = np.random.rand(nz, 1)
Vect = np.ones((nz, 1))
Vold = Vect
compt = 0

# --------------------------------------------------------- The power method ------------------------------------------

AA = la.splu(A)
eps = 1e-8
while (Tol > eps) & (compt < 5000):
    compt = compt + 1
    Vect = np.copy(AA.solve(Vect))
    Vect = np.copy(Vect / (dz * np.sum(Vect)))
    Tol = np.copy(next(norm_n(Vect - Vold, dz, 1)) / next(norm_n(Vect, dz, 1)))
    Vold = np.copy(Vect)
Vect = np.copy(AA.solve(Vect))
Lambda = np.copy(Vect.T.dot(Vold) / np.float64(Vold.T.dot(Vold)))
Lambda = np.copy(1. / Lambda)
Lambda = np.copy(Lambda - Coef)

# test_case_results = np.empty([], dtype=[('strength', np.float64),
#                                         # ('dz', np.float64),
#                                         # ('N', np.float64, (nz,)),
#                                         ('C', np.float64, (nVol,)),
#                                         ('phi', np.float64, (nVol,))])
# test_case_results = np.delete(test_case_results, 0)

eps = 1e-9
Tol = 10000

Strength = 0.
n = 0

# -------------------- Functions for testing the numerical methods ---------------------
# C_test = CellVariable(name='$C$', mesh=mesh, value=0.)
# C_test.setValue((1-mo(x))*numerix.exp(-mo(x)))
#
# phi_test = CellVariable(name='$\phi$', mesh=mesh, value=0.)
# phi_test.setValue(1 - 2*mo(x))
# phi_test.setValue(np.sin(xt)*np.sin(yt) - 4*np.sin(0.5)**4)
# phi_test.setValue((1./3.)*mo(x)**3 - 0.5*mo(x)**2)
# phi_test.setValue((1./36.)*mo(x)**4 + (1./3.)*mo(x)**3 - (10./9.)*mo(x))

# Source = CellVariable(name='$C$', mesh=mesh, value=0.)

# ---------------------------------------------- Newton's Method -----------------------------------------------

mu_new = 0.
mu_old = 100.
F_mu = 10.
F_mu_p = 0.
mus = []
# while np.abs(mu_new - mu_old) > eps:
while np.abs(F_mu) > eps:
    mus.append([mu_old])
    bE = (mesK * mu_old * pf * S.value).T
    E_new = la.spsolve(-D * (aDiffE.tocsc()) - params.chi * mu_old * aConVE + Id.multiply(mesK * gF), bE)
    E.setValue(E_new)
    E.updateOld()
    # -------------------------Solve the Tumor Growth equation---------------------------
    Strength = numerix.sum(mesK * delta.value * E.value)
    F_mu = Strength - np.abs(Lambda)

    bE_p = (mesK * mu_old * pf * S.value).T + aConVE.dot(E.value)
    E_p_new = la.spsolve(-D * (aDiffE.tocsc()) - params.chi * mu_old * aConVE + Id.multiply(mesK * gF), bE_p)
    E_p.setValue(E_p_new)
    E_p.updateOld()

    F_mu_p = numerix.sum(mesK * delta.value * E_p.value)
    # F_mu = Strength - Lambda

    mu_new = np.float64(mu_old - F_mu/F_mu_p)

    print('mu_k+1 - mu_k: ', np.abs(mu_new - mu_old), 'etape: ', n, 'mu_k: ', mu_old)
    mu_old = mu_new
    # test_case_results = np.append(test_case_results, np.asarray((Strength,
    #                                                              # dt,
    #                                                              # dz,
    #                                                              # mu0_new.value * Vect,
    #                                                              E.value,
    #                                                              phi.value), dtype=test_case_results.dtype))


plt.title("$\int_{\Omega} \delta(y) C_{\mu}(y)dy$")
plt.xlabel("$\mu$")
plt.grid()
plt.plot(mus)
plt.legend()

# last mu = 1.7554120763635088
# 1.15576025235



# ---------------------------------------------- Dichotomy Method -----------------------------------------------

mu_a = 0.
mu_b = 10000.

F_mu_m = 0.
F_mu_a = 0.

mus = []
# E_a = CellVariable(name="$E_a(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
# E_m = CellVariable(name="$\partial_{\mu_1} E_b(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
# Lambda = 0.0625
while mu_b - mu_a > eps:
    mu = (mu_b + mu_a)/2.
    mus.append([mu])
    bE = (mesK * mu * pf * S.value).T
    E_new = la.spsolve(-D * (aDiffE.tocsc()) - params.chi * mu * aConVE + Id.multiply(mesK * gF), bE)
    E.setValue(E_new)
    E.updateOld()

    F_mu_m = numerix.sum(mesK * delta.value * E.value) - np.abs(Lambda)

    bE = (mesK * mu_a * pf * S.value).T
    E_new = la.spsolve(-D * (aDiffE.tocsc()) - params.chi * mu_a * aConVE + Id.multiply(mesK * gF), bE)
    E.setValue(E_new)
    E.updateOld()

    F_mu_a = numerix.sum(mesK * delta.value * E.value) - np.abs(Lambda)

    if F_mu_m * F_mu_a <= 0:
        mu_b = mu
    else:
        mu_a = mu

    print (mu)


# 5.3063631592230722

