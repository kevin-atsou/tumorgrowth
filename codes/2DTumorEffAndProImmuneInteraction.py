#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fipy import *
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
import scipy.linalg as lin
from matplotlib import pyplot as plt
import sys
import os
import time
import parameterFunctions.divisionRate as ad
import parameterFunctions.immuneResponse as delt
import parameterFunctions.gammaF as gammaF
import parameterFunctions.pF as pF
import parameterFunctions.sigmaF as sigmaF
import parameterFunctions.kernels as kern
import parameterFunctions.divisionOperatorKernel as K
import parameterFunctions.TestFeature_ProTumor as TFP
# from fipy.solvers import _MeshMatrix
# from fipy.variables.cellVariable import CellVariable

# sys.setrecursionlimit(10000)


cellSize = .05
radius = 1.


def plus(z):
    return 0.5*(z+np.abs(z))


def psi(x, y):
    # return (Abs(x - y) + x - y) * Exp(-1. / (x - y))/ 2.
    # if x -y > 0.:
    #     return x-y
    # else:
    #     return 0.
    return (np.abs(x - y) + x - y) / 2.


def minus(z):
    return 0.5*(-z+np.abs(z))

# Stationary profile


def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)

#
# Lz = 6.0
# nz = 3334
# dz = Lz / nz

# nmax = 100000
# meshtemp = Grid1D(dz, nmax)
# ztmp, = meshtemp.cellCenters


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

# mesh = Gmsh2D('disk2.msh', communicator=serialComm)

# mesh = Gmsh2D('disk2.msh', communicator=serialComm)
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

# ----------------------Configuration the test cases and the results saving --------------------
# (self, a, V, delta, pf, S, chi, D, sF, gF, k)
# Parameters for the tumor growth
a = 0.8
# a = 1.
V = .616

# Parameters for the effector cells
delta = 1.
# delta = 0.01
pf = 0.25
# pf = 10.
S = 20.
chi = np.float64(8.64e-1)
# chi = 1000
D = 0.025
sF = np.float64(200.e-5)
gF = 0.18
# d = 2.16
d = 1.

# Parameters for the pro-tumor immune cells
# b1 = 1.
b1 = 1e-6
# b2 = 0.3
b2 = 0.0
eta = 10.
# sF2 = 1.12e-2
sF2 = 1.

# tau = 1.155e-2
tau = 1.
gammar = 0.18
kr = 1.
theta1 = 0.2
theta2 = 0.2
pfr = 0.20
kc = 1.
# mc = 9.1
mc = 2.0
init_value = 5.
alpha = 0.
beta = 1.
Sr = S
# tested_param = ''
tested_param_name = 'delta'
params = TFP.TestFeature(a, V, delta, pf, S, chi, D, sF, gF, d, b1, b2, eta, sF2, tau, gammar, kr, theta1, theta2, pfr, kc
                         , mc, tested_param_name, init_value, alpha, beta, Sr)


# Tbar = 0.
# for n in np.arange(20):
#     Tbar = Tbar + ((-1) ** n) * alphan(n) * numerix.exp(
#         -(2. ** (n + 1)) * params.a * ztmp / params.V)
# nztmp = np.size(Tbar[Tbar > 1e-04])

# Re-adjust The tumor size computational domain with the asymptotic profile

# pro-tumor immune cells parameters. We put it here because we it for the Re-adjustment
# of the tumor size computational domain

b1 = CellVariable(name="$b(x,y)$", mesh=mesh, value=0.)
b1.setValue(kern.GaussianKernel2D(params.b1, xt, yt, Ra=0.03))

b2 = CellVariable(name="$b(x,y)$", mesh=mesh, value=0.)
b2.setValue(kern.GaussianKernel2D(params.b2, xt, yt, Ra=0.03))

b1_max = np.max(b1.value)
b2_max = np.max(b2.value)


theta1 = CellVariable(name="$\\theta1(x,y)$", mesh=mesh, value=0.)
theta1.setValue(kern.GaussianKernel2D(params.theta1, xt, yt, Ra=0.3))

theta2 = CellVariable(name="$\\theta2(x,y)$", mesh=mesh, value=0.)
theta2.setValue(kern.GaussianKernel2D(params.theta2, xt, yt, Ra=0.3))

kr = params.kr * np.ones(nVol)
kc = params.kc * np.ones(nVol)
pfr = params.pfr * np.ones(nVol)
gammar = params.gammar * np.ones(nVol)

# -------------------------Computation of the solutions in time-------------------------

# Lz = 5.0
# nz1 = 5000
# dz = Lz / nz1

# dztmp = 0.003

# -------------------------Computation of the solutions in time-------------------------
CFL = .8

# THe time step
# dt = min(3e-3, CFL * dz / params.V);
# tt = 0.
# tps = np.arange(0., 1, 1)
# Tf = 200.
# n = 0

mesh_temp = 0.
n_max = 0
if params.a > params.V*(1.):
    # if params.a > params.V:
    n_max = 10000
    # dz = dz / (2. * params.a)
    Lz = 10.
    dz = Lz/n_max
    mesh_temp = Grid1D(dz, n_max)
# if params.a < params.V:
if params.a < params.V*(1.):
    n_max = 7000
    Lz = 10.
    dz = Lz / n_max
    mesh_temp = Grid1D(dz, n_max)


# z_tmp, = mesh_temp.cellCenters

# We add some points to be a little broad
# bonus = 1000

# t_bar = 0.
# for n in np.arange(20):
#     t_bar = t_bar + ((-1) ** n) * alphan(n) * numerix.exp(
#         -(2. ** (n + 1)) * params.a * z_tmp / (params.V*(1. + b1_max)))
# t_bar = t_bar / (dz * np.sum(t_bar))
# nz_tmp = np.size(t_bar[t_bar > 1e-06]) + bonus
# print ('taille: ' + str(nz_tmp))

# Re-adjust The tumor size computational domain
# nz = nz_tmp
nz = n_max
print ('taille: ' + str(nz))
# meshz = Grid1D(dz, nz)
# meshz = mesh_temp
z, = mesh_temp.cellCenters
z_max = np.max(z)

# Compute New t_bar with the new domain
t_bar = 0.
for n in np.arange(20):
    t_bar = t_bar + ((-1) ** n) * alphan(n) * numerix.exp(
        -(2. ** (n + 1)) * params.a * z / (params.V*(1. + b1_max)))
t_bar = t_bar / (dz * np.sum(t_bar))

# Compute the division kernel
try:
    if nz < 5000:
        G = sp.csc_matrix(K.gaussianKernel(z, dz))
    else:
        G = sp.csc_matrix(K.loopGaussianKernel(z, dz))
except MemoryError:
    print ('Memory Error, Using the loop...')
    G = sp.csc_matrix(K.loopGaussianKernel(z, dz))

# ----------------Variables and parameters for the Tumor growth equation--------------------
# const = 0.1
# Tau = 1.
T = CellVariable(name="$T(t,z)$", mesh=mesh_temp, value=0.0, hasOld=1)
# T.setValue(5.*(4./(z_max**2)), where=(z > 0.0625) & (z < z_max/2.))
T.setValue(5./((z_max/3.)**2 - 0.125**2), where=(z > 0.125) & (z < z_max/3.))

V = CellVariable(name="$V$", mesh=mesh_temp, value=params.V)
# V.setValue(params.V*z*np.log(z_max/z))
V1 = V.harmonicFaceValue.value

a = CellVariable(name="$a$", mesh=mesh_temp, value=params.a)
a.setValue(0., where=z <= z_max / 10.)
a.setValue(params.a, where=z > z_max / 10.)
# a.setValue(params.a * np.exp(-((z-np.max(z)/2)**2)/(2*1e-3))/np.sqrt(2*np.pi))
# Convection Matrix construction
vPlus = plus(V1)
vMinus = minus(V1)
GP = vPlus
GM = vMinus
GP2 = np.concatenate([[0], vPlus[1:nz]])
GM2 = np.concatenate([[0], vMinus[1:nz]])
conVT = sp.lil_matrix(sp.spdiags([-GP[1:nz + 1], GP[1:nz + 1] + GM[0:nz], -GM2], [-1, 0, 1], nz, nz))


# The 2z matrix
# G = sp.csc_matrix(K.gaussianKernel(z))
# G[G > 0.] = 1.

# ---------------Variables and parameters for the ChemoAttractant diffusion---------------
test = CellVariable(mesh=mesh, value = 0.)
phi = CellVariable(name="$\phi(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
sF = sigmaF.SigmaF2D(params.sF, xt, yt, Rs=0.05)
# Diffusion Matrix construction
# x = mesh.cellCenters
# xt, yt = mesh.cellCenters

# ------------------------------------------ THE CHEMICAL POTENTIAL ------------------------------
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
F = sF
extendedF = np.append(mesK *F, 0.)
# extendedF = np.append(F, 0.)
phi_new = la.spsolve(d * EaDiffP,  extendedF)

phi.setValue(phi_new[0:nVol])
phi.updateOld()

# -----------------------------------Dirichlet Boundary condition------------------------------------------
# aDiffP = aDiffP + sp.coo_matrix(-TKL[extF], (extFacesCells[0], extFacesCells[0]), shape=(nVol, nVol))


# ------------------------------------------ THE IMMUNE CELLS ------------------------------

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

# dPhi_ext = (phi.value[extFacesCells[1]] - phi.value[extFacesCells[0]])
# dPhi_ext = numerix.dot(phi.grad.harmonicFaceValue.value,(mesh.faceNormals))._get_data()[extF]
dPhi_ext = numerix.dot(phi.faceGrad.value, mesh.faceNormals)[extF]

aConVE = aConVE + sp.coo_matrix((mes_edge[extF] * plus(dPhi_ext), (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)
Idz = sp.spdiags(numerix.ones(nz), [0], nz, nz)

# ---------------Variables and parameters for the Immune Cells Displacement equation---------
E = CellVariable(name="$c(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
B = CellVariable(name="$c_r(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
# E.cacheMe(False)
# B.cacheMe(False)
# T.cacheMe(False)
# S = CellVariable(name="$S(x,y)$",mesh=mesh, value=0.0, hasOld=1)
S = CellVariable(mesh=mesh, value=0.0, hasOld=1)
# S = CellVariable(name="$S(t,x,y)$", mesh=mesh, value=1.35e-4, hasOld=1)
# Test values for S
# S.setValue(1.0, where=xt<-3)
# S.setValue(1.35)
# center_x1 = 0.207
# center_x2 = 0.907
# S.setValue(params.S * numerix.exp(-(((xt - 0.) ** 2)/0.5 +((yt - 1.) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt - center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt + center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125))
# S.setValue(S/numerix.sum(mesK*S))
S.setValue(params.S / np.sum(mesK))

Sr = CellVariable(mesh=mesh, value=0.0, name="$S_r(x,y)$", hasOld=1)
# S = CellVariable(name="$S(t,x,y)$", mesh=mesh, value=1.35e-4, hasOld=1)
# Test values for S
# S.setValue(1.0, where=xt<-3)
# S.setValue(1.35)
center_x1 = 0.207
center_x2 = 0.907
Sr.setValue(params.Sr*(numerix.exp(-(((xt - 0.) ** 2)/0.5 +((yt - 1.) ** 2))/0.125) +
            numerix.exp(-(((xt - center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125) +
            numerix.exp(-(((xt + center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125)))
Sr.setValue(Sr/numerix.sum(mesK*Sr))
# S.setValue(params.S / np.sum(mesK))


# nS = numerix.sum(mesK*S).value
# ring_vol = numerix.sum(mesK[numerix.where(((xt**2 + yt**2) > 0.5) & ((xt**2 + yt**2) < 0.75))])
# nS_per_vol = nS/ring_vol
#
# S.setValue(nS_per_vol, where=((xt**2 + yt**2) > 0.5) & ((xt**2 + yt**2) < 0.75))
# S.setValue(0., where=(xt**2 + yt**2) <= 0.5)
# S.setValue(0., where=(xt**2 + yt**2) >= 0.75)
Sr.updateOld()

delta = CellVariable(name="$\delta(x,y)$", mesh=mesh, value=0.)
# delta.setValue(params.delta)
delta.setValue(delt.GaussianImmuneResponse2D(params.delta, xt, yt, Ra=0.02))

# pf = pF.GaussianConversionRate2D(amppf, xt, yt, mean=0, Rp=1.)
pf = params.pf * np.ones(nVol)
# gF = gammaF.gaussianGammaF2D(ampgF, xt, yt, Rs=0.3)
gF = params.gf * np.ones(nVol)

# -------------------------The Chemokines Signal I -----------------------------------
I = 0.
I_new = 0.
Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# -------------------------Solutions representation -----------------------------------
# fig, axes = plt.subplots(1, 2)
# cmap = plt.get_cmap("Spectral_r")
# Tfigure = Viewer(vars=T, datamax=np.max(T.value))
# Tfigure = Matplotlib1DViewer(vars=T, axes=axes[0], datamin=0., figaspect='auto', cmap=cmap)
# Efigure = Matplotlib2DViewer(vars=E, axes=axes[1], figaspect='auto', cmap=cmap)
# Bfigure = Matplotlib2DViewer(vars=B, axes=axes[2], figaspect='auto', cmap=cmap)
# axes[0].set_title('The Tumor $T$',fontsize=10, fontweight='bold')
# axes[1].set_title('The Im. Cells',fontsize=10, fontweight='bold')
# axes[2].set_title('The Pro tum.',fontsize=10, fontweight='bold')
# viewers = MultiViewer(viewers=(Tfigure, Efigure))
# viewers = MultiViewer(viewers=(Tfigure))


# ----------------------Configuration the test case results saving --------------------
test_case_results = np.empty([], dtype=[(params.tested_param, np.float64),
                                        ('t', np.float64),
                                        ('dt', np.float64),
                                        ('dz', np.float64),
                                        ('T', np.float64, (nz,)),
                                        ('E', np.float64, (nVol,)),
                                        ('B', np.float64, (nVol,)),
                                        ('I', np.float64),
                                        ('phi', np.float64, (nVol,))])
test_case_results = np.delete(test_case_results, 0)
# TestCaseParam = np.empty([], dtype=[('V', np.float64),
#                                     ('dat', TestCaseResults.dtype, ())])
# if sys.argv[1] == 0:
# typeOfS = 'homo'
# else:
typeOfS = 'homo'
#
# filename = '2D'+typeOfS+'SCParams' + params.tested_param +'_1'
# dataPath = "data/" + filename + ".npz"

newpath = 'TestResults_New/NSINSC/' + typeOfS + 'geneousS/parameter' + params.tested_param
if not os.path.exists(newpath):
    os.makedirs(newpath)

# -------------------------Computation of the solutions in time-------------------------
CFL = .8

# THe time step
dt = min(5e-4, CFL * dz / V.value[0]);
tt = 0.
tps = np.arange(0., 1, 1)
Tf = 90.
n = 0

Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# compute the number of cells in the tumor
mu0_new = dz * numerix.sum(T)
# compute the volume of the tumor
mu1_new = dz * numerix.sum(z * T)
# compute the number of Immune cells
mu0E_new = numerix.sum(mesK * E)
mu1E_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * E)

# compute the number of pro-tumor Immune cells
mu0B_new = numerix.sum(mesK * B)
mu1B_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * B)
# Saving the moments at time t0
# mu1s.append(mu1_new.value)
# mu1Es.append(mu1E_new.value)

Strength = 0.0
#
# dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/with_growth_promotion/V_constant/with_Sr/parameter' + params.tested_param
dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/with_growth_promotion/V_gompertz/with_Sr/parameter' + params.tested_param
filename = '2D' + typeOfS + 'SCParams' + params.tested_param + str(
    params.getAttributeByName(params.tested_param)) \
           # + 'b1'
           # + 'g_mu1'
    # params.getAttributeByName(params.tested_param)) + 'b1'+str(params.b1)
dataPath = dirpath + '/' + filename + ".txt"

if not os.path.exists(dirpath):
    os.makedirs(dirpath)

outF = open(dataPath, "w")

eps = 1e-10
while tt < Tf:
    if np.isnan(mu1B_new.value):
        break
    if n % 10 == 0:
        test_case_results = np.append(test_case_results, np.asarray((params.getAttributeByName(params.tested_param),
                                                                     tt,
                                                                     dt,
                                                                     dz,
                                                                     T.value,
                                                                     E.value,
                                                                     B.value,
                                                                     I_new,
                                                                     phi.value), dtype=test_case_results.dtype))
        outF.write(str(tt)+' '+str(dt)+' '+str(dz)+' '+str(mu0_new.value)+' '+str(mu1_new.value)+' '+str(mu1E_new.value)
                   +' '+str(mu1B_new.value)+' '+str(Strength)+' '+str(I))
        # outF.writelines()
        outF.write('\n')

    # Bconversion = numerix.sum(mesK * theta1 * B)
    Bconversion = theta1 * E
    # viewers.plot()
    # -----------------Solve the Immune Cells displacement equation---------------------
    # bE = ((mesK / dt) * E.value + mesK * mu1_new.value * pf * S - mesK * params.kr * I * Bconversion.value -
    # kc*theta2*B.value*E.value)

    bE = ((mesK / dt) * E.value + mesK * (mu1_new.value) * pf * S)
    # bE = ((mesK/dt)*E.value - mesK*gF*E.value + mesK * mu1_old.value * pf * S.value).T
    # time_start = time.clock()
    E_new = la.spsolve((Id.multiply(mesK / dt) - D * aDiffE.tocsc() - params.chi * (mu1_new.value/(params.beta + mu1_new.value)) * aConVE
                        + Id.multiply(mesK * gF) + Id.multiply(mesK * kr * I * theta1) + Id.multiply(kc*mesK*B.value)),
                       bE)
    # E_new = la.inv(sp.coo_matrix(Id.tocsc() + dt*(aDiffE.tocsc()/mesK))).dot(bE)
    # print("--- %s seconds ---" % (time.time() - time_start))

    # if any(E_new < 0):
    #     break

    E.setValue(E_new)
    E.updateOld()

    # -----------------Solve the Pro-tumor Immune Cells displacement equation---------------------

    bB = ((mesK / dt) * B.value + mesK * I * pfr * Sr + mesK * params.kr * I * Bconversion.value)
    # bB = ((mesK / dt) * B + mesK * I_new * pfr * S + mesK * kr * I_new * Bconversion)
    # bE = ((mesK/dt)*E.value - mesK*gF*E.value + mesK * mu1_old.value * pf * S.value).T
    # time_start = time.clock()
    B_new = la.spsolve((Id.multiply(mesK / dt) - D * aDiffE.tocsc() - params.chi * (mu1_new.value/(params.beta + mu1_new.value)) * aConVE
                        + Id.multiply(mesK * gammar)),
                       bB)
    # E_new = la.inv(sp.coo_matrix(Id.tocsc() + dt*(aDiffE.tocsc()/mesK))).dot(bE)
    # print("--- %s seconds ---" % (time.time() - time_start))

    # if any(E_new < 0):
    #     break

    B.setValue(B_new)
    B.updateOld()

    # np.max(B.value)
    # -------------------------Solve the Tumor Growth equation---------------------------
    Strength = numerix.sum(mesK * delta * E).value
    B1effect = numerix.sum(mesK * b1 * B)
    # B1effect = numerix.sum(mesK * b1 * B/(params.eta + B))
    B2effect = numerix.sum(mesK * b2 * B/(params.eta + B))

    # AS = sp.lil_matrix((sp.spdiags(a.value*(1. + B2effect.value) + Strength.value, [0], nz, nz)))
    # AS = sp.lil_matrix((sp.spdiags(a.value + Strength.value, [0], nz, nz)))

    # Explicit scheme for the tumor growth
    AS = sp.lil_matrix((sp.spdiags(a.value + Strength, [0], nz, nz)))
    T_new = (Idz - (dt / dz) * (conVT.multiply(1. + B1effect.value).tocsc() + dz * AS.tocsc())).dot(T.value) + 4 * dt * dz * G.dot(
        a.value * T.value)
    # T_new = la.spsolve(Idz + (dt / dz) *conVT.tocsc(), (Idz - (dt / dz) * (dz * AS.tocsc())).dot(T.value) + 4 * dt * dz * G.dot(
    #     a.value * T.value))
    T.setValue(T_new)
    T.updateOld()

    # Implicit scheme for the tumor growth
    # BT = T.value + 4. * dt * dz * G.dot(a.value*(1. + B2effect.value) * T.value)
    # BT = T.value + 4. * dt * dz * G.dot(a.value * T.value)
    # T_new = sp.linalg.spsolve((Idz + (dt / dz) * conVT.multiply(1. + B1effect.value).tocsc() + (dt / dz) * (dz * AS.tocsc())), BT)
    # # T_new = sp.linalg.spsolve((Idz + (dt / dz) * conVT.tocsc() + (dt / dz) * (dz * AS.tocsc())), BT)
    #
    # T.setValue(T_new)
    # T.updateOld()

    # compute the number of cells in the tumor
    # mu0_new = dz * numerix.sum(T)
    # compute the volume of the tumor
    # mu1_new = dz * numerix.sum(z * T)
    # compute the number of Immune cells
    # mu0E_new = numerix.sum(mesK * E)
    mu1E_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * E)
    mu1B_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * B)

    # compute the volume of the tumor
    mu1_new = dz * numerix.sum(z * T)
    mu0_new = dz * numerix.sum(T)

    # -------------------------Solve the Chemokines signal equation---------------------------
    I_new = (1. - dt*params.tau)*I + dt*psi(mu1_new.value, mc)*params.sF2
    # I_new = (I + dt * psi(mu1_new.value, params.mc) * params.sF2) / (1. + dt * params.tau)
    # I_new = (I + dt * psi(mu1_new.value, params.mc) * params.sF2) / (1. + dt * params.tau)
    I = I_new
    # if np.abs(Strength-max(a)) < eps:
    #     break

    # Setting the CFL Condition
    # Flux = params.chi * (mu1_new.value/(params.beta + mu1_new.value)) * aConVE

    # dtE = 5e-4
    dtE = 5e-2
    # dtE = CFL * np.min(mesK)/np.max(Flux)

    # dtT = CFL * dz / np.max(np.abs(V.value*(1. + b1_max)))
    dtT = CFL * dz / np.max(np.abs(V.value*(1. + B1effect.value)))
    # dtT = 5e-2

    dtI = CFL*1./params.tau
    # dtT = 0.1
    # dtT = 3e-3

    dt = min(dtT, dtE, dtI)
    if n % 100 == 0:
        print('time: ', tt, 'etape: ', n)
        # print ("dtE: {0} dtT: {1} ".format(dtE, dtT))
        print ('dt = {0}'.format(dt))
        print ('mu1 = {0} mu1E = {1} mu1B = {2} strength = {3} I = {4}'.format(mu1_new, mu1E_new, mu1B_new,  Strength, I))

    tt = tt + dt
    # viewers.plot()
    # Tfigure.setLimits(datamax=np.max(T.value))
    # Tfigure.plot()
    n = n + 1
    # if mu1_new > 1e11:
    #     break
outF.close()

# dirpath = 'TestResults_Protumor/' + typeOfS + 'geneousS/deltaVariable/avec_masse_critique/parameter'
# + params.tested_param
# filename = '2D' + typeOfS + 'SCParams' + params.tested_param + str(
#     params.getAttributeByName(params.tested_param))
# dataPath = dirpath + '/' + filename + ".dat"
# if not os.path.exists(dirpath):
#     os.makedirs(dirpath)
# with open(dataPath, "wb") as fi:
#     # global test_case_results
#     print ('... saving ' + params.tested_param + ' in ' + filename)
#     # np.save(fi, TestCaseParam, test_case_results)
#     np.savez_compressed(fi, test_case_results)
#     fi.close()
# data = np.loadtxt(dataPath)
imStrength = np.sum(mesK*delta.value*test_case_results['E'], axis=1)
mu1 = np.sum(z*dz*test_case_results['T'], axis=1)
muB = np.sum((np.sqrt(xt ** 2 + yt ** 2)) * mesK * test_case_results['B'], axis=1)

ts = test_case_results['t']
Is = test_case_results['I']

n=0
n=n+1
fig, ax1 = plt.subplots()
# fig.set_size_inches(size_of_fig+2, size_of_fig)
color1 = 'tab:red'
ax1.set_xlabel('t')
ax1.set_ylabel('$\mu_1$', color=color1)
ax1.plot(ts, mu1, '-', color=color1, linewidth=2, label='$t \mapsto \mu_1(t)$')
# ax1.legend()
ax1.tick_params(axis='y', labelcolor=color1)

ax2 = ax1.twinx()

color2 = 'tab:blue'
color3 = 'tab:olive'
ax2.set_ylabel('$\overline{\mu_{c}}$', color=color2)
ax2.plot(ts, imStrength, '--', marker='o', markersize=2, color=color2, linewidth=4, label='$t \mapsto \overline{\mu_{c}}$')

# ax2.plot(ts, params.a*np.ones(ts.shape[0]), '--', marker='o', markersize=0, color=color2, linewidth=4)
ax2.plot(ts, params.a*np.ones(ts.shape[0]), '--', marker='o', markersize=0, color='k', linewidth=4)
ax2.plot(ts, muB, '--', marker='o', markersize=0, color='orange', linewidth=4)

ax2.tick_params(labelcolor=color2)
fig.tight_layout()

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams.update({'font.size': 18})
#
# plt.plot(imStrength/np.max(imStrength), '--', label='$\mu_c$', linewidth=4)
# plt.plot(mu1/np.max(mu1), label='$\mu_1$', linewidth=4)
# plt.plot(muB/np.max(muB), label='$\mu_{c_r}$', linewidth=4)
# plt.plot(Is/np.max(Is), label='$I$', linewidth=4)
#
# # plt.plot(imStrength, label='$\mu_c$')
# # plt.plot(mu1, label='$\mu_1$')
# # plt.plot(muB, label='$\mu_{c_r}$')
# # plt.plot(Is, label='$I$')
#
# plt.legend()
# dataPath = 'TestResults/New_tests/'
# filename = '2DHeteroLongTimeTest'
# with open(dataPath+filename, "wb") as fi:
#     # global test_case_results
#     print ('... saving ' + params.tested_param + ' in ' + filename)
#     # np.save(fi, TestCaseParam, test_case_results)
#     np.savez_compressed(fi, test_case_results)
#     fi.close()
#
#
# # data = np.load(dataPath+filename)
#
# def norm_n(V, dx, n):
#     if n == 0:
#         # maximum des valeurs absolues
#         c_max = np.max(np.abs(V))
#         return c_max
#     else:
#         norme = np.sum(np.abs(dx*V)**n)
#         # racine n-ieme
#         return norme**(1./n)
#
#
# Coef = 10
#
# div = sp.lil_matrix(sp.spdiags(params.a, [0], nz, nz))
# operatorMat = conVT.multiply(1./dz).tocsc() + div.tocsc() - 4*numerix.dot(G.multiply(dz), div)
# A = operatorMat.tocsc() + Idz.multiply(Coef)
#
#
# eps = 1e-6
# # x= np.random.random((nz,1))
# # x_init = np.ones((nz, 1))
#
# # x_init = np.eye(nz, 1)
# # e, v = power_method(A, x_init, eps)
#
# Tol = 10000
# Vect = np.random.rand(nz, 1)
# # Vect = np.ones((nz, 1))
# Vold = Vect
# compt = 0
#
# AA = lin.lu(A.todense())
# AA_L = sp.lil_matrix(AA[1]).tocsc()
# AA_R = sp.lil_matrix(AA[2]).tocsc()
#
# while (Tol > 10 **(-4)) & (compt < 5000):
#     compt = compt + 1
#     # plt.clf()
#     Vect_1 = la.spsolve(AA_L, Vect)
#     # plt.title('Vect_1')
#     # plt.plot(Vect_1)
#     # plt.grid()
#     # plt.pause(dt)
#
#     # plt.clf()
#     Vect = la.spsolve(AA_R, Vect_1)
#     # plt.title('Vect')
#     # plt.plot(Vect)
#     # plt.grid()
#     # plt.pause(dt)
#     # % Vect = Vect / L2norm(Vect);
#     Vect = Vect / (dz * np.sum(Vect))
#     # print (Tol)
#     Tol = norm_n(Vect - Vold, dz, 1)/norm_n(Vect, dz, 1)
#     # print(numerix.L2norm(Vect - Vold)/numerix.L2norm(Vect))
#     Vold = Vect
#
# compt = compt
# Tol = Tol
#
# Vect_1 = la.spsolve(AA_L, Vect)
# Vect = la.spsolve(AA_R, Vect_1)
# Lambda = Vect.dot(Vold)/(Vold.dot(Vold))
# Lambda = 1. / Lambda
# Lambda = Coef - Lambda