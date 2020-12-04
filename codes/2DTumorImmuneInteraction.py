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


# ----------------------------- The computational domain for the immune cells displacement --------------------------

cellSize = .1
radius = 1.


def plus(z):
    return 0.5*(z+np.abs(z))


def minus(z):
    return 0.5*(-z+np.abs(z))

def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)

# ---------------------------The immune cells displacement domain----------------------------
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

xt, yt = mesh.cellCenters
nVol = mesh.numberOfCells


# ----------------------Configuration the test cases and the results saving --------------------
a = .8
V = .616
delta = 1.
pf = 0.25
S = 20.
chi = 8.64e-1
D = 0.025
sF = 200.e-5
gF = 0.18
d = 1.
init_value = 5.
alpha = 0.5
beta = 0.
tested_param_name = 'a' # set the name of the parameter to be tested for the data file
params = TF.TestFeature(a, V, delta, pf, S, chi, D, sF, gF, d, tested_param_name, init_value, alpha, beta)


# Adjust the tumor size computational domain

mesh_temp = 0.
n_max = 0
if params.a > params.V:
    n_max = 10000
    Lz = 5.
    dz = Lz/n_max
    mesh_temp = Grid1D(dz, n_max)
if params.a < params.V:
    n_max = 7000
    Lz = 10.
    dz = Lz / n_max
    mesh_temp = Grid1D(dz, n_max)

nz = n_max
print ('grid size: ' + str(nz))
meshz = mesh_temp
z, = meshz.cellCenters
z_max = np.max(z)


# Compute the stationary tumor cells profile
t_bar = 0.
for n in np.arange(20):
    t_bar = t_bar + ((-1) ** n) * alphan(n) * numerix.exp(
        -(2. ** (n + 1)) * params.a * z / params.V)
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


# ----------------Variables and unknowns for the Tumor growth equation--------------------
# const = 0.1
Tau = 1.
T = CellVariable(name="$T(t,z)$", mesh=meshz, value=0.0, hasOld=1)
T.setValue(5.*(4./(z_max**2)), where=z < z_max/2.)
T.setValue(5./((z_max/3.)**2 - 0.125**2), where=(z > 0.125) & (z < z_max/3.))# pour a=16

# the tumor cells growth rate
V = CellVariable(name="$V$", mesh=meshz, value=params.V)
V1 = V.harmonicFaceValue.value

# the tumor cells division rate
a = CellVariable(name="$a$", mesh=meshz, value=params.a)

# --------------- Convection Matrix construction for the tumor growth ---------------
vPlus = plus(V1)
vMinus = minus(V1)
# conVT = sp.lil_matrix(sp.spdiags([-vPlus[1:nz+1], (vPlus[1:nz+1] + vMinus[0:nz]), -vMinus[1:nz+1]], [-1, 0, 1], nz, nz))
GP = vPlus
GM = vMinus
GP2 = np.concatenate([[0], vPlus[1:nz]])
GM2 = np.concatenate([[0], vMinus[1:nz]])
conVT = sp.lil_matrix(sp.spdiags([-GP[1:nz + 1], GP[1:nz + 1] + GM[0:nz], -GM2], [-1, 0, 1], nz, nz))


# --------------------- the ChemoAttractant diffusion equation ------------------

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

# ------------------------------------------ The Chemical Potential Diffusion Matrix ------------------------------
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


# ------------------------------------------ THE IMMUNE CELLS ------------------------------

# ----------------------------------------- The diffusion matrix -----------------------------------------

aDiffE = np.zeros((nVol, nVol))
aDiffE = sp.csc_matrix(aDiffE)
aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

# Test for the Dirichlet Laplacian
# la.spsolve(-aDiffE, mesK*(2*np.pi*np.sin((mo(x)**2)*np.pi/2.) + (np.pi**2)*(mo(x)**2)*np.cos((mo(x)**2)*np.pi/2.)))
# the solution should be (np.cos((mo(x)**2)*np.pi/2.))

# -----------------------------------Neumann Boundary condition------------------------------------------
# aDiffE = aDiffP + sp.coo_matrix((0.*TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

# e = np.ones((1, nVol))
# EaDiffE = sp.csc_matrix(np.concatenate((np.concatenate((d*aDiffE.T.todense(), e.T), axis=1),
#                                         np.array([np.append(e.T, 0.)])), axis=0))

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

# ---------------Variables and unknwons for the Immune Cells Displacement equation---------
E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)

# -------------------- Homogeneous distribution of naive immune cells ------------------------------------------
S = CellVariable(mesh=mesh, value=0.0, hasOld=1)
S.setValue(params.S / np.sum(mesK))

# -------------------- Heterogeneous distribution of naive immune cells ------------------------------------------
# center_x1 = 0.207
# center_x2 = 0.907
# Sh.setValue(params.S * numerix.exp(-(((xt - 0.) ** 2)/0.5 +((yt - 1.) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt - center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125) +
#            params.S * numerix.exp(-(((xt + center_x1) ** 2)/0.25 +((yt + center_x2) ** 2))/0.125))
# Sh.setValue(Sh/numerix.sum(mesK*Sh))

delta = CellVariable(name="$\delta(x,y)$", mesh=mesh, value=0.)
delta.setValue(delt.GaussianImmuneResponse2D(params.delta, xt, yt, Ra=0.02))

# pf = pF.GaussianConversionRate2D(amppf, xt, yt, mean=0, Rp=1.)
pf = params.pf * np.ones(nVol)
pfh = 0.0*params.pf * np.ones(nVol)
# gF = gammaF.gaussianGammaF2D(ampgF, xt, yt, Rs=0.3)
gF = params.gf * np.ones(nVol)


Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# -------------------------Solutions representation -----------------------------------
# fig, axes = plt.subplots(1, 3)
# cmap = plt.get_cmap("Spectral_r")
# Tfigure = Matplotlib1DViewer(vars=T, axes=axes[0], datamin=0., datamax=10., figaspect='auto', cmap=cmap)
# Efigure = Matplotlib2DViewer(vars=E, axes=axes[1], figaspect='auto', cmap=cmap)
# Phifigure = MatplotlibStreamViewer(vars=phi.grad, axes=axes[2], figaspect='auto', cmap=cmap)
# axes[0].set_title('The Tumor $T$',fontsize=10, fontweight='bold')
# axes[1].set_title('The Immune Cells',fontsize=10, fontweight='bold')
# axes[2].set_title('The chemo. pot.',fontsize=10, fontweight='bold')
# viewers = MultiViewer(viewers=(Tfigure, Efigure, Phifigure))


# ----------------------Configuration the test case results saving --------------------
test_case_results = np.empty([], dtype=[(params.tested_param, np.float64),
                                        ('t', np.float64),
                                        ('dt', np.float64),
                                        ('dz', np.float64),
                                        ('T', np.float64, (nz,)),
                                        ('E', np.float64, (nVol,)),
                                        ('phi', np.float64, (nVol,))])
test_case_results = np.delete(test_case_results, 0)

# typeOfS = 'homo'
typeOfS = 'hetero'


# -------------------------Computation of the solutions in time-------------------------
CFL = .8

# The time step settings
dt = min(5e-4, CFL * dz / V.value[0]);
tt = 0.
tps = np.arange(0., 1, 1)
Tf = 200.
n = 0

Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# compute the number of cells in the tumor
mu0_new = dz * numerix.sum(T)
# compute the volume of the tumor
mu1_new = dz * numerix.sum(z * T)
# compute the number of Immune cells
mu0E_new = numerix.sum(mesK * E)
mu1E_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * E)

# Saving the moments at time t0
# mu1s.append(mu1_new.value)
# mu1Es.append(mu1E_new.value)

eps = 1e-10
while tt < Tf:
    if n % 10 == 0:
        test_case_results = np.append(test_case_results, np.asarray((params.getAttributeByName(params.tested_param),
                                                                     tt,
                                                                     dt,
                                                                     dz,
                                                                     T.value,
                                                                     E.value,
                                                                     phi.value), dtype=test_case_results.dtype))
    # viewers.plot()
    # -----------------Solve the Immune Cells displacement equation---------------------
    bE = ((mesK / dt) * E.value + mesK * mu1_new.value * pf * S.value).T
    E_new = la.spsolve((Id.multiply(mesK / dt) - D * aDiffE.tocsc() - params.chi * mu1_new.value * aConVE + Id.multiply(mesK * gF )),
                       bE)

    E.setValue(E_new)
    E.updateOld()
    # -------------------------Solve the Tumor Growth equation---------------------------
    Strength = numerix.sum(mesK * delta * E)

    AS = sp.lil_matrix((sp.spdiags(a.value + Strength.value, [0], nz, nz)))

    # Implicit scheme for the tumor growth
    BT = T.value + 4. * dt * dz * G.dot(a.value * T.value)
    T_new = sp.linalg.spsolve((Idz + (dt / dz) * conVT.tocsc() + (dt / dz) * (dz * AS.tocsc())), BT)

    T.setValue(T_new)
    T.updateOld()

    # Explicit scheme for the tumor growth
    # T_new = (Idz - (dt / dz) * (conVT.tocsc() + dz * AS.tocsc())).dot(T.value) + 4 * dt * dz * G.dot(
    #     a.value * T.value)
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

    # compute the volume of the tumor
    mu1_new = dz * numerix.sum(z * T)

    # Setting the CFL Condition
    Flux = params.chi * aConVE

    # dtE = 5e-4
    dtE = 5e-2

    dtT = CFL * dz / np.max(np.abs(V.value))
    # dtT = 5e-2
    # dtT = 3e-3
    if n % 100 == 0:
        print('time: ', tt, 'etape: ', n)
        print ("dtE: {0} dtT: {1} ".format(dtE, dtT))
        print ('mu1 = {0} mu1E = {1} strength = {2}'.format(mu1_new, mu1E_new, Strength.value))
    dt = min(dtT, dtE)

    tt = tt + dt
    # viewers.plot()
    n = n + 1


# -----------------------------------------Save the test results ---------------------------------------------

dirpath = 'TestResults/' + typeOfS + 'geneousS/parameter' + params.tested_param
filename = '2D' + typeOfS + 'SCParams' + params.tested_param + str(params.getAttributeByName(params.tested_param))
dataPath = dirpath + '/' + filename + ".dat"
if not os.path.exists(dirpath):
    os.makedirs(dirpath)
with open(dataPath, "wb") as fi:
    # global test_case_results
    print ('... saving ' + params.tested_param + ' in ' + filename)
    # np.save(fi, TestCaseParam, test_case_results)
    np.savez_compressed(fi, test_case_results)
    fi.close()

# ------------------------------------------- Results visualization  -----------------------------------------
test_case_results = np.load(dataPath)
imStrength = np.sum(mesK*delta.value*test_case_results['arr_0.npy']['E'], axis=1)
mu1 = np.sum(dz*z*test_case_results['T'], axis=1)
ts = test_case_results['t']

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams.update({'font.size': 18})


# n=0
n=n+1
fig, ax1 = plt.subplots()
size_of_fig=8
fig.set_size_inches(size_of_fig+2, size_of_fig)
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

ax2.plot(ts, params.a*np.ones(ts.shape[0]), '--', marker='o', markersize=0, color=color2, linewidth=4)

ax2.tick_params(labelcolor=color2)
fig.tight_layout()
