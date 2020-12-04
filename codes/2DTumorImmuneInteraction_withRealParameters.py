#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fipy import *
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as la
from scipy import signal as sign
from scipy import optimize as opt
from matplotlib import pyplot as plt
import cma
import parameterFunctions.immuneResponse as delt
import parameterFunctions.sigmaF as sigmaF
import parameterFunctions.TestFeature as TF

def plus(z):
    return 0.5*(z+np.abs(z))


def minus(z):
    return 0.5*(-z+np.abs(z))


def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)


def decay_func(x,a,b,c,d):
    return a*x**3 + b*x**2 + c*x + d

# def decay_func(x,a,b,c,d,e,f):
#     return a*x**5 + b*x**4 + c*x**3 + d*x**2 + e*x + f

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
# tested_param = ''
tested_param_name = 'a' # the parameter to test

t_r = 0.

# mu0_tild = 1.41765e6
# mu1_tild = 24.96e9  # * (u.micrometer**3)

mu0_tild = 1.
mu1_tild = 4188.  # * (u.micrometer**3)

# mu0_tild = 1.
# mu1_tild = 4189 # * (u.micrometer**3)

# mu0_tild = 0.00028
# mu1_tild = 5.  # * (u.micrometer**3)

params = TF.TestFeature(a, V, delta, pf, S, chi, D, sF, gF, d, tested_param_name, init_value, alpha, beta)

# ------------------------------------ Dimensionless parameters ---------------------------
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
cellSize = radius/20.

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

# c_s = t_s*S*mu1_s


# --------------------------------- Diffusion Matrix construction ---------------------------
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

# --------------------------------------- THE CHEMICAL POTENTIAL --------------------------------

# Diffusion matrix
aDiffP = np.zeros((nVol, nVol))
aDiffP = sp.csc_matrix(aDiffP)
aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

# ----------------------------------- Neumann Boundary condition ------------------------------------------
aDiffP = aDiffP + sp.coo_matrix((0.*TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

e = np.ones((1, nVol))
EaDiffP = sp.csc_matrix(np.concatenate((np.concatenate((d*aDiffP.T.todense(), (mesK*e).T), axis=1),
                                        np.array([np.append((mesK*e).T, 0.)])), axis=0))

# ----------------------------------- Dirichlet Boundary condition ------------------------------------------
test = CellVariable(mesh=mesh, value = 0.)
phi = CellVariable(name="$\phi(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
# sF = sigmaF.SigmaF2D(params.sF, xt, yt, Rs=0.05)
sF = sigmaF.SigmaF2D(1./x_s, xt, yt, Rs=0.05/(x_s**2))

# ------------------------------------ Solution of the chemical potential equation ------------------------------------
F = sF
extendedF = np.append(mesK * U * F, 0.)
phi_new = la.spsolve(EaDiffP,  extendedF)
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

# ---------------Variables and unknowns for the Immune Cells Displacement equation---------
E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0., hasOld=1)

# delta_tild represents the dimensionless version of delta

delta_tild = CellVariable(name="$\delta_t(x,y)$", mesh=mesh, value=0.)
delta_tild.setValue(delt.GaussianImmuneResponse2D(1./x_s, xt, yt, Ra=0.02/x_s**2))


# pf = pF.GaussianConversionRate2D(amppf, xt, yt, mean=0, Rp=1.)
pf = params.pf * np.ones(nVol)
# pfh = 0.0*params.pf * np.ones(nVol)
# gF = gammaF.gaussianGammaF2D(ampgF, xt, yt, Rs=0.3)
gF = params.gf * np.ones(nVol)

Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# -------------------------Solutions representation -----------------------------------
# fig, axes = plt.subplots(1, 3)
# cmap = plt.get_cmap("Spectral_r")
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
                                        ('mu0', np.float64),
                                        ('mu1', np.float64),
                                        ('E', np.float64, (nVol,)),
                                        ('phi', np.float64, (nVol,))])
test_case_results = np.delete(test_case_results, 0)


# -------------------------Computation of the solutions in time-------------------------
CFL = .8
# The time step
dt = 0.001
tt = 0.
tps = np.arange(0., 1, 1)
Tf = 500.
n = 0

Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

# compute the number of cells in the tumor
mu0 = mu0_tild/mu0_s
# mu0 = 1.41765
# compute the volume of the tumor
mu1 = mu1_tild/mu1_s
# compute the number of Immune cells
mu0E = numerix.sum(mesK * E)
mu1E = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * E)


# Saving the moments at time t0
# mu0s = []
# mu1s = []
# mu1s.append(mu1_new.value)
# mu1Es.append(mu1E_new.value)


eps = 1e-10
while tt < Tf:
    # print('time: ', tt, 'etape: ', n)
    if dt > 1e-2:
        if n % 10 == 0:
            test_case_results = np.append(test_case_results, np.asarray((params.getAttributeByName(params.tested_param),
                                                                         tt,
                                                                         dt,
                                                                         mu0,
                                                                         mu1,
                                                                         E.value,
                                                                         phi.value), dtype=test_case_results.dtype))

    else :
        if n % 100 == 0:
            test_case_results = np.append(test_case_results, np.asarray((params.getAttributeByName(params.tested_param),
                                                                         tt,
                                                                         dt,
                                                                         mu0,
                                                                         mu1,
                                                                         E.value,
                                                                         phi.value), dtype=test_case_results.dtype))
    # viewers.plot()
    # -----------------Solve the Immune Cells displacement equation---------------------
    bE = ((mesK / dt) * E.value + mesK * Q * mu1).T
    E_new = la.spsolve((Id.multiply(mesK / dt) - aDiffE.tocsc() - (mu1) * aConVE + Id.multiply(mesK * gF*t_s)),
                       bE)
    E.setValue(E_new[0])
    E.updateOld()
    # -------------------------Solve the Tumor Growth equation---------------------------
    Strength = numerix.sum(mesK * delta_tild * E)

    mu1E_new = numerix.sum(mesK * numerix.L2norm(mesh.cellCenters) * E)

    # compute the volume of the tumor
    mu0_new = (1. + dt*(1. - Strength.value))*mu0
    mu1_new = dt*mu0_new + (1. - dt*Strength.value)*mu1

    mu0 = mu0_new
    mu1 = mu1_new

    # Setting the CFL Condition
    Flux = aConVE

    dtE = CFL * np.min(mesK) / (np.max(np.abs(Flux)));
    dtT = CFL / Strength.value
    # dtT = 5e-2
    # dtT = 3e-3
    if n % 100 == 0:
        # print ("dtE: {0} dtT: {1} ".format(dtE, dtT))
        print('time: ', tt, 'etape: ', n)
        print ('mu1 = {0} mu1E = {1} strength = {2}'.format(mu1_new, mu1E_new, Strength.value))
        print(dt)
    # dt = min(dtT, dtE, dt_tmp)
    dt = min(5e-2,dtT, dtE)

    tt = tt + dt
    # viewers.plot()
    n = n + 1

# -----------------------------------------Save the test results ---------------------------------------------

dataPath = 'TestResults_withRealParameters/evolution_realparams_a{0}.dat'.format(params.tested_param)
with open(dataPath, "wb") as fi:
    # np.save(fi, TestCaseParam, test_case_results)
    np.savez_compressed(fi, test_case_results)
    fi.close()
data = np.load(dataPath)


# ------------------------------------------- Results visualization  -----------------------------------------

imStrength = (x_s**2)*c_s*params.delta*np.sum(mesK*delta_tild.value*data['arr_0.npy']['E'], axis=1)
mu1s = data['arr_0.npy']['mu1']
ts = data['arr_0.npy']['t']

# time_thr = np.max(ts*t_s)
time_thr = 500
times = (ts*t_s)[(ts*t_s)<=time_thr]
Nt = times.shape[0]
diameters = 2.*np.power(3.*(mu1s[0:Nt]*mu1_s*1e-9)/(4.*np.pi), 1./3.)

# --------------------- We compute the peaks of the tumor mass oscillations---------------------
peaks = sign.find_peaks(diameters)
n_peaks = peaks[0].shape[0]

# --------------------- We use the cma-es algorithm to fit the peaks of the tumor mass oscillations---------------------
# def obj_func(x):
#     f = 0.
#     for i in range(n_peaks):
#         f = f + ((x[0]*np.exp(-x[1]*(times[peaks[0][i]] + x[2])) + x[3]) - diameters[peaks[0][i]])**2
#     return f

# def obj_func(x):
#     f = 0.
#     for i in range(n_peaks):
#         f = f + ((x[0]*times[peaks[0][i]]**3 + x[1]*times[peaks[0][i]]**2 + x[2]*times[peaks[0][i]] + x[3]) - diameters[peaks[0][i]])**2
#     return f

def obj_func(x):
    f = 0.
    for i in range(n_peaks):
        f = f + ((x[0]*times[peaks[0][i]]**5 + x[1]*times[peaks[0][i]]**4 + x[2]*times[peaks[0][i]]**3
                  + x[3]*times[peaks[0][i]]**2 + x[4]*times[peaks[0][i]] + x[5]) - diameters[peaks[0][i]])**2
    return f


es = cma.CMAEvolutionStrategy(6*[0.0], 0.5)
# es = cma.CMAEvolutionStrategy(4*[0.0], 0.5)
es.optimize(obj_func)
# mu1s_grad = np.gradient(mu1s[0:Nt])

# Visualization
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['font.sans-serif'] = ['Tahoma']
plt.rcParams.update({'font.size': 22})

fig, ax1 = plt.subplots()
size_of_fig=8
fig.set_size_inches(size_of_fig+2, size_of_fig)
color1 = 'tab:red'
ax1.set_xlabel('t in $day$')
ax1.set_ylabel('tumor size in $mm$', color='k')
# ax1.set_ylabel('$\mu_1$ in $mm^3$', color=color1)
# ax1.text(np.max(times) - 20.*(times[-1] - times[-2]), diameters[Nt-1]+0.00001, 'tum. diam.')
ax1.plot(times, diameters, '-', color='k', linewidth=4, label='tum. size')
# ax1.plot((ts*t_s + t_r)/secs, mu1*mu1_s*1e-9, '-', color=color1, linewidth=2, label='$t \mapsto \mu_1(t)$')
ax1.legend(loc="lower right")

# ----------------------------------------------- Peaks and limit sup -------------------------------------------
# ax1.scatter(times[peaks[0]], diameters[peaks[0]])
p1, p2 = opt.curve_fit(decay_func, times[peaks[0]], diameters[peaks[0]], method='trf', absolute_sigma=True)
ax1.plot(times[np.logical_and(times>=times[peaks[0]][0], times<=(times[peaks[0]][-1]+10))]
         , decay_func(times[np.logical_and(times>=times[peaks[0]][0], times<=(times[peaks[0]][-1]+10))], *p1)
         , '-'
         , color='black'
         , linewidth=2)

# -----------------------------------------------------------------------------------------------------------------

# ax1.plot(times[times>=times[peaks[0]][0]], decay_func(times[times>=times[peaks[0]][0]], *p1), '-', color='black', linewidth=2)
# ax1.plot(times[peaks[0][0]:Nt], decay_func(times[peaks[0][0]:Nt], *es.result[0]), '-', color='black', linewidth=2)
ax1.tick_params(axis='y', labelcolor='k')

ax2 = ax1.twinx()

color2 = 'tab:blue'
color3 = 'tab:olive'
ax2.set_ylabel('$\overline{\mu_{c}}$', color='k')
ax2.plot(times, imStrength[0:Nt], '--', marker='o', markersize=0, color='k', linewidth=2, label='Im. Strength $t \mapsto \overline{\mu_{c}}$')
# ax2.text(np.max(times) - 20*(times[-1] - times[-2]), params.a+0.001, '$t \mapsto \overline{\mu_{c}}$')
ax2.plot(times, params.a*np.ones(Nt), '--', marker='o', markersize=0, color='black', linewidth=2)

ax2.tick_params(labelcolor='k')
ax2.legend(loc="upper right")
fig.tight_layout()
