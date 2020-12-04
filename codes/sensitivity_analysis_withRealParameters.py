#!/usr/bin/env python
# -*- coding: utf-8 -*-

from fipy import *
from numpy import *
import scipy.sparse as sp
import scipy.sparse.linalg as la
import parameterFunctions.immuneResponse as delt
import parameterFunctions.sigmaF as sigmaF
import inspect
from collections import OrderedDict
import pygpc

from pygpc.sobol_saltelli import get_sobol_indices_saltelli
from pygpc.sobol_saltelli import saltelli_sampling


# This function is a modified version of the original pygpc function
def modified_get_sobol_indices(gpc_object, coeffs, n_samples=1e4):
    """
    Calculate the available sobol indices from the gPC coefficients (standard) or by sampling.
    In case of sampling, the Sobol indices are calculated up to second order.

    sobol, sobol_idx, sobol_idx_bool = SGPC.get_sobol_indices(coeffs, algorithm="standard", n_samples=1e4)

    Parameters
    ----------
    coeffs:  ndarray of float [n_basis x n_out]
        GPC coefficients
    algorithm : str, optional, default: "standard"
        Algorithm to determine the Sobol indices
        - "standard": Sobol indices are determined from the gPC coefficients
        - "sampling": Sobol indices are determined from sampling using Saltelli's Sobol sampling sequence [1, 2, 3]
    n_samples : int, optional, default: 1e4
        Number of samples to determine Sobol indices by sampling. The efficient number of samples
        increases to n_samples * (2*dim + 2) in Saltelli's Sobol sampling sequence.

    Returns
    -------
    sobol: ndarray of float [n_sobol x n_out]
        Normalized Sobol indices w.r.t. total variance
    sobol_idx: list of ndarray of int [n_sobol x (n_sobol_included)]
        Parameter combinations in rows of sobol.
    sobol_idx_bool: ndarray of bool [n_sobol x dim]
        Boolean mask which contains unique multi indices.

    Notes
    -----
    .. [1] Sobol, I. M. (2001).  "Global sensitivity indices for nonlinear
           mathematical models and their Monte Carlo estimates."  Mathematics
           and Computers in Simulation, 55(1-3):271-280,
           doi:10.1016/S0378-4754(00)00270-6.
    .. [2] Saltelli, A. (2002).  "Making best use of model evaluations to
           compute sensitivity indices."  Computer Physics Communications,
           145(2):280-297, doi:10.1016/S0010-4655(02)00280-1.
    .. [3] Saltelli, A., P. Annoni, I. Azzini, F. Campolongo, M. Ratto, and
           S. Tarantola (2010).  "Variance based sensitivity analysis of model
           output.  Design and estimator for the total sensitivity index."
           Computer Physics Communications, 181(2):259-270,
           doi:10.1016/j.cpc.2009.09.018.
    """

    if gpc_object.p_matrix is None:
        dim = gpc_object.problem.dim
    else:
        dim = gpc_object.problem_original.dim

    if gpc_object.problem_original is None:
        problem_original = gpc_object.problem
    else:
        problem_original = gpc_object.problem_original

    # generate uniform distributed sobol sequence (parameter space [0, 1])
    coords_norm_01 = saltelli_sampling(n_samples=n_samples, dim=dim, calc_second_order=True)
    coords_norm = zeros(coords_norm_01.shape)

    # transform to respective input pdfs using inverse cdfs
    for i_key, key in enumerate(problem_original.parameters_random.keys()):
        coords_norm[:, i_key] = problem_original.parameters_random[key].icdf(coords_norm_01[:, i_key])

    # run model evaluations
    res = gpc_object.get_approximation(coeffs=coeffs, x=coords_norm)

    # determine sobol indices
    sobol, sobol_idx, sobol_idx_bool = get_sobol_indices_saltelli(y=exp(res),
                                                                  dim=dim,
                                                                  calc_second_order=True,
                                                                  num_resamples=100,
                                                                  conf_level=0.95)

    # sort
    idx = flip(argsort(sobol[:, 0], axis=0))
    sobol = sobol[idx, :]
    sobol_idx = [sobol_idx[i] for i in idx]
    sobol_idx_bool = sobol_idx_bool[idx, :]


    return sobol, sobol_idx, sobol_idx_bool


# @wrap_non_picklable_objects
class MyModel(pygpc.AbstractModel):

    def __init__(self):
        self.fname = inspect.getfile(inspect.currentframe())
        # pass

    # def __reduce__(self):
    #     return (MyModel, (self.fname,))

    def validate(self):
        pass

    def dichotomy(self, mu_a, mu_b, eps, mesK, Q, aDiffE, aConVE, Id, gF, t_s, E, delta_tild):
        while mu_b - mu_a > eps:
            mu = (mu_b + mu_a) / 2.

            bE = (mesK * Q * mu).T
            E_new = la.spsolve((- aDiffE.tocsc() - mu * aConVE + Id.multiply(mesK * gF * t_s)), bE)
            E.setValue(E_new)
            E.updateOld()

            F_mu_m = numerix.sum(mesK * delta_tild.value * E.value) - 1.

            bE = (mesK * Q  * mu_a).T
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

    def not_converge(self, x, y):
        if (abs(x - y) / y) <= 1e-6:
            return True
        else:
            return False


    def simulate(self, process_id=None, matlab_engine=None):
        step = 0
        res = asarray([])
        print(self.p['a'].flatten())
        print("PARAM/LA TAILLE: {0}/{1}".format(self.p["a"], self.p["a"].shape))

        print('HIHIHIHIHIHI  ', float64(self.p["a"]))
        for idx in range(self.p["a"].shape[0]):
            # print(self.p["a"]*self.p['D']*self.p['sF'])
            t_s = 1. / self.p["a"][idx]
            x_s = sqrt(self.p["D"][idx] * t_s)
            c_s = 1. / (t_s * (x_s ** 2) * self.p["delta"][idx])
            # nu = (D/a)*sqrt(D/a)*(S*d*delta)/(chi*sF)
            # mu1_s = mu1_tild
            mu1_s = c_s / (self.p["R"][idx] * t_s)
            # mu0_s = a * mu1_s / V

            phi_s = (x_s ** 2) / (mu1_s * t_s * self.p["chi"][idx])

            # Q = S*mu1_s*t_s/c_s
            Q = 1.
            # print(self.p["sF"] )
            U = (self.p["sF"][idx] * x_s ** 2) / (self.p["K"][idx] * phi_s)
            # print('VOICI LA VALEUR DU PARAMETRE: {0}'.format(self.p["sF"][0]))

            radius = 1. / x_s
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
            # print('je suis ici')
            x = mesh.cellCenters
            xt, yt = mesh.cellCenters
            nVol = mesh.numberOfCells
            nFaces = mesh.numberOfFaces

            intF = mesh.interiorFaceIDs
            extF = arange(0, nFaces, 1)[array(mesh.exteriorFaces)]

            intFacesCells = mesh.faceCellIDs[:, intF]
            extFacesCells = mesh.faceCellIDs[:, extF]

            TKL = mesh._calcFaceAreas() / mesh._calcFaceToCellDistAndVec()[0].sum(axis=0)
            mes_edge = mesh._calcFaceAreas()
            mesK = mesh.cellVolumes
            # ------------------------------------------ The Chemical Potential ------------------------------
            aDiffP = zeros((nVol, nVol))
            aDiffP = sp.csc_matrix(aDiffP)
            aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
            aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
            aDiffP = aDiffP + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
            aDiffP = aDiffP + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

            # -----------------------------------Neumann Boundary condition------------------------------------------
            aDiffP = aDiffP + sp.coo_matrix((0. * TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

            e = ones((1, nVol))
            EaDiffP = sp.csc_matrix(concatenate((concatenate((aDiffP.T.todense(), (mesK * e).T), axis=1),
                                                    array([append((mesK * e).T, 0.)])), axis=0))

            # -----------------------------------Dirichlet Boundary condition------------------------------------------
            test = CellVariable(mesh=mesh, value=0.)
            phi = CellVariable(name="$\phi(t,x,y)$", mesh=mesh, value=0.0, hasOld=1)
            # sF = sigmaF.SigmaF2D(params.sF, xt, yt, Rs=0.05)
            sF = sigmaF.SigmaF2D(1. / x_s, xt, yt, Rs=0.05 / (x_s ** 2))

            F = sF
            extendedF = append(mesK * U * F, 0.)
            phi_new = la.spsolve(EaDiffP, extendedF)
            phi.setValue(phi_new[0:nVol])
            phi.updateOld()

            # ------------------------------------------ The Chemoattractant ------------------------------
            aDiffE = zeros((nVol, nVol))
            aDiffE = sp.csc_matrix(aDiffE)
            aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[0], intFacesCells[0])), shape=(nVol, nVol))
            aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[0], intFacesCells[1])), shape=(nVol, nVol))
            aDiffE = aDiffE + sp.coo_matrix((TKL[intF], (intFacesCells[1], intFacesCells[0])), shape=(nVol, nVol))
            aDiffE = aDiffE + sp.coo_matrix((-TKL[intF], (intFacesCells[1], intFacesCells[1])), shape=(nVol, nVol))

            # -----------------------------------Dirichlet Boundary condition------------------------------------------
            aDiffE = aDiffE + sp.coo_matrix((-TKL[extF], (extFacesCells[0], extFacesCells[0])), shape=(nVol, nVol))

            aConVE = zeros((nVol, nVol))
            aConVE = sp.csc_matrix(aConVE)

            dPhi_int = numerix.dot(phi.faceGrad.value, mesh.faceNormals)[intF]

            aConVE = aConVE + sp.coo_matrix((mes_edge[intF] * plus(dPhi_int), (intFacesCells[0], intFacesCells[0])),
                                            shape=(nVol, nVol))
            aConVE = aConVE + sp.coo_matrix((-mes_edge[intF] * minus(dPhi_int), (intFacesCells[0], intFacesCells[1])),
                                            shape=(nVol, nVol))
            aConVE = aConVE + sp.coo_matrix((-mes_edge[intF] * plus(dPhi_int), (intFacesCells[1], intFacesCells[0])),
                                            shape=(nVol, nVol))
            aConVE = aConVE + sp.coo_matrix((mes_edge[intF] * minus(dPhi_int), (intFacesCells[1], intFacesCells[1])),
                                            shape=(nVol, nVol))

            dPhi_ext = numerix.dot(phi.faceGrad.value, mesh.faceNormals)[extF]

            aConVE = aConVE + sp.coo_matrix((mes_edge[extF] * plus(dPhi_ext), (extFacesCells[0], extFacesCells[0])),
                                            shape=(nVol, nVol))

            Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

            # ---------------Variables and parameters for the Immune Cells Displacement equation---------
            # E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0.53235e6/c_s, hasOld=1)
            E = CellVariable(name="$E(t,x,y)$", mesh=mesh, value=0., hasOld=1)

            delta_tild = CellVariable(name="$\delta_t(x,y)$", mesh=mesh, value=0.)
            delta_tild.setValue(delt.GaussianImmuneResponse2D(1. / x_s, xt, yt, Ra=0.02 / x_s ** 2))

            gF = self.p["gF"][idx]

            Id = sp.spdiags(numerix.ones(nVol), [0], nVol, nVol)

            # ---------------------------------------------- Dichotomie Method --------------------------------------------

            mu_a = 0.
            mu_b = 1.

            F_mu_m = 0.
            F_mu_a = 0.

            eps = 1e-10

            mu = self.dichotomy(mu_a, mu_b, eps, mesK, Q,  aDiffE, aConVE, Id, gF, t_s, E, delta_tild)

            while self.not_converge(mu, mu_b):
                # print(mu, mu_b)
                mu_b = (mu_a + mu_b) / 2.
                mu = self.dichotomy(mu_a, mu_b, eps, mesK, Q, aDiffE, aConVE, Id, gF, t_s, E, delta_tild)
            print('Step -- > {0}'.format(step))
            step = step + 1
            res = append(res, [mu*mu1_s*1e-9])

            print('mu:{0}'.format(mu))
        # res = np.asarray([mu*mu1_s*1e-9])
        # res[:, newaxis]
        res = log(res[:, newaxis])
        return res


def norm_n(V, dx, n):
    if n == 0:
        c_max = max(abs(V))
        yield c_max
    else:
        norme = sum(abs(dx*V)**n)
        yield norme**(1./n)


def mo(x):
    return numerix.L2norm(x)


def plus(z):
    return 0.5*(z+abs(z))


def minus(z):
    return 0.5*(-z+abs(z))

def alphan(n):
    if n == 0:
        return 1
    return 2 * alphan(n - 1.) / ((2. ** n) - 1.)

def toss(deb, fin):
    return random.uniform(deb, fin)

# --------------- Sensitivity Analysis---------------------------
# Create the coffee cup model
# model = un.Model(run=evaluate_mu_un, labels=["tumor volume($mm^3$)"])

model = MyModel()

parameters = OrderedDict()

parameters["a"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[0.1, 0.5])
# parameters["a"] = pygpc.Norm(pdf_shape=[0.2, 0.09])
parameters["D"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[8.64e-5, 1e-3])
# parameters["D"] = [8.64e-5]
parameters["delta"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[1., 60.])
# parameters["R"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[7.573e-8, 1.231e-6])
# parameters["R"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[6.456e-8, 1.520e-6])#IC1
# parameters["R"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[5.5e-7, 1.036e-6])#IC3_99
parameters["R"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[6.11e-7, 9.74e-7])#IC4_95
# parameters["R"] = pygpc.Norm(pdf_shape=[7.923174114490609e-07, 7.945822739100839e-15])
parameters["chi"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[86.4, 86.4e5])
# parameters["chi"] = [86.4]
parameters["sF"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[5e-17, 0.625e-16])
# parameters["sF"] = [5e-17]
parameters["K"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[1e-2, 1.])
# parameters["K"] = [1e-2]
# parameters["gF"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[2e-2, 1.])
parameters["gF"] = pygpc.Beta(pdf_shape=[1, 1], pdf_limits=[2e-2, 1.])

interval = 'IC4_hetero'
# parameters["chi"] = 86.4
# parameters["sF"] = 5e-17
# parameters["K"] = 1e-2
# parameters["gF"] = 2e-2

problem = pygpc.Problem(model=model, parameters=parameters)

# basis = pygpc.Basis()
# basis.init_basis_sgpc(problem=problem,
#                       order=[5, 5, 5],
#                       order_max=15,
#                       order_max_norm=1,
#                       interaction_order=3)

# basis.plot_basis(dims=[0, 1, 2])

#
fn_results = 'Sensitivity_data/PCE_data'.format(interval)
save_session_format = ".hdf5"
# ---------------------------------- Personnalized Options ------------------------------
options = dict()
options["method"] = "reg"
# options["method"] = "quad"
options["solver"] = "Moore-Penrose"
# options["solver"] = "OMP"
options["settings"] = None
options["order"] = [5] * problem.dim # The univariate polynomials expansion orders
options["order_max"] = 5
options["order_max_norm"] = 0.7
# options["order_max_norm"] = 1.
options["interaction_order"] = 2
# options["interaction_order"] = 2
options["matrix_ratio"] = 2
# options["error_type"] = "nrmsd"
options["error_type"] = "loocv"
options["n_samples_validation"] = 1e3
options["n_cpu"] = 2
options["fn_results"] = fn_results
options["save_session_format"] = save_session_format
options["gradient_enhanced"] = False
options["gradient_calculation"] = "FD_1st2nd"
options["gradient_calculation_options"] = {"dx": 0.001, "distance_weight": -2}
options["backend"] = "omp"
# options["grid"] = pygpc.Random
# options["grid"] = pygpc.LHS(parameters_random=problem.parameters_random, seed=1)
options["grid_options"] = None

n_coeffs = pygpc.get_num_coeffs_sparse(order_dim_max=options["order"],
                                       order_glob_max=options["order_max"],
                                       order_inter_max=options["interaction_order"],
                                       dim=problem.dim)
# problem.dim
grid = pygpc.LHS(parameters_random=problem.parameters_random,
                    n_grid=options["matrix_ratio"] * n_coeffs,
                    seed=1)
# grid = pygpc.Random(parameters_random=problem.parameters_random,
#                     n_grid=options["matrix_ratio"] * n_coeffs,
#                     seed=1)
# print('taille grille', grid.n_grid)
# options["fn_results"] = 'Sensitivity_data/PCE_data_{0}'.format(grid.n_grid)

algorithm = pygpc.Static(problem=problem, options=options, grid=grid)
#
# gpc, coeffs, results = algorithm.run()
session = pygpc.Session(algorithm=algorithm)
# # session.grid = algorithm.grid
#
# # #
# # # # # run gPC session
session, coeffs, results = session.run()


dataPath = 'Sensitivity_data/Pygpc_Sobol_idx.txt'.format(interval)

outF = open(dataPath, "w")

mean = session.gpc[0].get_mean(coeffs)
# outF.write('Mean: '+mean)
print("Mean: {}".format(mean))

std = session.gpc[0].get_std(coeffs)
# outF.write('Std: '+std)
print("Std: {}".format(std))

sobol, sobol_idx, sobol_idx_bool = modified_get_sobol_indices(session.gpc[0], coeffs, n_samples=10)


n_idx = len(sobol_idx)
for i in range(n_idx):
    print("Parameter x{}: {}".format(sobol_idx[i]+1, sobol[i][0]))
    str_tmp = ''
    for k in range(problem.dim):
        if len(sobol_idx[i])==k+1:
            if k+1==1:
                str_tmp = str(sobol_idx[i][k] + 1)
            elif k+1>1:
                for m in range(k):
                    str_tmp = str_tmp + str(sobol_idx[i][m] + 1)+' '
                str_tmp = str_tmp + str(sobol_idx[i][k] + 1)
            outF.write(str_tmp +','+str(sobol[i][0]))
            outF.write('\n')


print(sobol_idx_bool)

outF.close()


pygpc.validate_gpc_plot(session=session,
                        coeffs=coeffs,
                        random_vars=["a", "delta"],
                        n_grid=[25, 25],
                        output_idx=0,
                        fn_out=session.fn_results+'plot',
                        folder="gpc_vs_original_plot",
                        n_cpu=options["n_cpu"])

# Validate gPC approximation vs original model function using Monte Carlo simulation
nrmsd = pygpc.validate_gpc_mc(session=session,
                              coeffs=coeffs,
                              n_samples=1e3,
                              fn_out=session.fn_results+'mc',
                              n_cpu=options["n_cpu"])


