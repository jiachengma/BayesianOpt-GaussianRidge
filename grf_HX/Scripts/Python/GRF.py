import autograd.numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as st
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from pymanopt.manifolds import Stiefel
from pymanopt.optimizers import TrustRegions
from pymanopt import Problem
import pymanopt
from scipy.linalg import cholesky, solve_triangular, cho_solve

class GRF(object):
    """
    Gaussian ridge function for reduced dimension Gaussian Process
    """
    def __init__(self, X_train, X_test, y_train, y_test, m_ridge, n_restart=20, tol=1e-3):
        """
        X_train    -- input train data
        X_test     -- input test data
        Y_train    -- output train data
        Y_test     -- output test data
        m_ridge    -- ridge function input dimension
        n_restart  -- number of times to restart fitting and pick the model with the lowest objective function value
        tol        -- error tolerance of cost to stop iteration
        """
        self.X_train = X_train
        self.X_test = X_test
        self.y_train = y_train
        self.y_test = y_test
        self.m_ridge = m_ridge
        self.n_restart = n_restart
        self.tol = tol
        dim = X_train.shape[1] # original dimension
        self.manifold = Stiefel(dim, m_ridge)
        self.M = np.random.rand(dim, m_ridge) # initialize projection matrix M
        self.kernel = 1.0 * RBF(length_scale=[1 for _ in range(m_ridge)]) + WhiteKernel(noise_level=1.0) # initialize covariance kernel

    def create_cost_fun(self):
        """
        create cost function for optimizing ridge matrix M
        """
        @pymanopt.function.autograd(self.manifold)
        def cost(M):
            U_train = self.X_train @ M
            U_test = self.X_test @ M

            N_train, m = self.X_train.shape

            lengthscales = self.kernel.get_params()['k1__k2__length_scale'] # rbf lengthscale
            sigma2_f = self.kernel.get_params()['k1__k1__constant_value'] # rbf variance
            sigma2_n = self.kernel.get_params()['k2__noise_level'] # noise variance
            L_inv = np.diag(1. / lengthscales.reshape(-1))
            dim = U_train.shape[1] # dimension of ridge function space

            U_train_tilde = U_train @ L_inv 
            # covariance on training data
            G = sigma2_f * np.exp(-0.5*(np.sum(U_train_tilde**2,1).reshape(-1,1) + np.sum(U_train_tilde**2,1) - 
                                       2 * np.dot(U_train_tilde, U_train_tilde.T)))

            G = G + sigma2_n * np.eye(N_train)
            b = np.linalg.solve(G, self.y_train)

            N_test = self.X_test.shape[0]
            U_test_tilde = U_test @ L_inv
            # covariance of testing and training data K(U_test, U_train)
            K_test = sigma2_f * np.exp(-0.5*(np.sum(U_test_tilde**2,1).reshape(-1,1) + np.sum(U_train_tilde**2,1) - 
                                       2 * np.dot(U_test_tilde, U_train_tilde.T)))
            g_test = K_test @ b
            r = 0.5 * np.linalg.norm(self.y_test - g_test) ** 2 / N_test
            return r
        return cost
    
    def pred(self, X_test_pred, return_var=False):
        """ 
        X_test_pred: test points to evaluate ridge function outputs
        return_var: flag to return variance at test points
        
        Return:
        g_test: predictions of posterior mean using ridge function
        var_test: posterior variance at test points if return_var=True
        """
        U_train = self.X_train @ self.M
        U_test = X_test_pred @ self.M

        G = self.kernel(U_train)
        L_ = cholesky(G, lower=True, check_finite=False) # lower triangular
        b = cho_solve((L_, True), self.y_train, check_finite=False)

        K_test = self.kernel(U_test, U_train) # covariance of testing and training data
        g_test = K_test @ b # predicted posterior mean
        
        if not return_var:
            return g_test
        else:  
            v = solve_triangular(L_, K_test.T, lower=True, check_finite=False)
            var_test = self.kernel.diag(U_test).copy() - np.einsum("ij,ji->i", v.T, v)
            return g_test.squeeze(), var_test.squeeze()
    
    def set_XY(self, X_new, Y_new):
        """
        Update GPR model dataset
        """
        self.X_train = np.vstack((self.X_train, X_new))
        self.y_train = np.hstack((self.y_train, Y_new.squeeze()))

    def grf_fit(self):
        last_r =1e10
        err = np.inf
        d, m = self.M.shape
        n_iter = 0
        
        # re-initialize projection matrix M
        V = np.random.randn(d, m)
        q = np.linalg.qr(V)[0]
        self.M = q.copy()
        
        while err > self.tol:
            self.M_pred = self.M.copy()
            n_iter += 1
            U_train = self.X_train @ self.M_pred
            
            # prior covariance
            ker = 1.0 * RBF(length_scale=[1.0 for _ in range(m)], length_scale_bounds=(1e-5, 1e5)) \
              + WhiteKernel(noise_level=1e-4, noise_level_bounds=(1e-6,1e2)) # noise_level: iid noise variance
            
            gpr = GaussianProcessRegressor(kernel=ker, n_restarts_optimizer=20, alpha=1e-10, normalize_y=False) 
            # alpha: adding to diagonal of covariance matrix to prevent numerical issue during fitting

            gpr.fit(U_train, self.y_train)
            g_test_pred, std_test = gpr.predict(self.X_test@self.M_pred, return_std=True)
            E_var = np.mean(std_test ** 2) # expectation of variance
            self.kernel = gpr.kernel_ # posterior kernel
            my_cost = self.create_cost_fun()
            problem = Problem(manifold=self.manifold, cost=my_cost)
            optimizer = TrustRegions(verbosity=0)
            M_new = optimizer.run(problem).point
        
            r = my_cost(self.M)
            err = np.abs(last_r - r) / last_r
            last_r = r
            self.M = M_new.copy()

        return self.M_pred, gpr, r, n_iter, E_var
    
    def run(self):
        E_var_min = np.inf
        r_min = np.inf
        for i in range(self.n_restart):
            M, gpr, r, n_iter, E_var = self.grf_fit();
            print(f'iteration {i}:, error={r}, E_var={E_var}')
            if r < r_min:
                M_opt = M
                gpr_opt = gpr
                r_min = r
                E_var_min = E_var
        self.M = M_opt
        self.kernel = gpr_opt.kernel_
        return M_opt, gpr_opt, r_min, E_var_min