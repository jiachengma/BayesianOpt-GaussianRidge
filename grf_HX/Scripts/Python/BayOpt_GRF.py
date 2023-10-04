import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from GRF import GRF
from mpl_toolkits.mplot3d import Axes3D
import scipy.stats as st

# maximum upper interval
def mui(m, sigma, ymax, psi=1.96):
    return m + psi * sigma

# Expected improvement acquisition function
def ExpectedImprovement(m, sigma, ymax):
  """Return the expected improvement.

  Arguments
  m       -- The predictive mean at the test points.
  sigma   -- The predictive standard deviation at the test points.
  ymax    -- The maximum observed value so far.
  """
  u = (m - ymax) / sigma
  ei = sigma * (u * st.norm.cdf(u) + st.norm.pdf(u))
  ei[sigma <= 0.] = 0.
  return ei

class BayOpt_GRF(object):
    """
    Bayesian optimization with Gaussian ridge function for calibrating transient vapor compression cycle
    """
    def __init__(self, Xdata, Ydata, simulationModel, alpha_fun=mui, alpha_fun_params={}):
        """
        Xdata              -- training samples of calibration parameters
        Ydata              -- training samples of objective function values
        simulationModel    -- Dymola model
        alpha_fun          -- information acquisition function
        alpha_fun_params   -- acquisition function parameters
        """
        self.Xdata = Xdata
        self.Ydata = Ydata
        self.simulationModel = simulationModel
        self.alpha_fun = alpha_fun
        self.alpha_fun_params = alpha_fun_params
        self.Y_normalizer = StandardScaler()
        self.Y_normalizer.fit(Ydata.reshape(-1,1))

    def train_grf(self, m_max, n_restart=20, tol=1e-4, plotResults=1):
        """
        determine ridge space
        m_max          -- maximum dimension to evaluate
        n_restart      -- number of times to repeat fitting to pick the optimal ridge subspace for each dimension
        tol            -- error tolerance for each fitting
        """
        # split data to training and testing 
        Y_scaled = self.Y_normalizer.transform(self.Ydata.reshape(-1,1))
        X_scaled = self.Xdata # Xdata is already scaled
        X_train, X_test, y_train, y_test = train_test_split(X_scaled, Y_scaled, test_size=0.5, random_state=10)
        # run over all possible dimensions to determine ridge space dimension
        grf_all = []
        M_all = []
        gpr_all = []
        r_all = []
        E_var_all = []
        for i in range(1, m_max+1):
            grf_i = GRF(X_train, X_test, y_train[:,0], y_test[:,0], i, n_restart, tol)
            M_final, gpr_final, r_final, E_var_final = grf_i.run()
            grf_all.append(grf_i)
            M_all.append(M_final)
            gpr_all.append(gpr_final)
            r_all.append(r_final)
            E_var_all.append(E_var_final)
            print(f'dimension m={i} prediction error r={r_final}, expectation of posterior variance on test data: {E_var_final}')
            print('-' * 50)
        
        # plot results for all dimensions
        if plotResults:
            # error and variance expectation for test data
            fig, ax = plt.subplots(figsize=(19.2,14.4), dpi=100)
            x = range(1, len(r_all)+1)
            ln1 = ax.plot(x, r_all, '-s', markersize=10, label='Prediction error $r$')
            ax.set_ylabel('$r$', fontsize=18, fontweight='bold')
            ax.set_xlabel('Ridge subspace dimension', fontsize=15, fontweight='bold')
            plt.xticks(x, fontsize=12, fontweight='bold')
            plt.yticks(fontsize=12, fontweight='bold')
            ax.grid(alpha=0.8)
            ax1 = ax.twinx()
            ax1.set_ylabel('$\mathbb{E}[V]$', fontsize=18, fontweight='bold')
            ln2 = ax1.plot(x, E_var_all, 'r-o', markersize=8, label='Variance expectation $\mathbb{E}[V]$')
            lns = ln1+ln2
            labs = [l.get_label() for l in lns]
            ax.legend(lns, labs, loc=9, fontsize=14)
            plt.yticks(fontsize=12, fontweight='bold')

        return grf_all, M_all, gpr_all, r_all, E_var_all

    def runBO(self, gpr, X_design, n_iter, plotResults=1):
        """
        run Bayesian optimization to optimize (maximize) the objective function

        gpr: The final Gaussian process regression model chosen based on Gaussian ridge dimension
        X_design: The set of candidate points
        n_iter: Nmber of iterations to run BO
        """
        af_all = [] # values of acquisition function
        x_all = []
        y_all = []
        self.simulationModel.instantiate_dymola()
        for count in range(n_iter):
            g, var = gpr.pred(X_design, return_var=True)
            af_values = self.alpha_fun(g, np.sqrt(var.squeeze()), gpr.y_train.max(), **self.alpha_fun_params)
            i_opt = np.argmax(af_values)
            x_new = X_design[i_opt]
            y_new = self.simulationModel.J_calib(x_new)

            if not y_new:
                X_design = np.delete(X_design, i_opt, axis=0)
                print(count+1, 'failed', i_opt)
            else:
                x_all.append(x_new)
                y_all.append(y_new)
                af_all.append(af_values[i_opt])
                # update gpr
                y_new_scaled = self.Y_normalizer.transform(y_new.reshape(-1,1))
                gpr.set_XY(x_new, y_new_scaled.squeeze())
                print(count+1, y_new)
        if plotResults:
            fig, ax = plt.subplots(figsize=(12.8, 9.6))
            ax.plot(range(1, n_iter+1), y_all, '-*', markersize=14, markeredgewidth=2, linewidth=3)
            # ax.set_xticks(range(1,n_iter+1,5))
            ax.set_xlabel('Iterations', fontsize=30, fontweight='bold')
            ax.set_ylabel('$f(x)$', fontsize=30, fontweight='bold')
            ax.grid(alpha=0.7)
        self.simulationModel.close_dymola()
        return gpr, af_all, x_all, y_all






