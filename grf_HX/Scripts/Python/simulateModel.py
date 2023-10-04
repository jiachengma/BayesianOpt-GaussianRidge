import numpy as np
import matplotlib.pyplot as plt
from dymola.dymola_interface import DymolaInterface
from pyDOE import lhs
import math


class simulateModel(object):
    """
    Run Dymola to simulate transient model
    """
    def __init__(self, model_path, model, startTime, stopTime, outputInterval, lb, ub, outputWeight):
        """
        model_path       -- Modelica library path 
        model            -- Dymola model name
        startTime        -- simulation start time
        stopTime         -- simulation stop time
        outputInterval   -- output data storing interval
        lb               -- list of lower bound of calibration parameters
        ub               -- list of upper bound of calibration parameters
        outputWeight     -- list of output weights in loss 
        """
        self.model_path = model_path
        self.model = model
        self.startTime = startTime
        self.stopTime = stopTime
        self.outputInterval = outputInterval
        self.lb = lb
        self.ub = ub
        self.outputWeight = outputWeight
        self.n_output = len(outputWeight)
        self.dymola = None

    def instantiate_dymola(self):
        self.dymola = DymolaInterface(showwindow=True)
        self.dymola.openModel(path=self.model_path, changeDirectory=False)
        self.dymola.openModel(path='C:\Jiacheng Ma\Modelica libraries\ExternalMedia-4.0.0\ExternalMedia\Package.mo') # ExternalMedia library

    def close_dymola(self):
        if self.dymola is not None:
            self.dymola.close()
            self.dymola = None

    def run_dymola(self, theta_in):
        """
        Run simulation using Dymola-Python interface

        theta_in: list of model parameters

        Outputs:
        y_pred: time series data of model predictions
        y_mea: time series data of measurements
        """
        assert(len(theta_in) == len(self.lb))
        initialNames = ['u[{}]'.format(i) for i in range(1,len(theta_in)+1)]
        initialValues = theta_in
        self.dymola.experimentSetupOutput(derivatives=False, events=False, auxiliaries=False)
        result, finalVar = self.dymola.simulateExtendedModel(problem=self.model,
                                            startTime=self.startTime, 
                                            stopTime=self.stopTime,
                                            outputInterval=self.outputInterval,
                                            method='Dassl',
                                            tolerance=0.0001,
                                            initialNames=initialNames,
                                            initialValues=initialValues)
        if not result:
            print(theta_in)
            print("Simulation failed. Below is the translation log.")
            log = self.dymola.getLastErrorLog()
            print(log)
            return None, None
        else:
            Nrows = self.dymola.readTrajectorySize("dsres.mat")
            outputNames = ['y[{}]'.format(i) for i in range(1,self.n_output+1)] + ['y_mea[{}]'.format(i) for i in range(1,self.n_output+1)]
            outputVar = self.dymola.readTrajectory("dsres.mat", outputNames, Nrows)
            y_pred = np.array(outputVar[:self.n_output])
            y_mea = np.array(outputVar[self.n_output:])
            ner = np.linalg.norm(y_pred[:,10:] - y_mea[:,10:], axis=1) / np.linalg.norm(y_mea[:,10:], axis=1) # omit some initialization points
            W = np.diag(self.outputWeight)
            cost = np.dot(ner.T,W.dot(ner))
            return cost, outputVar

    def J_calib(self, u):
        """
        calibration objective function

        Input: 
        u: scaled parameters
        Output:
        negative log loss to maximize
        """
        u_truescale = np.round(self.lb + u * (self.ub - self.lb),3)
        cost, outputVar = self.run_dymola(list(u_truescale))
        if not cost:
            return math.nan
        else:
            return -np.log(cost)
        
    def generate_samples(self, n, save_path):
        """
        run simulations to generate samples of calibration objective function values
        """
        np.random.seed(123)
        # Generate scaled samples of the input space
        X_normalize = lhs(len(self.lb), n, 'c')
        # Get corresponding results at function space
        Y = np.zeros(n)
        for i in range(n):
            Y[i] = self.J_calib(X_normalize[i,:])
            print(i+1, Y[i])
        X_normalize = X_normalize[~np.isnan(Y),:]
        Y = Y[~np.isnan(Y)]

        # Plot objective funciton values
        fig, ax = plt.subplots()
        ax.plot(Y,'kx',markersize=10, markeredgewidth=2)
        ax.set_xlabel('$n$')
        ax.set_ylabel('$J(u)$')
        ax.grid(alpha=0.7)
        np.savez(save_path, X_normalize=X_normalize, Y=Y, lb=self.lb, ub=self.ub)
        


        


