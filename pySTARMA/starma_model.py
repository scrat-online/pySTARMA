"""
Project: 
File: STARIMA

Created by Scrat on 02.03.2017
"""

import numpy as np
from prettytable import PrettyTable

from pySTARMA import utils
from pySTARMA.utils import set_stationary


class STARMA:
    def __init__(self, p, q, ts_matrix, wa_matrices, iterations=2):
        """
        Initialising object/Instance of type STARMA
        :param p: Number or list of auto-regressive-parameters
        :param q: Number or list of moving-average-parameters
        :param ts_matrix: Time series matrix
        :param wa_matrices: List of adjacency matrices
        :param iterations: Number of iterations for kalman filter
        :param cls_name: Name of class or user specified model name
        """
        self._p = p
        self._q = q
        self._wa_matrices = wa_matrices
        self._ts_matrix = ts_matrix
        self._iter = iterations
        self._max_p_tlag = self._get_max_tlag(p)
        self._max_q_tlag = self._get_max_tlag(q)
        self._max_tlag = max(self._max_p_tlag, self._max_q_tlag)
        self._model = None

    def __str__(self):
        if self._model is not None:
            return 'No model fitted yet'
        else:
            return 'Object of class STARMA/STARIMA ' \
                   '\n\t AR-Orders: %s ' \
                   '\n\t MA-Orders: %s ' \
                   % (self._p, self._q)

    @staticmethod
    def _get_max_tlag(x):
        """
        Get maximum time lag of ar- or ma-parameters
        :param x: Number or list of ar- or ma-parameters 
        :return: maximum time lag
        """
        if type(x) is list:
            return max(x) + 1
        else:
            return x

    @staticmethod
    def _get_order_matrix(tlag, slag):
        """
        Generates an matrix containing the indices of the parameters to estimate
        :param tlag: time lag
        :param slag: spatial lag
        :return: matrix with time and spatial orders
        """
        iterate = 0
        if type(tlag) is list:
            order_matrix = np.empty([2, len(tlag) * slag], dtype=int)
            for i, t in enumerate(tlag):
                for s in range(0, slag):
                    order_matrix[0, iterate] = t
                    order_matrix[1, iterate] = s
                    iterate += 1
        else:
            order_matrix = np.empty([2, tlag * slag], dtype=int)
            for t in range(tlag):
                for s in range(slag):
                    order_matrix[0, iterate] = t
                    order_matrix[1, iterate] = s
                    iterate += 1

        return order_matrix

    def _get_ma_matrix(self):
        """
        Get matrix of ma-order with indices of time-lag and spatial lag
        :return: matrix of ma-order-indices
        """
        return self._get_order_matrix(self._q, len(self._wa_matrices))

    def _get_ar_matrix(self):
        """
        Get matrix of ar-order with indices of time-lag and spatial lag
        :return: matrix of ar-order-indices
        """
        return self._get_order_matrix(self._p, len(self._wa_matrices))

    def _get_total_parameter(self):
        """
        Get the total count of parameter to estimate
        :return: Number of total parameter
        """
        return len(self._get_ar_matrix()) + len(self._get_ma_matrix())

    def _fit_model(self, ts_matrix):
        """
        Implementation of the Kalman Filter by Cipra & Motykova - Study on Kalman filter in time series anlysis (1987) and Cheysson - starma: Modelling Space Time AutoRegressive Moving Average (STARMA) Processes
        """
        print('Model fitting in progress')
        print(self.__str__())

        # run kalman filter
        # first iteration
        eps = ts_matrix.copy()
        self._model = utils.kalmanfilter_estimation(ts_matrix,
                                                    self._wa_matrices,
                                                    eps,
                                                    self._get_ar_matrix(),
                                                    self._get_ma_matrix(),
                                                    self._max_p_tlag,
                                                    self._max_q_tlag,
                                                    self._max_tlag)

        # if ma orders present, do further iteration
        if self._q > 0:
            count = 0
            while self._iter > count:
                eps[0: self._max_tlag] = utils.residuals_estimation(ts_matrix[0:self._max_tlag],
                                                                    self._wa_matrices,
                                                                    self._model['phi'],
                                                                    self._model['theta'])

                self._model = utils.kalmanfilter_estimation(ts_matrix,
                                                            self._wa_matrices,
                                                            eps,
                                                            self._get_ar_matrix(),
                                                            self._get_ma_matrix(),
                                                            self._max_p_tlag,
                                                            self._max_q_tlag,
                                                            self._max_tlag)

                count += 1

        # write information to model
        self._model['timeseries'] = ts_matrix
        self._model['residuals'] = utils.residuals_estimation(ts_matrix, self._wa_matrices, self._model['phi'],
                                                              self._model['theta'])
        self._model['sigma2'] = np.trace(self._model['sigma2VarianceMatrix']) / len(self._model['sigma2VarianceMatrix'])
        self._model['llh'] = utils.loglikelihood(self._model) + np.log(
            ts_matrix.size * ((self._get_total_parameter()) * len(self._wa_matrices)))
        self._model['bic'] = self._get_total_parameter() * np.log(ts_matrix.size) - 2 * self._model['llh']
        self._model['phi_tvalue'] = self._model['phi'] / self._model['phi_sd']
        self._model['theta_tvalue'] = self._model['theta'] / self._model['theta_sd']
        self._model['phi_pvalue'] = self._p_value(self._model['phi_tvalue'])
        self._model['theta_pvalue'] = self._p_value(self._model['theta_tvalue'])

    def _p_value(self, t_value):
        """
        TODO check if calculation is correct with tvalue = paramter / std and pavalue with df as totalparameter and not
        as total observations - totalparameter
        :param t_value:
        :return: p-value of parameter
        """
        from scipy.stats import t as t_dist
        df = (self._ts_matrix.size) - (self._get_total_parameter() * len(self._wa_matrices))
        # TODO check calculation of p-values with degrees of freedom
        # p_value = tdist.pdf(abs(tvalue), self._total_parameter())
        return t_dist.pdf(abs(t_value), df)

    def fit(self):
        """
        Estimate parameter for model
        """
        self._fit_model(self._ts_matrix)

    def predict(self, ts_matrix, t_lags):
        """
        Estimate forecasting for model
        :param ts_matrix: Time series matrix
        :param t_lags: Maximum time lags in future 
        :return: Matrix with predictions 
        """
        return utils.prediction(ts_matrix
                                , self._wa_matrices
                                , self._model['phi']
                                , self._model['theta']
                                , t_lags)

    def print_results(self):
        """
        Print model results to screen
        """
        if self._model != 0:
            #table = PrettyTable(['coefficients', 'parameter', 'std deviation', 't-value', 'p-value'])
            #for i in range(0, len(self._model['phi'])):
            #        table.add_row(i)
            #for i in range(0, len(self._model['theta'])):
            #        table.add_row(i)
            #print table
            print('sigma2 is estimated as:\t\t %s' % self._model['sigma2'])
            print('standard error is estimated as:\t\t %s' % np.sqrt(self._model['sigma2']))
            print('BIC is estimated as:\t\t %s' % self._model['bic'])
            print('LogLikelihood is estimated as:\t\t %s' % self._model['llh'])

    def get_model(self):
        """
        Get the model
        :return: model
        """
        return self._model

    def get_item(self, item):
        """
        Get specific item of model
        :param item: Key for  item
        :return: Value
        """
        return self._model[item]


class STARIMA(STARMA):
    def __init__(self, p, q, d, ts_matrix, wa_matrices, iterations=2):
        """
        Initialising object/Instance of type STARIMA
        :param p: Number or list of auto-regressive-parameters
        :param q: Number or list of moving-average-parameters
        :param d: Number or list of differencing 
        :param ts_matrix: Time series matrix
        :param wa_matrices: List of adjacency matrices
        :param iterations: Number of iterations for kalman filter
        """
        STARMA.__init__(self, p, q, ts_matrix, wa_matrices, iterations)
        self._ts = ts_matrix
        self._d = d  # Number of Differencing

    def __str__(self):
        return STARMA.__str__(self) + \
               '\n\t Difference: %s ' \
               % (self._d,)

    def fit(self):
        """
        Estimate parameter for model
        """
        self._ts = set_stationary(self._ts, self._d)
        self._fit_model(self._ts)
