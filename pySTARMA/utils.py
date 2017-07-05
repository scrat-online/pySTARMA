"""
Project: 
File: diagnostic

Created by Scrat on 04.04.2017
"""
import numpy as np
import math

import pandas as pd
from numpy.linalg import inv


def loglikelihood(model):
    """
    Calculates Log-Likelihood of model
    :param ts_matrix: Matrix of timeseries
    :param model: Estimated model
    :return: log-likelihood
    """
    res = model['residuals']
    ts_cols = len(res[0])
    ts_rows = len(res)
    sigma2 = model['sigma2']

    llh = ts_cols * ts_rows * (np.log(2 * math.pi) + np.log(sigma2))

    for t in range(ts_rows):
        llh += res[[t]].dot(1 / sigma2).dot(res[[t]].T)
    return -llh / 2


def prediction(ts_matrix, wa_matrices, phi, theta, prediction_lag=1):
    """
    Calcuates predictions for time series matrix 
    :param ts_matrix: time series matrix
    :param wa_matrices: list of adjacency matrices
    :param phi: parameters phi
    :param theta: parameters theta
    :param prediction_lag: maximum prediction lag
    :return: 
    """
    predictions = ts_matrix[:ts_matrix.shape[0] - prediction_lag, :].copy()
    residuals = residuals_estimation(predictions, wa_matrices, phi, theta)

    for h in range(prediction_lag):
        predict = 0
        for tlag in range(len(phi)):
            for slag in range(len(phi[0])):
                predict += (predictions[[- tlag - 1]] * phi[tlag, slag]).dot(wa_matrices[slag].T)
        for tlag in range(len(theta)):
            for slag in range(len(theta[0])):
                predict -= (residuals[[- tlag - 1]] * theta[tlag, slag]).dot(wa_matrices[slag].T)
            # TODO implement recursive innovation algorithm to update innovations/errors in prediction
            residuals = np.concatenate((residuals, np.zeros([1, len(residuals[0])])), axis=0)
        predictions = np.concatenate((predictions, predict), axis=0)

    return predictions


def kalmanfilter_estimation(ts_matrix, wa_matrices, residuals, ar_matrix, ma_matrix, p_lag, q_lag, max_t_lag):
    """
    Estimate parameters and variance with kalman filtering. After implementation of the Kalman Filter 
    by Cipra & Motykova - Study on Kalman filter in time series anlysis (1987) and 
    Cheysson - starma: Modelling Space Time AutoRegressive Moving Average (STARMA) Processes
    :param p_lag: maximum auto regressive parameter
    :param max_t_lag: maximum time lag
    :param ts_matrix: time series matrix
    :param wa_matrices: list of adjacency matrices
    :param residuals: calculated residuals
    :param ar_matrix: matrix with indexes for auto regressive parameters to estimate
    :param ma_matrix: matrix with indexes for moving average parameters to estimate
    :param q_lag: maximum moving average parameter
    :return: ar-parameters, ar_std, ma-parameters, ma-std, variance matrix
    """
    # get variables
    ar = len(ar_matrix[0])
    ma = len(ma_matrix[0])
    dim = ar + ma

    # initialise kalman filter
    h = np.zeros([dim, ts_matrix.shape[1]])
    ksi = np.zeros(dim, )
    p = 100000. * np.eye(dim, dim)
    sigma2 = (1 / 100000.) * np.eye(ts_matrix.shape[1], ts_matrix.shape[1])

    # run the filter
    for t in range(max_t_lag, ts_matrix.shape[0]):
        # Update the observation matrix
        # Fill for 'phi', AR parameters
        for it in range(ar):
            weights = wa_matrices[ar_matrix[1, it]]
            h[it] = (ts_matrix[[t - 1 - ar_matrix[0, it]], :]).dot(weights.T)

        # Fill for 'theta', MA parameters
        for it in range(ar, dim):
            weights = wa_matrices[ma_matrix[1, it - ar]]
            h[it] = (residuals[[t - 1 - ma_matrix[0, it - ar]], :]).dot(weights.T)

        # Create
        nm1 = inv(h.T.dot(p).dot(h) + np.eye(ts_matrix.shape[1], ts_matrix.shape[1]))
        nu = ts_matrix[t].T - h.T.dot(ksi)
        # Prediction & update equations all - in -ones
        ksi += p.dot(h).dot(nm1).dot(nu)  # 2.28 Cipra & Motykova 1987
        p -= p.dot(h).dot(nm1).dot(h.T).dot(p)  # 2.29 Cipra & Motykova 1987
        sigma2 = (sigma2 * (t + 1 - max_t_lag) + nu.T * nu) / (t + 2 - max_t_lag)  # 2.30 Cipra & Motykova 1987
        # Estimate the residual
        residuals[[t]] = ts_matrix[[t]] - (ksi.T.dot(h))  # 2.31 Cipra & Motykova 1987

    # Get estimated standard deviation of the parameters
    sd = np.sqrt(np.trace(sigma2) * np.diag(p) / ts_matrix.shape[1])

    # Rename and reshape
    phi = np.zeros([p_lag, len(wa_matrices)])
    phi_sd = np.zeros([p_lag, len(wa_matrices)])
    theta = np.zeros([q_lag, len(wa_matrices)])
    theta_sd = np.zeros([q_lag, len(wa_matrices)])

    for it in range(ar):
        phi[ar_matrix[0, it], ar_matrix[1, it]] = ksi[it];
        phi_sd[ar_matrix[0, it], ar_matrix[1, it]] = sd[it];
        pass

    for it in range(ar, dim):
        theta[ma_matrix[0, it - ar], ma_matrix[1, it - ar]] = ksi[it];
        theta_sd[ma_matrix[0, it - ar], ma_matrix[1, it - ar]] = sd[it];
        pass

    return {'phi': phi,
            'phi_sd': phi_sd,
            'theta': theta,
            'theta_sd': theta_sd,
            'sigma2VarianceMatrix': sigma2, }


def residuals_estimation(ts_matrix, wa_matrices, phi, theta):
    """
    Calculation of residuals for model 
    :param ts_matrix: time series matrix
    :param wa_matrices: list of adjacency matrices
    :param phi: auto regressive parameters
    :param theta: moving average parameters
    :return: residual matrix
    """
    residuals = ts_matrix.copy()
    for t in range(ts_matrix.shape[0]):
        t_lim = min([t, len(phi)])
        for t_lag in range(t_lim):
            for slag in range(0, len(phi[0])):
                weights = wa_matrices[slag]
                residuals[[t]] -= (ts_matrix[[t - t_lag - 1]] * phi[t_lag, slag]).dot(weights.T)

        t_lim = min([t, len(theta)])
        for t_lag in range(t_lim):
            for slag in range(0, len(theta[0])):
                weights = wa_matrices[slag]
                residuals[[t]] -= (residuals[[t - t_lag - 1]] * theta[t_lag, slag]).dot(weights.T)
    return residuals


def set_stationary(ts_matrix, lags):
    """
    Differencing of time series
    :param ts_matrix: time series matrix
    :param lags: list of differencing 
    :return: difference time series matrix
    """
    time_series_matrix = pd.DataFrame(ts_matrix).copy()
    for t_lag in lags:
        time_series_matrix -= time_series_matrix.shift(t_lag)
    time_series_matrix.dropna(inplace=True)
    return time_series_matrix.as_matrix()

