"""
Project: 
File: stacf_stpacf

Created by Scrat on 04.04.2017
"""


import math
import numpy as np


class Stacf:
    def __init__(self, ts_matrix, wa_matrices, t_lags):
        self._ts_matrix = ts_matrix
        self._wa_matrices = wa_matrices
        self._t_lags = t_lags
        self._s_lags = len(self._wa_matrices)
        self._stacf = None

    @staticmethod
    def _st_cov(ts_matrix, w1, w2, t_lag):
        """
        
        :param ts_matrix: 
        :param w1: 
        :param w2: 
        :param t_lag: 
        :return: 
        """
        max_lags = ts_matrix.shape[0] - t_lag
        gamma = 0
        power_wn = w2.T.dot(w1)

        for t in range(0, max_lags):
            gamma += (power_wn.dot(ts_matrix[[t + t_lag]].T.dot(ts_matrix[[t]]))).trace()

        gamma = gamma / max_lags * ts_matrix.shape[1]
        return gamma

    def _st_acf(self, ts_matrix, wa_matrices, t_lags):
        """
        
        :param ts_matrix: 
        :param wa_matrices: 
        :param t_lags: 
        :return: 
        """
        st_acf = np.zeros((t_lags, len(wa_matrices)))

        cov000 = self._st_cov(ts_matrix, wa_matrices[0], wa_matrices[0], 0)
        for idx, w in enumerate(wa_matrices):
            covss0 = self._st_cov(ts_matrix, w, w, 0)
            for t in range(1, t_lags + 1):
                covs0t = self._st_cov(ts_matrix, w, wa_matrices[0], t)
                st_acf[t - 1, idx] = covs0t / math.sqrt(covss0 * cov000)
        return st_acf

    def estimate(self):
        self._stacf = self._st_acf(self._ts_matrix, self._wa_matrices, self._t_lags)
        return self._stacf

    def get(self):
        return self._stacf


class Stpacf(Stacf):
    def __init__(self, ts_matrix, wa_matrices, t_lags):
        Stacf.__init__(self, ts_matrix, wa_matrices, t_lags)

    def _st__mat(self, t_lag):
        """

        :param t_lag: 
        :return: Matrix
        """
        stmat = np.zeros((self._s_lags, self._s_lags))

        for idx, wm in enumerate(self._wa_matrices):
            for jdx, wn in enumerate(self._wa_matrices):
                stmat[idx, jdx] = self._st_cov(self._ts_matrix, wm, wn, t_lag)
        return stmat

    def _st_mat(self):
        """

        :param t_lags: 
        :return: Matrix
        """
        s_lag = self._s_lags
        t_lag = self._t_lags

        slideye = np.eye(t_lag, 2 * t_lag - 1)
        stmat = np.zeros((s_lag * t_lag, s_lag * t_lag))
        for t in range(1, t_lag - 1):
            stmat = stmat + np.kron(slideye[0:t_lag, t:t + t_lag], self._st__mat(t))

        stmat = stmat + np.transpose(stmat)
        stmat += np.kron(np.eye(t_lag, t_lag), self._st__mat(0))
        return stmat

    def _st_vec(self):
        """

        :return: Vector
        """
        st_vec = np.zeros((self._s_lags * self._t_lags))
        for t in range(1, self._t_lags + 1):
            for s, wm in enumerate(self._wa_matrices):
                st_vec[(t - 1) * self._s_lags + s] = self._st_cov(self._ts_matrix, wm, self._wa_matrices[0], t)
        return st_vec

    def _st_pacf(self):
        """

        :return: Matrix
        """
        print('create yule-walker matrix')
        YWmat = self._st_mat()
        print('create yule-walker vector')
        YWvec = self._st_vec()

        s_lag = self._s_lags
        t_lag = self._t_lags

        st_pacf = np.zeros((t_lag, s_lag))

        print('solve yule-walker equitation')
        for t in range(0, t_lag):
            for s in range(0, s_lag):
                index = t * s_lag + s
                sol = np.linalg.solve(YWmat[:index + 1, :index + 1], YWvec[:index + 1])
                st_pacf[t, s] = sol[index]
        return st_pacf

    def estimate(self):
        self._stacf = self._st_pacf()
        return self._stacf
