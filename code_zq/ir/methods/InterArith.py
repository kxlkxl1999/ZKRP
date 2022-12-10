import numpy as np
from scipy.optimize import minimize


# Minimum weighted L2 distance
class WL2:
    def __init__(self, xc, xr, yc, yr):
        """
        # Inputs:
        # xc: the center of input interval Xs
        # yc: the center of output interval Ys
        # xr: the range of input interval Xs
        # yr: the range of output interval Ys
        """
        assert xc.shape[0] == yc.shape[0] and xr.shape[0] == yr.shape[0]
        self.xc = xc
        self.yc = yc
        self.xr = xr
        self.yr = yr

    def fit(self, theta):
        def wl2_obj(para, xc, xr, yc, yr, theta):
            a_c = para[0]
            b_c = para[1]
            a_r = para[2]
            b_r = para[3]

            yc_hat = a_c + b_c * xc
            yr_hat = a_r + b_r * xr
            return np.mean((yc_hat - yc) ** 2 + theta * (yr_hat - yr) ** 2)

        result = minimize(
            wl2_obj,
            x0=np.array([1, 1, 1, 1]),
            args=(self.xc, self.xr, self.yc, self.yr, theta),
            method="Nelder-Mead",
        )
        self.centerbeta = np.array([result["x"][0], result["x"][1]])
        self.rangebeta = np.array([result["x"][2], result["x"][3]])
        return self.centerbeta, self.rangebeta

    def predict(self, xl, xu, cr=False):
        if not cr:
            xc = (xl + xu) / 2
            xr = (xu - xl) / 2
        else:
            xc = xl
            xr = xu

        pred_c = self.centerbeta[0] + xc * self.centerbeta[1]
        pred_r = self.rangebeta[0] + xr * self.rangebeta[1]

        if cr:
            return pred_c, pred_r
        else:
            return pred_c - pred_r, pred_c + pred_r
