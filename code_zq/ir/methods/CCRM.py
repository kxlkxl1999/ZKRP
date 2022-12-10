import numpy as np
from ..utils.utils import simplelstsq, constrained_simplelstsq


# Constrained center and range method
class CCRM:
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

    def fit(self):
        # fit two separate models
        self.centerbeta = simplelstsq(self.xc, self.yc)
        # for the range, we have constraint so that the
        # regression coefficients are nonnegative
        self.rangebeta = constrained_simplelstsq(self.xr, self.yr)
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
