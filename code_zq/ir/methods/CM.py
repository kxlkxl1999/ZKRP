import numpy as np
from ..utils.utils import simplelstsq

# Center method
class CM:
    def __init__(self, xc, yc):
        """
        # Inputs:
        # xc: the center of input interval Xs
        # yc: the center of output interval Ys
        """
        assert xc.shape[0] == yc.shape[0]
        self.xc = xc
        self.yc = yc

    def fit(self):
        self.hatbeta0, self.hatbeta1 = simplelstsq(self.xc, self.yc)
        self.hatbeta = np.array([self.hatbeta0, self.hatbeta1])

        return self.hatbeta

    def predict(self, xl, xu, cr=False):
        if cr:
            left = xl - xu
            right = xl + xu
        else:
            left = xl
            right = xu

        # print(left.shape, self.hatbeta0, self.hatbeta1)

        pred_l = self.hatbeta0 + left * self.hatbeta1
        pred_u = self.hatbeta0 + right * self.hatbeta1

        if cr:
            return (pred_l + pred_u) / 2, (pred_u - pred_l) / 2
        else:
            return pred_l, pred_u
