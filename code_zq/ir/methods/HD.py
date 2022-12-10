import numpy as np
from scipy.optimize import minimize
from ..utils.distances import hausdorff_dist_vec


# Minimum Hausdorff distance based methods
class HD:
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
        def hd_obj(para, xc, xr, yc, yr):
            a_c = para[0]
            b_c = para[1]
            a_r = para[2]
            b_r = para[3]

            yc_hat = a_c + b_c * xc
            yr_hat = a_r + b_r * xr
            return hausdorff_dist_vec(
                np.vstack([yc, yr]).T, np.vstack([yc_hat, yr_hat]).T, cr=True
            ).mean()

        # def hd_obj(para, xc, xr, yc, yr):
        #     a = para[0]
        #     b = para[1]

        #     yc_hat = a + b * xc
        #     yr_hat = a + b * xr

        #     # return np.mean(
        #     #     np.min(np.vstack([abs(yc - yc_hat), abs(yr - yr_hat)]), 0) ** 2
        #     # )
        #     return np.mean(
        #         np.min(
        #             np.vstack(
        #                 [
        #                     (yc - yr - (yc_hat - yr_hat)) ** 2,
        #                     (yc + yr - (yc_hat + yr_hat)) ** 2,
        #                 ]
        #             ),
        #             0,
        #         )
        #     )

        result = minimize(
            hd_obj,
            x0=np.array([1, 1, 1, 1]),
            args=(self.xc, self.xr, self.yc, self.yr),
            method="Nelder-Mead",
        )
        # print(result)
        self.centerbeta = np.array([result["x"][0], result["x"][1]])
        self.rangebeta = np.array([result["x"][2], result["x"][3]])
        # self.hatbeta = np.array([result["x"][0], result["x"][1]])
        # print(self.hatbeta)
        # return self.hatbeta
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

    # def predict(self, xl, xu, cr=False):
    #     if cr:
    #         left = xl - xu
    #         right = xl + xu
    #     else:
    #         left = xl
    #         right = xu

    #     # print(left.shape, self.hatbeta0, self.hatbeta1)

    #     pred_l = self.hatbeta[0] + left * self.hatbeta[1]
    #     pred_u = self.hatbeta[0] + right * self.hatbeta[1]

    #     if cr:
    #         return (pred_l + pred_u) / 2, (pred_u - pred_l) / 2
    #     else:
    #         return pred_l, pred_u
