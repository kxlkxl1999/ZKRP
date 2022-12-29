import numpy as np

from .ir.methods.CM import CM
from .ir.methods.CRM import CRM
from .ir.methods.CCRM import CCRM
from .ir.methods.InterArith import WL2
from .ir.methods.HD import HD
# from .ir.utils.utils import evaluation


# def cr_indep_data_generation(n, a, b, c, d, e, f, g, h, i, j, seed):
#     rng = np.random.default_rng(seed)
#     xc = rng.uniform(low=a, high=b, size=n)
#     beta0 = rng.uniform(c, d)
#     beta1 = rng.uniform(c, d)
#     epsilon = rng.uniform(low=e, high=f, size=n)
#     yc = beta0 + beta1 * xc + epsilon

#     # betastar = rng.uniform(g, h)
#     # epsilonr = rng.uniform(i, j, size=n)

#     xr = rng.uniform(g, h, size=n)
#     yr = rng.uniform(i, j, size=n)
#     # xr = rng.uniform(i, j, size=n)

#     n1 = 250

#     train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
#     test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]
#     # print(beta0, beta1, betastar)

#     return train, test, beta0, beta1


def data_generation(n, a, b, c, d, e, f, g, h, i, j, seed):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr

    # xr = rng.uniform(i, j, size=n)

    n1 = 250

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]
    # print(beta0, beta1, betastar)

    return train, test, beta0, beta1, betastar


def kangxinlai(n, scenario, seed):
    rng = np.random.default_rng(seed)
    xr = abs(rng.random(n)) * 2
    xc = rng.random(n) * 10
    # xx = np.vstack((xx_c - xx_r, xx_c + xx_r)).T
    # xx_r = xx_r.reshape((n, 1))
    # xx_c = xx_c.reshape((n, 1))
    if scenario == 1:
        yr = abs(xr + 0.25 * rng.standard_normal(n))
        yc = 0.5 * xc + 1 + 0.15 * rng.standard_normal(n)
    elif scenario == 2:
        yr = abs(xr + 0.25 * rng.standard_normal(n))
        yc = 0.5 * xc + 1 + 0.5 * rng.standard_normal(n)
    elif scenario == 3:
        yr = abs(xr + 0.25 * rng.standard_normal(n))
        yc = 0.5 * xc + 1 + 2 * rng.standard_normal(n)
    elif scenario == 4:
        yr = abs(xr + 0.5 * rng.standard_normal(n))
        yc = 0.5 * xc + 1 + 0.5 * rng.standard_normal(n)
    elif scenario == 5:
        xc = rng.standard_normal(n)
        epsilonc = rng.standard_normal(n)
        xr = rng.chisquare(1, size=n)
        epsilonr = rng.chisquare(1, size=n)
        yl = 2 * (xc - xr) + (epsilonc - epsilonr)
        yu = 2 * (xc + xr) + (epsilonc + epsilonr)
        yc = (yl + yu) / 2
        yr = (yu - yl) / 2
    elif scenario == 6:
        xc = rng.standard_normal(n)
        epsilonc = 4 * rng.standard_normal(n)
        xr = rng.chisquare(1, size=n)
        epsilonr = rng.chisquare(1, size=n)
        yl = 2 * (xc - xr) + (epsilonc - epsilonr)
        yu = 2 * (xc + xr) + (epsilonc + epsilonr)
        yc = (yl + yu) / 2
        yr = (yu - yl) / 2

    n1 = 1000
    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]

    return train, test


# def comparison(a, b, c, d, e, f, g, h, i, j):
#     print(a, b, c, d, e, f, g, h, i, j)
#     # def comparison(scenario):
#     # print(a, b, c, d, e, f, g, h, i, j)
#     cm_rmsel, cm_rmseu = [], []
#     crm_rmsel, crm_rmseu = [], []
#     ccrm_rmsel, ccrm_rmseu = [], []
#     wl2_rmsel, wl2_rmseu = [], []
#     hd_rmsel, hd_rmseu = [], []
#
#     for iter in range(100):
#         # train, test, _, _ = cr_indep_data_generation(
#         #     375, a, b, c, d, e, f, g, h, i, j, iter
#         # )
#         train, test, beta0, beta1, betastar = data_generation(375, a, b, c, d, e, f, g, h, i, j, iter)
#         # train, test = kangxinlai(2000, scenario, iter)
#         # CM
#         cm = CM(train[0], train[2])
#         hatbeta_cm = cm.fit()
#         hatyc, hatyr = cm.predict(test[0], test[1], cr=True)
#         rmsel, rmseu = evaluation(test[2], test[3], hatyc, hatyr)
#         cm_rmsel.append(rmsel)
#         cm_rmseu.append(rmseu)
#
#         # CRM
#         crm = CRM(train[0], train[1], train[2], train[3])
#         hatbeta_c_crm, hatbeta_r_crm = crm.fit()
#         hatyc, hatyr = crm.predict(test[0], test[1], cr=True)
#         rmsel, rmseu = evaluation(test[2], test[3], hatyc, hatyr)
#         # print("CRM")
#         # print(crm.centerbeta, crm.rangebeta)
#         crm_rmsel.append(rmsel)
#         crm_rmseu.append(rmseu)
#
#         # CCRM
#         ccrm = CCRM(train[0], train[1], train[2], train[3])
#         hatbeta_c_ccrm, hatbeta_r_ccrm = ccrm.fit()
#         hatyc, hatyr = ccrm.predict(test[0], test[1], cr=True)
#         rmsel, rmseu = evaluation(test[2], test[3], hatyc, hatyr)
#         # print("CCRM")
#         # print(ccrm.centerbeta, ccrm.rangebeta)
#         ccrm_rmsel.append(rmsel)
#         ccrm_rmseu.append(rmseu)
#
#         # WL2
#         wl2 = WL2(train[0], train[1], train[2], train[3])
#         wl2.fit(0.5)
#         hatyc, hatyr = wl2.predict(test[0], test[1], cr=True)
#         rmsel, rmseu = evaluation(test[2], test[3], hatyc, hatyr)
#         # print("WL2")
#         # print(wl2.centerbeta, wl2.rangebeta)
#         wl2_rmsel.append(rmsel)
#         wl2_rmseu.append(rmseu)
#
#         # HD
#         hd = HD(train[0], train[1], train[2], train[3])
#         hatbeta_c_hd, hatbeta_r_hd = hd.fit()
#         hatyc, hatyr = hd.predict(test[0], test[1], cr=True)
#         rmsel, rmseu = evaluation(test[2], test[3], hatyc, hatyr)
#         # print("HD")
#         # print(hd.centerbeta, hd.rangebeta)
#         hd_rmsel.append(rmsel)
#         hd_rmseu.append(rmseu)
#
#     print(np.array(cm_rmsel).mean(), np.array(cm_rmseu).mean())
#     print(np.array(crm_rmsel).mean(), np.array(crm_rmseu).mean())
#     print(np.array(ccrm_rmsel).mean(), np.array(ccrm_rmseu).mean())
#     print(np.array(wl2_rmsel).mean(), np.array(wl2_rmseu).mean())
#     print(np.array(hd_rmsel).mean(), np.array(hd_rmseu).mean())


# CRM better than CCRM
# comparison(20, 40, 0, 1, -20, 20, 20, 40, 20, 40)
# comparison(20, 40, 0, 1, -5, 5, 20, 40, 20, 40)
# comparison(20, 40, 0, 1, -20, 20, 1, 5, 1, 5)
# comparison(20, 40, 0, 1, -5, 5, 1, 5, 1, 5)

# CCRM better than CRM
# comparison(20, 40, 0, 1, -20, 20, 20, 40, 1, 5)
# comparison(20, 40, 0, 1, -20, 20, 20, 40, 10, 20)
# comparison(20, 40, 1, 10, 1, 5, 1, 5, 1, 5)
# comparison(20, 40, 0, 1, -5, 5, 1, 5, 10, 20)
# comparison(10, 20, 0, 1, -1, 1, 1, 2, 1, 2)
# for i in range(5, 7):
#     # for i in [3]:
#     print("=============================")
#     comparison(i)
#     print("=============================")
