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


def data_generation2(n, a, b, c, d, e, f, g, h, i, j, seed):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(e, f) * rng.chisquare(df=1, size=n)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j) * rng.chisquare(df=1, size=n)

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


def data_generation_outlier(n, a, b, c, d, e, f, g, h, i, j, k, l, m, o, seed=0, alpha=0.1):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n+125)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n+125)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n+125)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr
    # xr = rng.uniform(i, j, size=n)
    n1 = n
    alpha = alpha

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    indices_c = rng.choice(n1, size=int(n1 * alpha), replace=False)
    train[2][indices_c] = train[2][indices_c] + rng.uniform(k, l)
    indices_r = rng.choice(n1, size=int(n1 * alpha), replace=False)
    train[3][indices_r] = train[3][indices_r] + rng.uniform(m, o)

    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]
    # print(beta0, beta1, betastar)

    return train, test, beta0, beta1, betastar


def data_generation_outlier1(n, a, b, c, d, e, f, g, h, i, j, seed=0, alpha=0.1, outlier='Central'):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n+100)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n+100)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n+100)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr
    # xr = rng.uniform(i, j, size=n)
    n1 = n
    alpha = alpha

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    indices = rng.choice(n1, size=int(n1 * alpha), replace=False)

    if outlier == 'Central':
        # train[0][indices] = 1
        train[2][indices] = 5
    elif outlier == 'Radius':
        # train[1][indices] = 1
        train[3][indices] = 5
    else:
        # train[0][indices] = 1
        train[2][indices] = 5
        # train[1][indices] = 1
        train[3][indices] = 5

    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]

    return train, test, beta0, beta1, betastar


def data_generation_outlier2(n, a, b, c, d, e, f, g, h, i, j, seed=0, alpha=0.1, outlier='Central'):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n + 100)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n + 100)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n + 100)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr
    # xr = rng.uniform(i, j, size=n)
    n1 = n
    alpha = alpha

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    indices = rng.choice(n1, size=int(n1 * alpha), replace=False)

    if outlier == 'Central':
        train[0][indices] = 1
        train[2][indices] = 1000
    elif outlier == 'Radius':
        train[1][indices] = 1
        train[3][indices] = 1000
    else:
        train[0][indices] = 1
        train[2][indices] = 1000
        train[1][indices] = 1
        train[3][indices] = 1000

    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]

    return train, test, beta0, beta1, betastar


def data_generation_outlier3(n, a, b, c, d, e, f, g, h, i, j, seed=0, alpha=0.1, outlier='Central'):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n + 100)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n + 100)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n + 100)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr
    # xr = rng.uniform(i, j, size=n)
    n1 = n
    alpha = alpha

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    indices = rng.choice(n1, size=int(n1 * alpha), replace=False)

    if outlier == 'Central':
        train[0][indices] = rng.uniform(8, 10, size=int(n1 * alpha))
        train[2][indices] = rng.uniform(8, 10, size=int(n1 * alpha))
    elif outlier == 'Radius':
        train[1][indices] = rng.uniform(2, 4, size=int(n1 * alpha))
        train[3][indices] = rng.uniform(2, 4, size=int(n1 * alpha))
    else:
        train[0][indices] = rng.uniform(8, 10, size=int(n1 * alpha))
        train[2][indices] = rng.uniform(8, 10, size=int(n1 * alpha))
        train[1][indices] = rng.uniform(2, 4, size=int(n1 * alpha))
        train[3][indices] = rng.uniform(2, 4, size=int(n1 * alpha))

    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]

    return train, test, beta0, beta1, betastar


def data_generation_outlier4(n, a, b, c, d, e, f, g, h, i, j, seed=0, alpha=0.1, outlier='Central'):
    rng = np.random.default_rng(seed)
    xc = rng.uniform(low=a, high=b, size=n + 100)
    beta0 = rng.uniform(c, d)
    beta1 = rng.uniform(c, d)
    epsilon = rng.uniform(low=e, high=f, size=n + 100)
    yc = beta0 + beta1 * xc + epsilon

    betastar = rng.uniform(g, h)
    epsilonr = rng.uniform(i, j, size=n + 100)

    xr = betastar * xc + epsilonr
    yr = betastar * yc + epsilonr

    # yr = betastar * (beta0 + beta1 * (xr - epsilonr) / betastar) + epsilonr
    #  = betastar * beta0 + beta1*xr +(1-beta1)*epsilonr
    # xr = rng.uniform(i, j, size=n)
    n1 = n
    alpha = alpha

    train = [xc[:n1], xr[:n1], yc[:n1], yr[:n1]]
    indices = rng.choice(n1, size=int(n1 * alpha), replace=False)

    if outlier == 'Central':
        # train[0][indices] = train[0][indices] + 1 * np.random.standard_t(df=1, size=int(n1 * alpha))
        train[2][indices] = train[2][indices] + 1 * np.random.standard_t(df=3, size=int(n1 * alpha))
    elif outlier == 'Radius':
        # train[1][indices] = train[1][indices] + 1 * abs(np.random.standard_t(df=1, size=int(n1 * alpha)))
        train[3][indices] = train[3][indices] + 1 * abs(np.random.standard_t(df=3, size=int(n1 * alpha)))
    else:
        # train[0][indices] = train[0][indices] + 1 * np.random.standard_t(df=1, size=int(n1 * alpha))
        train[2][indices] = train[2][indices] + 1 * np.random.standard_t(df=3, size=int(n1 * alpha))
        # train[1][indices] = train[1][indices] + 1 * abs(np.random.standard_t(df=1, size=int(n1 * alpha)))
        train[3][indices] = train[3][indices] + 1 * abs(np.random.standard_t(df=3, size=int(n1 * alpha)))

    test = [xc[n1:], xr[n1:], yc[n1:], yr[n1:]]

    return train, test, beta0, beta1, betastar
