import numpy as np
from .distances import *
from cvxopt import matrix
from cvxopt import solvers
from scipy.stats import pearsonr
from sklearn.linear_model import LinearRegression


def simplelstsq(x, y):
    hatbeta1 = pearsonr(x, y)[0] * np.sqrt(np.cov(y) / np.cov(x))
    hatbeta0 = y.mean() - x.mean() * hatbeta1
    return np.array([hatbeta0, hatbeta1])
    # print("simple linear regression sanity check")
    # reg = LinearRegression().fit(x.reshape((-1, 1)), y.reshape((-1, 1)))
    # print(np.array([hatbeta0, hatbeta1]))
    # print(np.array([reg.intercept_[0], reg.coef_[0, 0]]))


def medianregression(x, y, positivebta=False):
    """
    # CVXOPT linear programming standard form
    # min c^Tx
    # subject to Gx <=h, Ax=b

    # median regression as linear programming
    # min \sum_{i}(ei+ + ei-)
    # subject to ei+, ei- >=0
    # ei+ - ei- = yi - (a+b*xi)
    """
    n = x.shape[0]
    # equality constraints
    A = np.vstack([np.ones_like(x), x, np.eye(n), -np.eye(n)]).T
    b = y
    # the first two elements corresponds to the regression coefficient
    # the next n elements corresponds to the positive part of residual
    # the last n elements corresponds to the negative part of residual
    c = np.concatenate([np.zeros(2), np.ones(2 * n)])

    # inequality constraints
    if positivebta:
        G = -np.eye(2 + 2 * n)
        h = np.zeros(2 + 2 * n)
    else:
        G = np.concatenate(
            [
                np.zeros((2 * n, 2)),
                np.concatenate([np.diag(-np.ones(n)), np.zeros((n, n))], 0),
                np.concatenate([np.zeros((n, n)), np.diag(-np.ones(n))], 0),
            ],
            1,
        )
        h = np.zeros(2 * n)

    # solving the model
    sol = solvers.lp(
        matrix(c),
        matrix(G),
        matrix(h),
        matrix(A),
        matrix(b),
        options={"show_progress": False},
    )
    z = sol["x"]

    beta = np.array(z[:2])
    return beta


def HF_medianregression(x, y, positivebta=False):
    """
    # CVXOPT linear programming standard form
    # min c^Tx
    # subject to Gx <=h, Ax=b

    # median regression as linear programming
    # min \sum_{i}(eci+ + eci- + eri+ + eri-)
    # subject to eci+, eci-, eri+, eri- >=0
    # eci+ - eci- = yci - (ac+bc*xci)
    # eri+ - eri- = yri - (ar+br*xri)
    # return beta_c, beta_r
    """
    n = x.shape[0]
    x_l = x[:, 0]
    x_u = x[:, 1]
    x_c = x_u / 2 + x_l / 2
    x_r = x_u / 2 - x_l / 2
    y_l = y[:, 0]
    y_u = y[:, 1]
    y_c = y_u / 2 + y_l / 2
    y_r = y_u / 2 - y_l / 2
    # equality constraints
    Ac = np.vstack([np.ones_like(x_c), x_c, np.eye(n), -np.eye(n)]).T
    Ar = np.vstack([np.ones_like(x_r), x_r, np.eye(n), -np.eye(n)]).T
    A = np.concatenate([
        np.concatenate([Ac, np.zeros_like(Ac)], 0),
        np.concatenate([np.zeros_like(Ar), Ar], 0)
    ],1)
    b = np.concatenate([y_c, y_r])
    # the first two elements corresponds to the regression coefficient
    # the next n elements corresponds to the positive part of residual
    # the last n elements corresponds to the negative part of residual
    c = np.concatenate([np.zeros(2), np.ones(2 * n), np.zeros(2), np.ones(2 * n)])

    # inequality constraints
    if positivebta:
        G = np.diag(np.concatenate([np.zeros(2), -np.ones(2 + 4 * n)]))
        h = np.zeros(4 + 4 * n)
    else:
        G0 = np.concatenate(
            [
                np.zeros((2 * n, 2)),
                np.concatenate([np.diag(-np.ones(n)), np.zeros((n, n))], 0),
                np.concatenate([np.zeros((n, n)), np.diag(-np.ones(n))], 0)
            ],
            1
        )
        G = np.concatenate([
            np.concatenate([G0, np.zeros_like(G0)], 0),
            np.concatenate([np.zeros_like(G0), G0], 0)
        ], 1)
        h = np.zeros(4 * n)

    # solving the model
    sol = solvers.lp(
        matrix(c),
        matrix(G),
        matrix(h),
        matrix(A),
        matrix(b),
        options={"show_progress": False},
    )
    z = sol["x"]

    beta_c = np.array(z[:2])
    beta_r = np.array(z[(2*n+2):(2*n+4)])
    return beta_c, beta_r


def HF_quadraticregression(x, y, positivebta=False):
    """
    # CVXOPT Quadratic programming standard form
    # min 1/2 * x^TPx + q^Tx
    # subject to Gx <=h, Ax=b

    # median regression as linear programming
    # min \sum_{i}(eci+ + eci- + eri+ + eri-)^2
    # subject to eci+, eci-, eri+, eri- >=0
    # eci+ - eci- = yci - (ac+bc*xci)
    # eri+ - eri- = yri - (ar+br*xri)
    # return np.array([ac, bc, ar, br])
    """
    n = x.shape[0]
    x_l = x[:, 0]
    x_u = x[:, 1]
    x_c = x_u / 2 + x_l / 2
    x_r = x_u / 2 - x_l / 2
    y_l = y[:, 0]
    y_u = y[:, 1]
    y_c = y_u / 2 + y_l / 2
    y_r = y_u / 2 - y_l / 2
    # equality constraints
    Ac = np.vstack([np.ones_like(x_c), x_c, np.eye(n), -np.eye(n)]).T
    Ar = np.vstack([np.ones_like(x_r), x_r, np.eye(n), -np.eye(n)]).T
    A = np.concatenate([
        np.concatenate([Ac, np.zeros_like(Ac)], 0),
        np.concatenate([np.zeros_like(Ar), Ar], 0)
    ], 1)
    b = np.concatenate([y_c, y_r])
    # the first two elements corresponds to the central regression coefficient
    # the next n elements corresponds to the positive part of central residual
    # the next n elements corresponds to the negative part of central residual
    # the next two elements corresponds to the radius regression coefficient
    # the next n elements corresponds to the positive part of radius residual
    # the next n elements corresponds to the negative part of radius residual
    P0 = np.concatenate([
        np.concatenate([np.diag(np.zeros(2)), np.zeros((2 * n, 2))]),
        np.concatenate([np.zeros((2, n)), np.eye(n), np.eye(n)]),
        np.concatenate([np.zeros((2, n)), np.eye(n), np.eye(n)])
    ], 1)
    P = np.concatenate([
        np.concatenate([P0, P0]),
        np.concatenate([P0, P0])
    ], 1)
    q = np.zeros(4 * n + 4)

    # inequality constraints
    if positivebta:
        G = np.diag(np.concatenate([np.zeros(2), -np.ones(2 + 4 * n)]))
        h = np.zeros(4 + 4 * n)
    else:
        G0 = np.concatenate(
            [
                np.zeros((2 * n, 2)),
                np.concatenate([np.diag(-np.ones(n)), np.zeros((n, n))], 0),
                np.concatenate([np.zeros((n, n)), np.diag(-np.ones(n))], 0)
            ],
            1
        )
        G = np.concatenate([
            np.concatenate([G0, np.zeros_like(G0)], 0),
            np.concatenate([np.zeros_like(G0), G0], 0)
        ], 1)
        h = np.zeros(4 * n)

    # solving the model
    sol = solvers.qp(
        matrix(P),
        matrix(q),
        matrix(G),
        matrix(h),
        matrix(A),
        matrix(b),
        options={"show_progress": False},
    )
    z = sol["x"]
    beta_c = np.array(z[:2])
    beta_r = np.array(z[(2 * n + 2):(2 * n + 4)])
    return beta_c, beta_r


def constrained_simplelstsq(x, y):
    """
    The objective function is qudratic in beta0 and beta1
    \beta_0^2 + \beta_1^2 \bar{x^2} - 2\beta_0\bar{y}
    -2\beta_1\bar{xy} + 2\beta_0\beta_1\bar{x}
    """
    P = matrix(np.array([[1, x.mean()], [x.mean(), np.mean(x**2)]]), tc="d")
    q = matrix(np.array([-y.mean(), -np.mean(x * y)]), tc="d")
    G = matrix(np.array([[-1, 0], [0, -1]]), tc="d")
    h = matrix(np.array([0, 0]), tc="d")
    sol = solvers.qp(P, q, G, h, options={"show_progress": False})

    return np.array(sol["x"])


if __name__ == "__main__":
    n = 100
    x = np.random.randn(n)
    y = 2 + 1.5 * x + 0.01 * np.random.randn(n)
    print(medianregression(x, y))