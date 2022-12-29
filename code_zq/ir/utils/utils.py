import numpy as np
from distances import *
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
    x = sol["x"]

    beta = x[:2]
    return np.array(beta)


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