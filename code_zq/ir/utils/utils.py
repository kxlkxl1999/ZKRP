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


# performance metrics
def evaluation(yc, yr, hatyc, hatyr, criterion="RMSE"):
    """
    performance metric for interval linear regression

    :param yc: numpy-array; center of response interval
    :param yr: numpy-array; range of response interval
    :param hatyc: numpy-array; predicted center of response interval
    :param hatyr: numpy-array; predicted range of response interval
    :param criterion:
    'RMSE': root-mean-square error
    'MHD': mean Hausdorff Distance
    'IOR': interval overlap ratio
    'DC': determination coefficient

    :return:
    'RMSE': (float, float); root-mean-square error of central and range
    'MHD': float; mean Hausdorff Distance of real and estimated data
    'IOR': float;
    'DC': float;
    """

    if criterion == "RMSE":
        rmsel = np.sqrt(pow(yc - yr - hatyc + hatyr, 2).mean())
        rmseu = np.sqrt(pow(yc + yr - hatyc - hatyr, 2).mean())
        # rmsel = np.sqrt(pow(yc - hatyc, 2).mean())
        # rmseu = np.sqrt(pow(yr - hatyr, 2).mean())
        return rmsel, rmseu
    elif criterion == "MHD":
        return (
            hausdorff_dist_vec(
                np.vstack([yc, yr]).T, np.vstack([hatyc, hatyr]).T, cr=True
            ).mean(),
            hausdorff_dist_vec(
                np.vstack([yc, yr]).T, np.vstack([hatyc, hatyr]).T, cr=True
            ).mean(),
        )
    elif criterion == "IOR":
        ior = []
        for i in range(n):
            if (y_l[i] - yh_l[i]) * (y_u[i] - yh_u[i]) >= 0:
                if abs(y_c[i] - yh_c[i]) >= (abs(y_r[i]) + abs(yh_r[i])):
                    ior.append(0)
                else:
                    ior.append(
                        min(abs(y_l[i] - yh_u[i]), abs(y_u[i] - yh_l[i]))
                        / (2 * abs(y_r[i]))
                    )
            else:
                ior.append(min(abs(y_r[i]), abs(yh_r[i])) / abs(y_r[i]))
        return np.mean(ior)

    elif criterion == "DC":
        Rc2 = np.mean((hatyc - yc.mean()) ** 2) / np.mean((yc - yc.mean()) ** 2)
        Rr2 = np.mean((hatyr - yr.mean()) ** 2) / np.mean((yr - yr.mean()) ** 2)
        return Rc2, Rr2

    else:
        raise Exception(
            "Wrong criterion. Criterion can only be chosen from 'RMSE', 'MHD', 'IOR', 'DC'. "
        )


def show(data, path, result=None):
    """
    show the distribution of interval data
    :param data: numpy[n,2]: interval data, each line indicates an interval
    :param path: path of saving fig
    :param result: dict: the frechet_mean of the interval dataset
    :return: none
    """
    n = np.shape(data)[0]
    y = range(n)
    for i in range(n):
        plt.plot(data[i, :], np.array([y[i], y[i]]), alpha=1, color="black")
    if result:
        plt.plot(
            result["interval"],
            np.array([round(n / 2), round(n / 2)]),
            alpha=1,
            color="red",
            linewidth=2,
        )
    plt.title("plot of interval")  # 折线图标题
    plt.xlabel("Range")  # x轴标题
    plt.ylabel("Number")  # y轴标题
    plt.grid(ls="-.")  # 绘制背景线
    plt.tight_layout()
    plt.savefig(path)
    plt.show()


def show2(x, y, path, Cov=None):
    """
    show the distribution od two interval datasets
    :param x: numpy[n,2]: interval data, each line indicates an interval
    :param y: numpy[n,2]: interval data, each line indicates an interval
    :param path: path of saving fig
    :param Cov: double: frechet_cov of x and y
    :return: none
    """
    n = np.shape(x)[0]
    y1 = np.array(range(n))
    y2 = y1 + 0.1
    a = None
    b = None
    for i in range(n):
        (a,) = plt.plot(
            x[i, :], np.array([y1[i], y1[i]]), alpha=1, color="red", label="X"
        )
        (b,) = plt.plot(
            y[i, :], np.array([y2[i], y2[i]]), alpha=1, color="blue", label="Y"
        )
    if Cov:
        plt.title(
            "plot of " + r"$X_1$" + " and " + r"$X_2$" + " with Cov %.4f" % Cov
        )  # 折线图标题
    else:
        plt.title("plot of " + r"$X_1$" + " and " + r"$X_2$")
    plt.xlabel("Range")  # x轴标题
    plt.ylabel("Number")  # y轴标题
    plt.grid(ls="-.")  # 绘制背景线
    plt.legend(loc="best", handles=[a, b])
    plt.tight_layout()
    plt.savefig(path)
    plt.show()


def show3(x, y, samples=25, path=None, real=None, Cor=None):
    np.random.seed(0)
    n = np.shape(x)[0]
    plt.figure(dpi=200, figsize=(5, 5))
    for i in np.random.choice(n, samples, replace=False):
        plt.fill_between(
            x[i, :],
            np.array([y[i, 1], y[i, 1]]),
            np.array([y[i, 0], y[i, 0]]),
            alpha=0.25,
            linewidth=0,
            color="steelblue",
        )
        if real is not None:
            plt.scatter(real[i, 0], real[i, 1], alpha=1, color="black", s=2)

    if Cor is not None:
        plt.title(
            "plot of " + r"$X_1$" + " and " + r"$X_2$" + " with Cor %.4f" % Cor
        )  # 折线图标题
    else:
        plt.title("plot of " + r"$X_1$" + " and " + r"$X_2$")
    plt.xlabel(r"$X_1$")  # x轴标题
    plt.ylabel(r"$X_2$")  # y轴标题
    plt.grid(ls="-.")  # 绘制背景线
    # plt.legend(loc='best', handles=[a, b])
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()
