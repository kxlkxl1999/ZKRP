from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from zkrp import *


# Center method
def CM_Method(input_x_center, input_y_center):
    n = input_x_center.shape[0]
    vec_one = np.ones(n)
    X_c = np.insert(input_x_center, 0, vec_one, axis=1)
    X_cT = X_c.transpose()
    A = np.dot(X_cT, X_c)
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, X_cT)
    hatBeta_c = np.dot(B, input_y_center)

    return hatBeta_c


# Center and range method
def CRM_Method(input_x_center, input_y_center, input_x_width, input_y_width):
    n = input_x_center.shape[0]
    vec_one = np.ones(n)

    X_c = np.insert(input_x_center, 0, vec_one, axis=1)  # 第一列插入1
    X_cT = X_c.transpose()
    A = np.dot(X_cT, X_c)
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, X_cT)
    hatBeta_c = np.dot(B, input_y_center)

    X_r = np.insert(input_x_width, 0, vec_one, axis=1)
    X_rT = X_r.transpose()
    A2 = np.dot(X_rT, X_r)
    A2_inv = np.linalg.inv(A2)
    B2 = np.dot(A2_inv, X_rT)
    hatBeta_r = np.dot(B2, input_y_width)

    return hatBeta_c, hatBeta_r


# Constrained center and range method
def CCRM_Method(input_x_center, input_y_center, input_x_width, input_y_width):
    n, p = input_x_center.shape
    vec_one = np.ones(n)

    X_c = np.insert(input_x_center, 0, vec_one, axis=1)
    X_cT = X_c.transpose()
    A = np.dot(X_cT, X_c)
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, X_cT)
    hatBeta_c = np.dot(B, input_y_center)

    # Initialization for beta_r
    hatBeta_r = np.zeros(p + 1).reshape(p + 1, 1)
    P = np.array([], dtype=np.int64)
    Z = np.arange(p + 1)
    X_r = np.insert(input_x_width, 0, vec_one, axis=1)
    X_rT = X_r.transpose()

    while True:
        w = np.dot(X_rT, input_y_width - np.dot(X_r, hatBeta_r))
        neg_Z = np.array([], dtype=np.int64)
        for j in Z:
            if w[j] <= 0:
                neg_Z = np.append(neg_Z, j)

        if Z.size == 0 or Z.size == neg_Z.size:
            break

        t = Z[np.argmax(w[Z])]
        Z = np.delete(Z, np.argwhere(Z == t))
        P = np.append(P, t)

        while True:
            X_p = X_r.copy()
            for j in Z:
                X_p[:, j] = 0

            z = np.zeros(p + 1).reshape(p + 1, 1)
            X_pT = X_p.transpose()
            C = np.dot(X_pT, X_p)
            D = np.dot(X_pT, input_y_width)
            C2 = C[P][:, P]
            D2 = D[P]
            C2_inv = np.linalg.inv(C2)
            z2 = np.dot(C2_inv, D2)
            z[P] = z2

            pos_P = np.array([], dtype=np.int64)
            for j in P:
                if z[j] > 0:
                    pos_P = np.append(pos_P, j)

            if P.size == pos_P.size:
                hatBeta_r = z
                break

            neg_P = np.setdiff1d(P, pos_P)
            k = neg_P[0]
            alpha = hatBeta_r[k] / (hatBeta_r[k] - z[k])
            for j in neg_P:
                if hatBeta_r[j] / (hatBeta_r[j] - z[j]) < alpha:
                    alpha = hatBeta_r[j] / (hatBeta_r[j] - z[j])
                    k = j

            hatBeta_r += alpha * (z - hatBeta_r)

            for j in P:
                if hatBeta_r[j] == 0:
                    P = np.delete(P, np.argwhere(P == j))
                    Z = np.append(Z, j)

    return hatBeta_c, hatBeta_r


# Parametrized method
def PM_Method(input_x_center, input_y_center, input_x_width, input_y_width):
    input_x_lower = input_x_center - input_x_width / 2
    input_x_upper = input_x_center + input_x_width / 2
    input_y_lower = input_y_center - input_y_width / 2
    input_y_upper = input_y_center + input_y_width / 2

    n, p = input_x_center.shape
    X_lu = np.ones(n).reshape(n, 1)
    for j in range(p):
        X_lu = np.insert(X_lu, X_lu.shape[1], input_x_lower[:, j], axis=1)
        X_lu = np.insert(X_lu, X_lu.shape[1], input_x_upper[:, j], axis=1)

    X_luT = X_lu.transpose()
    A = np.dot(X_luT, X_lu)
    A_inv = np.linalg.inv(A)
    B = np.dot(A_inv, X_luT)
    hatBeta_l = np.dot(B, input_y_lower)
    hatBeta_u = np.dot(B, input_y_upper)

    return hatBeta_l, hatBeta_u


# Mean magnitude of error relative to the estimate index
def MMER(input_y_center, estimate_y_center, input_y_width, estimate_y_width):
    n = input_y_center.shape[0]

    input_y_center = input_y_center.reshape(n)
    estimate_y_center = estimate_y_center.reshape(n)
    input_y_width = input_y_width.reshape(n)
    estimate_y_width = estimate_y_width.reshape(n)

    input_y_lower = input_y_center - input_y_width / 2
    input_y_upper = input_y_center + input_y_width / 2
    estimate_y_lower = estimate_y_center - estimate_y_width / 2
    estimate_y_upper = estimate_y_center + estimate_y_width / 2

    ratio1 = abs((input_y_lower - estimate_y_lower) / estimate_y_lower).sum()
    ratio2 = abs((input_y_upper - estimate_y_upper) / estimate_y_upper).sum()

    return (ratio1 + ratio2) / (2 * n)


# The lower bound mean-squared error and upper bound mean-squared error
def RMSE(input_y_center, estimate_y_center, input_y_width, estimate_y_width):
    n = input_y_center.shape[0]

    input_y_center = input_y_center.reshape(n)
    estimate_y_center = estimate_y_center.reshape(n)
    input_y_width = input_y_width.reshape(n)
    estimate_y_width = estimate_y_width.reshape(n)

    input_y_lower = input_y_center - input_y_width / 2
    input_y_upper = input_y_center + input_y_width / 2
    estimate_y_lower = estimate_y_center - estimate_y_width / 2
    estimate_y_upper = estimate_y_center + estimate_y_width / 2

    rmse_l2 = pow(input_y_lower - estimate_y_lower, 2).mean()
    rmse_u2 = pow(input_y_upper - estimate_y_upper, 2).mean()

    return np.sqrt(rmse_l2), np.sqrt(rmse_u2)


def predict(x, beta):
    n = x.shape[0]
    vec_one = np.ones(n)
    x_1 = np.insert(x, 0, vec_one, axis=1)
    return x_1.dot(beta)


def eval(y1, y2, yh1, yh2, datatype='CR', method='RMSE'):
    """
    evaluation of teh regression effect.
    :param y1: numpy-array; real data
    :param y2: numpy-array; real data
    :param yh1: numpy-array; estimated data
    :param yh2: numpy-array; estimated data
    :param datatype: string;
    'CR': central and range with respect to y1 and y2.
    'LU': lower and upper with respect to y1 and y2
    :param method:
    'RMSE': root-mean-square error
    'MHD': mean Hausdorff Distance
    'IOR': interval overlap ratio
    'DC': determination coefficient

    :return:
    'RMSE': (float, float); root-mean-square error of central and range
    'NHD': float; mean Hausdorff Distance of real and estimated data
    'IOR': float;
    'DC': float;
    """
    n = y1.shape[0]
    if datatype == 'CR':
        y_c = y1.reshape(n)
        y_r = y2.reshape(n)
        y_l = y_c - y_r/2
        y_u = y_c + y_r/2

        yh_c = yh1.reshape(n)
        yh_r = yh2.reshape(n)
        yh_l = yh_c - yh_r / 2
        yh_u = yh_c + yh_r / 2
    elif datatype == 'LU':
        y_l = y1.reshape(n)
        y_u = y2.reshape(n)
        y_c = (y_l + y_u) / 2
        y_r = (y_u - y_l) / 2

        yh_l = yh1.reshape(n)
        yh_u = yh2.reshape(n)
        yh_c = (yh_l + yh_u) / 2
        yh_r = (yh_u - yh_l) / 2
    else:
        raise Exception('Wrong datatype. datatype can only be chosen from \'CR\' and \'LU\'. ')

    if method == 'RMSE':
        rmse_c2 = pow(y_r - yh_r, 2).mean()
        # rmse_u2 = pow(y_u - yh_u, 2).mean()
        return np.sqrt(rmse_c2)
    elif method == 'NHD':
        return np.mean([hausdorff_distance(np.array([y_l[i], y_u[i]]), np.array([yh_l[i], yh_u[i]])) for i in range(n)])
    elif method == 'IOR':
        return np.mean([max(0, y_r[i] + yh_r[i] - abs(y_c[i] - yh_c[i]))/(2*abs(y_r[i])) for i in range(n)])
    elif method == 'DC':
        # return 1 - np.sum((y_r - yh_r) ** 2) / np.sum((y_r - np.mean(y_r)) ** 2)
        return np.corrcoef(y_r, yh_r)[0][1] ** 2 /2 + np.corrcoef(y_c, yh_c)[0][1] ** 2 /2
        # return np.corrcoef(y_c, yh_c)[0][1] ** 2
    # elif method == 'DC' and datatype == 'LU':
    #     return 1 - np.sum((y_l - yh_l) ** 2) / np.sum((y_l - np.mean(y_l)) ** 2), 1 - np.sum((y_u - yh_u) ** 2) / np.sum((y_u - np.mean(y_u)) ** 2)
    else:
        raise Exception('Wrong method. method can only be chosen from \'RMSE\', \'LU\', \'NHD\', \'DC\'. ')
