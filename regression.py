from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from zkrp import *
from scipy.optimize import minimize
from collections import OrderedDict
from code_zq.ir.utils import utils
from time import process_time


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
    loop = 0

    while True and loop < 1000:
        loop += 1
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
        loop2 = 0

        while True and loop2 < 1000:
            loop2 += 1
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


def predict(x, beta):
    n = x.shape[0]
    vec_one = np.ones(n)
    x_1 = np.insert(x, 0, vec_one, axis=1)
    return x_1.dot(beta)


def eval(y1, y2, yh1, yh2, datatype='CR', method='RMSE', detail=False):
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
        y_l = y_c - y_r
        y_u = y_c + y_r

        yh_c = yh1.reshape(n)
        yh_r = yh2.reshape(n)
        yh_l = yh_c - yh_r
        yh_u = yh_c + yh_r
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
        rmse_l2 = pow(y_l - yh_l, 2).mean()
        rmse_u2 = pow(y_u - yh_u, 2).mean()
        rmse_r2 = pow(y_r - yh_r, 2).mean()
        rmse_c2 = pow(y_c - yh_c, 2).mean()
        if detail:
            return np.sqrt(rmse_l2), np.sqrt(rmse_u2), np.sqrt(rmse_c2), np.sqrt(rmse_r2)
        else:
            return np.sqrt(rmse_c2)
    elif method == 'NHD':
        return np.mean([hausdorff_distance(np.array([y_l[i], y_u[i]]), np.array([yh_l[i], yh_u[i]])) for i in range(n)])
    elif method == 'IOR':
        ior = []
        for i in range(n):
            if (y_l[i] - yh_l[i]) * (y_u[i] - yh_u[i]) > 0:
                if abs(y_c[i] - yh_c[i]) >= (abs(y_r[i])+abs(yh_r[i])):
                    ior.append(0)
                else:
                    ior.append(min(abs(y_l[i] - yh_u[i]), abs(y_u[i] - yh_l[i])) / (2 * abs(y_r[i])))
            else:
                if y_r[i] == 0:
                    return 1
                else:
                    ior.append(min(abs(y_r[i]), abs(yh_r[i])) / abs(y_r[i]))
        return np.mean(ior)

    elif method == 'DC':
        # return 1 - np.sum((y_r - yh_r) ** 2) / np.sum((y_r - np.mean(y_r)) ** 2)
        if detail:
            return np.corrcoef(y_c, yh_c)[0][1] ** 2, np.corrcoef(abs(y_r), abs(yh_r))[0][1] ** 2
        else:
            return np.corrcoef(abs(y_r), abs(yh_r))[0][1] ** 2 /2 + np.corrcoef(y_c, yh_c)[0][1] ** 2 /2
        # return np.corrcoef(y_c, yh_c)[0][1] ** 2
    # elif method == 'DC' and datatype == 'LU':
    #     return (1 - np.sum((y_l - yh_l) ** 2) / np.sum((y_l - np.mean(y_l)) ** 2))/2 + (1 - np.sum((y_u - yh_u) ** 2) / np.sum((y_u - np.mean(y_u)) ** 2))/2
    else:
        raise Exception('Wrong method. method can only be chosen from \'RMSE\', \'LU\', \'NHD\', \'DC\'. ')


def reg_obj1(para, x, y):
    assert x.shape[0] == y.shape[0]
    n = x.shape[0]
    a_c = para[0]
    b_c = para[1]
    a_r = para[2]
    b_r = para[3]
    dist1 = 0
    for i in range(n):
        y_c_hat = a_c * (x[i, 1] / 2 + x[i, 0] / 2) + b_c
        y_r_hat = a_r * (x[i, 1] / 2 - x[i, 0] / 2) + b_r
        # dist2 += max((y_c_hat - y_r_hat - y[i, 0])**2, (y_c_hat + y_r_hat - y[i, 1])**2)
        dist1 += max(abs(y_c_hat - y_r_hat - y[i, 0]), abs(y_c_hat + y_r_hat - y[i, 1]))
    return dist1 / n


def reg_obj2(para, x, y):
    assert x.shape[0] == y.shape[0]
    n = x.shape[0]
    a_c = para[0]
    b_c = para[1]
    a_r = para[2]
    b_r = para[3]
    dist2 = 0
    for i in range(n):
        y_c_hat = a_c * (x[i, 1] / 2 + x[i, 0] / 2) + b_c
        y_r_hat = a_r * (x[i, 1] / 2 - x[i, 0] / 2) + b_r
        dist2 += max((y_c_hat - y_r_hat - y[i, 0])**2, (y_c_hat + y_r_hat - y[i, 1])**2)
        # dist2 += max(abs(y_c_hat - y_r_hat - y[i, 0]), abs(y_c_hat + y_r_hat - y[i, 1]))
    return dist2 / n


def data_expansion(x1, x2):
    """
    :param x2: numpy (n,)
    :param x1: numpy (n,)
    :return: list of numpy array
    """
    assert len(x1) == len(x2)
    n = len(x1)
    res_list = []
    deta_list = []
    for i in range(n):
        deta_list.append(abs(x1[i] - x2[i]) / 2)

    res_list.append(x1)
    res_list.append(x2)
    for i in range(n):
        deta = np.zeros(n)
        deta[i] = deta_list[i]
        res_list.append(x1 + deta)
        res_list.append(x2 + deta)
        res_list.append(x1 - deta)
        res_list.append(x2 - deta)

    return res_list


def Medreg_Method(x, y):
    n = x.shape[0]
    x_l = x[:, 0]
    x_u = x[:, 1]
    x_c = x_u / 2 + x_l / 2
    x_r = x_u / 2 - x_l / 2
    y_l = y[:, 0]
    y_u = y[:, 1]
    y_c = y_u / 2 + y_l / 2
    y_r = y_u / 2 - y_l / 2
    beta_c = utils.medianregression(x_c, y_c)
    beta_r = utils.medianregression(x_r, y_r, positivebta=True)
    return beta_c[1][0], beta_c[0][0], beta_r[1][0], beta_r[0][0]


def HF_Med_Method(x, y):
    beta_c, beta_r = utils.HF_medianregression(x, y, positivebta=True)
    return beta_c[1][0], beta_c[0][0], beta_r[1][0], beta_r[0][0]


def HF_Qd_Method(x, y):
    beta_c, beta_r = utils.HF_quadraticregression(x, y, positivebta=True)
    return beta_c[1][0], beta_c[0][0], beta_r[1][0], beta_r[0][0]


def HF_Method1(x, y, method='Nelder-Mead'):
    n = x.shape[0]
    x_l = x[:, 0].reshape((n, 1))
    x_u = x[:, 1].reshape((n, 1))
    x_c = x_u / 2 + x_l / 2
    x_r = x_u / 2 - x_l / 2
    y_l = y[:, 0].reshape((n, 1))
    y_u = y[:, 1].reshape((n, 1))
    y_c = y_u / 2 + y_l / 2
    y_r = y_u / 2 - y_l / 2
    # beta_0 = CM_Method(x_c, y_c)
    # print('=')
    beta_c1, beta_r1 = CRM_Method(x_c, y_c, x_r, y_r)
    # print('==')
    beta_c2, beta_r2 = CCRM_Method(x_c, y_c, x_r, y_r)
    # print('===')
    # x0_list = data_expansion(np.array([beta_c1[1], beta_c1[0], beta_r1[1], beta_r1[0]]),
    #                          np.array([beta_c2[1], beta_c2[0], beta_r2[1], beta_r2[0]]))
    x0_list = [np.array([beta_c1[1], beta_c1[0], beta_r1[1], beta_r1[0]]),
               np.array([beta_c2[1], beta_c2[0], beta_r2[1], beta_r2[0]])]

    result = [minimize(reg_obj1, x0=i, args=(x, y), method=method) for i in x0_list]
    num_min = np.argmin([i['fun'] for i in result])

    return result[num_min]['x'][0], result[num_min]['x'][1], result[num_min]['x'][2], result[num_min]['x'][3]


def HF_Method2(x, y, method='Nelder-Mead'):
    n = x.shape[0]
    x_l = x[:, 0].reshape((n, 1))
    x_u = x[:, 1].reshape((n, 1))
    x_c = x_u / 2 + x_l / 2
    x_r = x_u / 2 - x_l / 2
    y_l = y[:, 0].reshape((n, 1))
    y_u = y[:, 1].reshape((n, 1))
    y_c = y_u / 2 + y_l / 2
    y_r = y_u / 2 - y_l / 2
    beta_0 = CM_Method(x_c, y_c)
    beta_c1, beta_r1 = CRM_Method(x_c, y_c, x_r, y_r)
    beta_c2, beta_r2 = CCRM_Method(x_c, y_c, x_r, y_r)
    # x0_list = data_expansion(np.array([beta_c1[1], beta_c1[0], beta_r1[1], beta_r1[0]]),
    #                          np.array([beta_c2[1], beta_c2[0], beta_r2[1], beta_r2[0]]))
    x0_list = [np.array([beta_c1[1], beta_c1[0], beta_r1[1], beta_r1[0]]),
               np.array([beta_c2[1], beta_c2[0], beta_r2[1], beta_r2[0]])]


    result = [minimize(reg_obj2, x0=i, args=(x, y), method=method) for i in x0_list]
    num_min = np.argmin([i['fun'] for i in result])

    return result[num_min]['x'][0], result[num_min]['x'][1], result[num_min]['x'][2], result[num_min]['x'][3]


def show4(x, y, yhat, samples=25, path=None, real=None, Cor=None):
    np.random.seed(0)
    n = np.shape(x)[0]
    plt.figure(dpi=200, figsize=(5, 5))
    for i in np.random.choice(n,samples,replace=False):
        plt.fill_between(x[i, :], np.array([y[i, 1], y[i, 1]]), np.array([y[i, 0], y[i, 0]]), alpha=0.15, linewidth=0,
                         color='darkorange', label=r'$y$')
        plt.fill_between(x[i, :], np.array([yhat[i, 1], yhat[i, 1]]), np.array([yhat[i, 0], yhat[i, 0]]), alpha=0.3, linewidth=0,
                         color='steelblue', label=r'$\hat{y}$')
        if real is not None:
            plt.scatter(real[i, 0], real[i, 1], alpha=1, color="black", s=2)

    if Cor is not None:
        plt.title('plot of ' + r'$X$'+' and ' + r'$Y$' + ', ' + r'$\hat{y}$' + ' with Cor %.4f' % Cor)  # 折线图标题
    else:
        plt.title('plot of ' + r'$X$'+' and ' + r'$Y$' + ', ' + r'$\hat{Y}$')
    plt.xlabel(r'$X$')  # x轴标题
    plt.ylabel(r'$Y$')  # y轴标题
    plt.grid(ls='-.')  # 绘制背景线

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc='best')

    # plt.legend(loc='best')
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()


