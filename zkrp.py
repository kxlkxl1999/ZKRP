# -*- coding: utf-8 -*-
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


# from scipy.stats import wasserstein_distance


def hausdorff_distance(x, y, sgn=False):
    """
    分别匹配x,y的端点就行四次比较，取最小的里的最大的

    :param x: array-like(2): interval data
    :param y: array-like(2): interval data
    :param sgn: bool: 返回的值是否带符号，x的点在y的点左边为正
    :return:
    """

    a = (1 if (x[0] - y[0]) * (x[0] - y[1]) > 0 else 0) * \
        (abs(x[0] - y[0]) if abs(x[0] - y[0]) < abs(x[0] - y[1]) else abs(x[0] - y[1]))
    b = (1 if (x[1] - y[0]) * (x[1] - y[1]) > 0 else 0) * \
        (abs(x[1] - y[0]) if abs(x[1] - y[0]) < abs(x[1] - y[1]) else abs(x[1] - y[1]))
    c = (1 if (x[0] - y[0]) * (x[1] - y[0]) > 0 else 0) * \
        (abs(x[0] - y[0]) if abs(x[0] - y[0]) < abs(x[1] - y[0]) else abs(x[1] - y[0]))
    d = (1 if (x[0] - y[1]) * (x[1] - y[1]) > 0 else 0) * \
        (abs(x[0] - y[1]) if abs(x[0] - y[1]) < abs(x[1] - y[1]) else abs(x[1] - y[1]))
    m = max(a, b, c, d)
    if sgn:
        return m if sum(x) < sum(y) else -1 * m
    else:
        return m


def hausdorff_distance2(x, y):
    abs_a = 0
    for w in np.hstack((np.linspace(start=x[0], stop=x[1], num=1000), np.linspace(start=y[0], stop=y[1], num=1000))):
        a0 = abs((dist(w, x) - dist(w, y)))
        if abs(a0) > abs_a:
            abs_a = abs(a0)
    return abs_a


def dist(x, y):
    """
    计算一个点到一个区间的距离
    :param x: double: 点的取值
    :param y: numpy[2]: 区间的左右端点
    :return: double:
    """
    return 0 if (y[0] < x < y[1]) else min(abs(x - y[0]), abs(x - y[1]))


def middle_distance(x, y):
    return abs((x[0] + x[1]) / 2 - (y[0] + y[1]) / 2)


def wasserstein_distance2(u_values, v_values):
    """
    wasserstein distance between intervals with squared loss
    :param u_values: numpy(n): empirical data in the interval
    :param v_values: numpy(n): empirical data in the interval
    :return: double
    """
    u_sorter = np.argsort(u_values)
    v_sorter = np.argsort(v_values)

    all_values = np.concatenate((u_values, v_values))
    all_values.sort(kind='mergesort')

    # Compute the differences between pairs of successive values of u and v.
    deltas = np.diff(all_values)

    # Get the respective positions of the values of u and v among the values of
    # both distributions.
    u_cdf_indices = u_values[u_sorter].searchsorted(all_values[:-1], 'right')
    v_cdf_indices = v_values[v_sorter].searchsorted(all_values[:-1], 'right')

    # Calculate the CDFs of u and v using their weights, if specified.

    u_cdf = u_cdf_indices / u_values.size
    v_cdf = v_cdf_indices / v_values.size
    # Compute the value of the integral based on the CDFs.
    # If p = 1 or p = 2, we avoid using np.power, which introduces an overhead
    # of about 15%.
    return np.sum(np.multiply(pow((u_cdf - v_cdf), 2), deltas))


def arithmetic_based_distance(x, y, theta=1):
    return np.sqrt(pow(np.mean(x)-np.mean(y), 2) + theta * pow(x[1]/2-x[0]/2-y[1]/2+y[0]/2, 2))


def fdsum(midpont, data, inter_l, method='hausdorff'):
    inter_len = inter_l
    interval = [midpont - inter_len / 2, midpont + inter_len / 2]
    dsum = 0
    n = np.shape(data)[0]
    for j in range(n):
        if method == 'hausdorff':
            dsum += pow(hausdorff_distance(interval, data[j, :]), 2)
        elif method == 'wasserstein':
            dsum += wasserstein_distance2(np.array(interval), data[j, :])
        else:
            dsum += pow(middle_distance(interval, data[j, :]), 2)
    return dsum


def frechet_mean(data, method='hausdorff', inter_l=0):
    """
    :param data: numpy[n,2]: interval data, each line indicates an interval
    :param method: string：'hausdorff': hausdorff distance;
                             'midpoint': middle point;
                             'symbolic': symbolic sample method
                             'arithmetic-based': [E(inf X), E(supX)]
    :param inter_l: bool:  interval length of frechet_mean, None for default, double for specific length
    :return: numpy[2]: frechet_mean
    """
    lower = min(data[:, 0])
    upper = max(data[:, 1])
    inter_len = inter_l if inter_l else np.mean(data[:, 1] - data[:, 0])

    if method == 'hausdorff' or method == 'wasserstein' or method == 'middle':
        output = optimize.fminbound(fdsum, lower, upper, args=(data, inter_len, method), full_output=True)
        return {'min_d': output[1],
                'midpoint': output[0],
                'interval': [output[0] - inter_len / 2, output[0] + inter_len / 2]}

    elif method == 'symbolic' or 'midpoint':
        return np.mean(data)
    elif method == 'arithmetic-based':
        return [np.mean(data[:, 0]), np.mean(data[:, 1])]


def frechet_variance(data, method='hausdorff', theta=1):
    """
    :param data: numpy[n,2]: interval data, each line indicates an interval
    :param method: string：'hausdorff': hausdorff method;
                            'midpoint': middle point;
                            'symbolic': symbolic sample method
                            'arithmetic-based': [E(inf X), E(supX)]
    :param theta: only use in arithmetic-based method. default is 1
    :return: double
    """
    n = data.shape[0]
    if method == 'hausdorff':
        m = frechet_mean(data)
        return m['min_d'] / n
    elif method == 'symbolic':
        s_sum = 0
        for i in range(n):
            a = data[i, 0]
            b = data[i, 1]
            s_sum += (pow(a, 2) + a * b + pow(b, 2))
        return s_sum / (3 * n) - pow(np.sum(data), 2) / (4 * pow(n, 2))
    elif method == 'midpoint':
        return np.var(np.mean(data, axis=1))
    elif method == 'arithmetic-based':
        return np.var(data[:, 1]/2 + data[:, 0]/2) + theta * np.var(data[:, 1]/2 - data[:, 0]/2)


def fcov(w, x, y, x_mean, y_mean):
    return -1 * abs((dist(w, x) - dist(w, x_mean)) * (dist(w, y) - dist(w, y_mean)))


def ffcov(w, x, y, x_mean, y_mean):
    return (dist(w, x) - dist(w, x_mean)) * (dist(w, y) - dist(w, y_mean))


def frechet_covariance(x, y, method='hausdorff', theta=1):
    """
    calculate frechet covariance of random interval x and y using formulation 2

    :param x: numpy[n,2]: interval data, each line indicates an interval
    :param y: numpy[n,2]: interval data, each line indicates an interval
    :param method: string：'hausdorff': hausdorff method;
                           'wasserstein': wasserstein method;
                           'symbolic': symbolic sample method
    :return: double
    """
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    assert x.shape[0] == y.shape[0]
    n = x.shape[0]
    if method == 'hausdorff':
        x_mean = frechet_mean(x)['interval']
        y_mean = frechet_mean(y)['interval']
        cov_sum = 0
        for i in tqdm(range(n)):

            output1 = optimize.brute(fcov, ranges=((x[i, 0], x[i, 1]),), args=(x[i, :], y[i, :], x_mean, y_mean), full_output=True, finish=optimize.fmin)
            output2 = optimize.brute(fcov, ranges=((x_mean[0], x_mean[1]),), args=(x[i, :], y[i, :], x_mean, y_mean), full_output=True, finish=optimize.fmin)
            output3 = optimize.brute(fcov, ranges=((y[i, 0], y[i, 1]),), args=(x[i, :], y[i, :], x_mean, y_mean), full_output=True, finish=optimize.fmin)
            output4 = optimize.brute(fcov, ranges=((y_mean[0], y_mean[1]),), args=(x[i, :], y[i, :], x_mean, y_mean), full_output=True, finish=optimize.fmin)
            output = [output1, output2, output3, output4]
            m_out = max([-1 * i[1] for i in output])
            for j in range(4):
                if output[j][1] * -1 == m_out:
                    max_output = output[j]
            cov_sum += ffcov(max_output[0][0], x[i, :], y[i, :], x_mean, y_mean)
            # a = 0
            # abs_a = 0
            # for w in np.concatenate((np.linspace(start=x[i, 0], stop=x[i, 1], num=1000),
            #                          np.linspace(start=y[i, 0], stop=y[i, 1], num=1000),
            #                          np.linspace(start=x_mean[0], stop=x_mean[1], num=1000),
            #                          np.linspace(start=y_mean[0], stop=y_mean[1], num=1000)
            #                          )):
            #     a0 = (dist(w, x[i, :]) - dist(w, x_mean)) * (dist(w, y[i, :]) - dist(w, y_mean))
            #     if abs(a0) > abs_a:
            #         abs_a = abs(a0)
            #         a = a0
            # cov_sum += abs_a if a >= 0 else (-1 * abs_a)
        return cov_sum / n
    elif method == 'symbolic':
        mean1 = frechet_mean(x, method='symbolic')
        mean2 = frechet_mean(y, method='symbolic')
        c_sum = 0
        for i in range(n):
            a1 = x[i, 0]
            b1 = x[i, 1]
            a2 = y[i, 0]
            b2 = y[i, 1]
            c_sum += (2 * (a1 - mean1) * (a2 - mean2) + (a1 - mean1) * (b2 - mean2) + (b1 - mean1) * (a2 - mean2) + 2 *
                      (b1 - mean1) * (b2 - mean2))
        return c_sum / (6 * n)
    elif method == 'arithmetic-based':
        return np.cov(x[:,0]/2 + x[:,1]/2, y[:,0]/2 + y[:,1]/2)[0][1] + theta * np.cov(x[:,1]/2 - x[:,0]/2, y[:,1]/2 - y[:,0]/2)[0][1]

    # elif method == 'wasserstein':
    #     mean1 = np.array(frechet_mean(x, method='wasserstein')['interval'])
    #     mean2 = np.array(frechet_mean(y, method='wasserstein')['interval'])
    #     m1_values = np.linspace(start=mean1[0], stop=mean1[1], num=10000)
    #     m2_values = np.linspace(start=mean2[0], stop=mean2[1], num=10000)
    #     m1_sorter = np.argsort(m1_values)
    #     m2_sorter = np.argsort(m2_values)
    #     c_sum = 0
    #     for i in range(n):
    #         x_values = np.linspace(start=x[i, 0], stop=x[i, 1], num=10000)
    #         y_values = np.linspace(start=y[i, 0], stop=y[i, 1], num=10000)
    #         x_sorter = np.argsort(x_values)
    #         y_sorter = np.argsort(y_values)
    #
    #         all_values = np.concatenate((m1_values, m2_values, x_values, y_values))
    #         all_values.sort(kind='mergesort')
    #         all_values.sort(kind='mergesort')
    #
    #         # Compute the differences between pairs of successive values of u and v.
    #         deltas = np.diff(all_values)
    #
    #         # Get the respective positions of the values of u and v among the values of
    #         # both distributions.
    #         m1_cdf_indices = m1_values[m1_sorter].searchsorted(all_values[:-1], 'right')
    #         m2_cdf_indices = m2_values[m2_sorter].searchsorted(all_values[:-1], 'right')
    #         x_cdf_indices = x_values[x_sorter].searchsorted(all_values[:-1], 'right')
    #         y_cdf_indices = y_values[y_sorter].searchsorted(all_values[:-1], 'right')
    #
    #         # Calculate the CDFs of u and v using their weights, if specified.
    #
    #         m1_cdf = m1_cdf_indices / m1_values.size
    #         m2_cdf = m2_cdf_indices / m2_values.size
    #         x_cdf = x_cdf_indices / x_values.size
    #         y_cdf = y_cdf_indices / y_values.size
    #         # Compute the value of the integral based on the CDFs.
    #         # If p = 1 or p = 2, we avoid using np.power, which introduces an overhead
    #         # of about 15%.
    #         c_sum += np.sum(np.multiply((x_cdf - m1_cdf) * (y_cdf - m2_cdf), deltas))
    #     return c_sum / n


def frechet_correlation(x, y, method='hausdorff'):
    """
    return the correlation coefficient of x and y using method
    :param x: numpy[n,2]: interval data, each line indicates an interval
    :param y: numpy[n,2]: interval data, each line indicates an interval
    :param method: 'hausdorff': hausdorff method;
                   'wasserstein': wasserstein method;
                   'symbolic': symbolic sample method
    :return: double
    """
    if method == 'midpoint':
        return np.corrcoef(x[:, 1]/2 + x[:, 0]/2, y[:, 1]/2 + y[:, 0]/2)[0][1]
    else:
        return frechet_covariance(x, y, method=method) / \
               np.sqrt(frechet_variance(x, method=method) * frechet_variance(y, method=method))



# def frechet_covariance(x, y, method='hausdorff'):
#     """
#     calculate frechet covariance of random interval x and y
#
#     :param x: numpy[n,2]: interval data, each line indicates an interval
#     :param y: numpy[n,2]: interval data, each line indicates an interval
#     :param method: string：'hausdorff': hausdorff method; 'middle': middle point
#     :return: double
#     """
#     assert x.shape[0] == y.shape[0]
#     n = x.shape[0]
#     if method == 'hausdorff':
#         x_mean = frechet_mean(x)['interval']
#         y_mean = frechet_mean(y)['interval']
#         cov_sum = 0
#         for i in range(n):
#             cov_sum += hausdorff_distance(x[i,], x_mean, sgn=True) * hausdorff_distance(y[i,], y_mean, sgn=True)
#         return {'location covariance': cov_sum / n,
#                 'scale covariance': np.cov(x[:, 1] - x[:, 0], y[:, 1] - y[:, 0])[0][1]}
#     else:
#         return 0


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
        plt.plot(result['interval'], np.array([round(n / 2), round(n / 2)]), alpha=1, color="red", linewidth=2)
    plt.title('plot of interval')  # 折线图标题
    plt.xlabel('Range')  # x轴标题
    plt.ylabel('Number')  # y轴标题
    plt.grid(ls='-.')  # 绘制背景线
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
        a, = plt.plot(x[i, :], np.array([y1[i], y1[i]]), alpha=1, color="red", label='X')
        b, = plt.plot(y[i, :], np.array([y2[i], y2[i]]), alpha=1, color="blue", label='Y')
    if Cov:
        plt.title('plot of ' + r'$X_1$'+' and ' + r'$X_2$' + ' with Cov %.4f' % Cov)  # 折线图标题
    else:
        plt.title('plot of ' + r'$X_1$'+' and ' + r'$X_2$')
    plt.xlabel('Range')  # x轴标题
    plt.ylabel('Number')  # y轴标题
    plt.grid(ls='-.')  # 绘制背景线
    plt.legend(loc='best', handles=[a, b])
    plt.tight_layout()
    plt.savefig(path)
    plt.show()


def show3(x, y, samples=25, path=None, real=None, Cor=None):
    np.random.seed(0)
    n = np.shape(x)[0]
    plt.figure(dpi=200, figsize=(5, 5))
    for i in np.random.choice(n,samples,replace=False):
        plt.fill_between(x[i, :], np.array([y[i, 1], y[i, 1]]), np.array([y[i, 0], y[i, 0]]), alpha=0.25, linewidth=0, color='steelblue')
        if real is not None:
            plt.scatter(real[i, 0], real[i, 1], alpha=1, color="black", s=2)

    if Cor is not None:
        plt.title('plot of ' + r'$X_1$'+' and ' + r'$X_2$' + ' with Cor %.4f' % Cor)  # 折线图标题
    else:
        plt.title('plot of ' + r'$X_1$'+' and ' + r'$X_2$')
    plt.xlabel(r'$X_1$')  # x轴标题
    plt.ylabel(r'$X_2$')  # y轴标题
    plt.grid(ls='-.')  # 绘制背景线
    # plt.legend(loc='best', handles=[a, b])
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()
