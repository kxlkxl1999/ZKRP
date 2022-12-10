# -*- coding: utf-8 -*-
import numpy as np

# from tqdm import tqdm
# import matplotlib.pyplot as plt


def hausdorff_dist(x, y, cr=False):
    """
    Hausdorff distance between two intervals
    :param x: array-like(2): interval data
    :param y: array-like(2): interval data
    :cr: whether input is center+range,
        if false then inputs the lower and upper bound
    """
    if cr:
        xl = x[0] + x[1]
        xu = x[0] - x[1]
        yl = y[0] + y[1]
        yu = y[0] - y[1]
    else:
        xl, xu = x[0], x[1]
        yl, yu = y[0], y[1]

    return min(abs(xl - yl), abs(xu - yu))


def hausdorff_dist_vec(x, y, cr=False):
    """
    Hausdorff distance between two sets of intervals
    which is the vectorized version of hausdorff_dist
    :param x: array-like (n, 2)
    :param y: array-like (n, 2)
    :cr: whether input is center+range,
        if false then inputs the lower and upper bound
    """
    if cr:
        xl = x[:, 0] - x[:, 1]
        xu = x[:, 0] + x[:, 1]
        yl = y[:, 0] - y[:, 1]
        yu = y[:, 0] + y[:, 1]
    else:
        xl, xu = x[:, 0], x[:, 1]
        yl, yu = y[:, 0], y[:, 1]

    return np.min(np.vstack([abs(xl - yl), abs(xu - yu)]), 0)


def arithmetic_based_distance(x, y, theta=1):
    return np.sqrt(
        pow(np.mean(x) - np.mean(y), 2)
        + theta * pow(x[1] / 2 - x[0] / 2 - y[1] / 2 + y[0] / 2, 2)
    )
