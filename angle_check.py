# -*- coding: utf-8 -*-
import numpy as np
# from math import *
from cmath import *
# from sympy import *
w = exp(2 * pi * 1.0j / 3.0)


# 全部转换为横向array，归一化，第三元素正数化


def transform(vector):
    if np.size(vector.dot(vector.conjugate().T)) != 1:
        vector = np.array([vector[i][0] for i in range(len(vector))])
    vector = vector / np.sqrt(vector.dot(vector.conjugate().T))
    fai = np.angle(vector[2])
    vector = vector * exp(-1.0j * fai)
    return vector


def pre_angle_check(vector):

    vector = transform(vector)
    [[th2, phi2, b], [th1, phi1, a]] = pre_angle(vector)
    z = np.array([[0], [0], [1]])
    R13 = np.matrix([[cos(th2 / 2.0), 0, 1.0j * (cos(phi2) - 1.0j * sin(phi2)) * sin(th2 / 2.0)],
                     [0, 1, 0], [1.0j * (cos(phi2) + 1.0j * sin(phi2)) * sin(th2 / 2.0), 0, cos(th2 / 2.0)]])
    R23 = np.matrix([[1, 0, 0], [0, cos(th1 / 2.0), 1.0j * (cos(phi1) - 1.0j * sin(phi1)) * sin(th1 / 2.0)],
                     [0, 1.0j * (cos(phi1) + 1.0j * sin(phi1)) * sin(th1 / 2.0), cos(th1 / 2.0)]])
    zer0_check = R23 * R13 * z - np.matrix(vector).T
    check = np.abs(zer0_check > 1e-8) * 1.0
    return check.T


def pre_angle(vector):
    ''' get rotation angles from a state vector

    [1]
    [2]
    [3]


    作用顺序为R23(a, phi1)*R13(b, phi2)*[0,0,1]=vector
    '''
    # 变形，归一化，第三分量去掉相位
    vector = transform(vector)

    # vector格式为1*3的态矢矩阵,能够解第三分量为实数的矢量
    if np.abs(np.abs(vector[0]) - 1.0) < 1e-14:
        # 如果矢量为[1，0，0]，直接给出结果，避免奇异
        [a, phi1, b, phi2] = [0, 0, pi, pi / 2.0]
    else:
        # 以前的版本，其实是一样的
        # a = 2 * asin(np.abs(vector[1]) / sqrt(1 - np.abs(vector[0]) ** 2))
        # phi1 = pi / 2.0 - np.angle(vector[1])
        # b = 2 * asin(np.abs(vector[0]))
        # phi2 = pi / 2.0 - np.angle(vector[0])
        a = 2 * acos(np.abs(vector[2]) / sqrt(1 - np.abs(vector[0]) ** 2))
        phi1 = pi / 2.0 - np.angle(vector[1])
        b = 2 * asin(np.abs(vector[0]))
        phi2 = pi / 2.0 - np.angle(vector[0])

    return [(b, phi2, 0), (a, phi1, 1)]


def GenStablizerVectors(X):
    ''' 生成X的本征向量

    按照 1, w， w＊w 本征值排列本征向量
    '''
    w, v = np.linalg.eig(X)
    a = np.exp(-2.0j * pi / 3)
    index = [i for iw in [w, w * a, w * a * a]
             for i, e in enumerate(iw) if np.abs(e - 1) < 1e-15]
    print index
    vs = [np.array(v[:, i].T)[0] for i in index]
    return vs


if __name__ == "__main__":
    X = np.matrix([[0, 0, 1.0], [1.0, 0, 0], [0, 1.0, 0]])
    Z = np.matrix([[1, 0, 0], [0, np.exp(2.0j * pi / 3.0), 0], [0, 0, np.exp(4.0j * pi / 3.0)]])
    a = np.array([[1.0j - 1], [1 - 2j], [3 + 4j]])
    a = GenStablizerVectors(X) + GenStablizerVectors(Z)
    print [pre_angle_check(i) for i in a]
    print np.angle(-1)
