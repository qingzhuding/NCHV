# -*- coding: utf-8 -*-

import numpy as np
from sympy import *
import pickle as p
from itertools import chain

# 归一化


def normalize(v):
    norm = sqrt(N(np.conj(v.T) * v, chop=True)[0])
    if (norm) == 0:
        return v
    else:
        return v / S(norm)


if __name__ == "__main__":

    wignerf = r'C:\Users\R.He\Documents\testing\IonControl\NCHV\ent_proj_decomposition-%s.pkl' % 216
    with open(wignerf, 'r') as f:
        (III, X, Z, D, Clif3, Clifford3, psi3, list_entproj, eig,
            project1, project2, Basis_1, Basis_2, Basis, decompose_all, rou_n, ind_k) = p.load(f)

    # input state
    state11 = Matrix([1, -2.0j, 0])
    state11 = Matrix([0.77036436, 0.10669539 - 0.21612418j, 0.29514623 - 0.51120827j])
    # state11=Matrix([cos(pi/4.0),0,sin(pi/4.0)])
    state22 = Matrix([0.77036436, 0.10669539 - 0.21612418j, 0.29514623 - 0.51120827j])
    state11 = normalize(state11)
    state22 = normalize(state22)
    state33 = Matrix(np.kron(state11, state22))

    # get correspongding probabilities
    proba1 = np.array([abs((state11.conjugate().T * p1)[0])**2 for p1 in project1]).reshape(3, 8)
    proba2 = np.array([abs((state22.conjugate().T * p1)[0])**2 for p1 in project2]).reshape(3, 4)

    proba1 = np.array([[0.60033974, 0.60033974, 0.53991597, 0.53991597, 0.47829132, 0.47959184, 0.04121649, 0.28571429],
                       [0.05492197, 0.34552621, 0.37434974, 0.06372549, 0.47959184, 0.47829132, 0.6552621, 0.6552621],
                       [0.3452621, 0.05492197, 0.06372549, 0.37434974, 0.0120048, 0.0120048, 0.28571429, 0.04121649]])

    proba2 = np.array([[0.6007403, 0.53841537, 0.04251701, 0.47659064],
                       [0.05222089, 0.06592637, 0.65506202, 0.47909164],
                       [0.34372949, 0.3787515, 0.30751501, 0.0110044]])
    proba1[2, :] = 1 - proba1[1, :] - proba1[0, :]
    proba2[2, :] = 1 - proba2[1, :] - proba2[0, :]

    print 'probability of 1 in theory'
    print N(Matrix(proba1), 5, chop=True)
    print 'probability of 2 in theory'
    print N(Matrix(proba2), 5, chop=True)
    # proba simulation check 按照输入态所对应的的概率制备0，1，带入前9个ent_proj，结果应该为1./////////////////////////////
    BES = []
    record = []
    for j in range(5000):
        RE1 = [[np.random.choice([0, 1], p=[1 - proba1[a][b], proba1[a][b]]) for b in range(8)] for a in range(3)]
        RE2 = [[np.random.choice([0, 1], p=[1 - proba2[a][b], proba2[a][b]]) for b in range(4)] for a in range(3)]
        record.append(list(chain.from_iterable(RE2)))
        # 所有基单次测量值，0、1，3*8矩阵，3*4矩阵
        single9 = []  # 9个ent_proj的单次值
        for n in range(216):  # 这里换成rou_n，全部为216
            # 下面一行可获得单次测量结果值，需要累积平均
            a = np.array(([-1] + [RE1[(ind_k[n][M][1] - i) % 3][ind_k[n][M][0] - 1] * RE2[i][M]
                                  for M in range(4) for i in range(3)])) / 3.0
            single9.append(np.sum(a))  # 也可以输出其他形式，这里直接加起来得到最终值
        BES.append(single9)
    print 'sampling average simulation'
    print np.average(record, axis=0)
    # rank—1投影的概率值
    prob_simulation = N(Matrix(np.average(BES, axis=0)).T, 5, chop=True)

    # trace 直接用ent_proj求的概率
    prob_list = Matrix([N(((state33 * state33.conjugate().T) *
                           list_entproj[i]).trace(), 5, chop=True) for i in range(9)]).T
    print 'prob_simulation'
    print prob_simulation
    print 'sum'
    print np.sum(prob_simulation)
    print prob_list, np.sum(prob_list)
    # output all result to a pkl file for future use
