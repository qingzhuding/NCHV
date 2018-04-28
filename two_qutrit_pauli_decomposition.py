# -*- coding: utf-8 -*-

import numpy as np
from sympy import *
import pickle as p

init_printing()
w1 = exp(2.0j / 3.0 * pi)
III = eye(3)
X = Matrix([[0, 0, 1.0], [1.0, 0, 0], [0, 1.0, 0]])
# v=Matrix([[0],[0],[1]])
Z = Matrix([[1.0, 0, 0], [0, w1, 0], [0, 0, w1**2]])

# F matrix//////////////////////////////////////////////////////////////////////////////////////////////////////
f1 = np.array([[0, 1], [2, 0]])
f2 = np.array([[0, 1], [2, 1]])
f3 = np.array([[0, 1], [2, 2]])

f4 = np.array([[0, 2], [1, 0]])
f5 = np.array([[0, 2], [1, 1]])
f6 = np.array([[0, 2], [1, 2]])

f7 = np.array([[1, 1], [2, 0]])
f8 = np.array([[1, 2], [1, 0]])

f9 = np.array([[1, 0], [0, 1]])
f10 = np.array([[1, 0], [1, 1]])
f11 = np.array([[1, 0], [2, 1]])

f12 = np.array([[1, 1], [0, 1]])
f13 = np.array([[1, 2], [0, 1]])

f14 = np.array([[1, 1], [1, 2]])
f15 = np.array([[1, 2], [2, 2]])

f16 = np.array([[2, 1], [2, 0]])
f17 = np.array([[2, 2], [1, 0]])

f18 = np.array([[2, 1], [1, 1]])
f19 = np.array([[2, 2], [2, 1]])

f20 = np.array([[2, 0], [0, 2]])
f21 = np.array([[2, 0], [1, 2]])
f22 = np.array([[2, 0], [2, 2]])

f23 = np.array([[2, 1], [0, 2]])
f24 = np.array([[2, 2], [0, 2]])

# thisSL(2,Z)群, F[n]表示f(n+1)群元素。F[n][a,b]表示群元素f[n+1]的矩阵中的矩阵元，如F[2][0,0]
F = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17, f18, f19, f20, f21, f22, f23, f24]
UF = [np.array(zeros(3)) for i in range(24)]

# t=S(((-1+I*sqrt(3))/2)**2)

for i in range(24):
    if(F[i][0, 1]) == 0:
        for k in range(3):
            alpha = F[i][0, 0]
            gama = F[i][1, 0]
            ark = ((alpha * gama * k * k) * 2) % 3
            alpha = (alpha * k) % 3
            UF[i][alpha, k] = w1**ark

    if(F[i][0, 1]) != 0:
        cons = 1 / sqrt(3)
        for j in range(3):
            for k in range(3):
                alpha = F[i][0, 0]
                beta = F[i][0, 1]
                delta = F[i][1, 1]
                ark = ((beta * (alpha * k * k - 2 * j * k + delta * j * j)) * 2) % 3
                UF[i][j, k] = w1**ark
        UF[i] = cons * UF[i]


# entangled projectors///////////////////////////////////////////////////////////////////////////////////////
D = [(N(w1**((2 * i * j) % 3)) * (X**i) * (Z**j)) for i in range(3) for j in range(3)]
for i in range(9):
    D[i] = np.array(D[i])
# print(X.dot(Z))
Clif3 = [D[i].dot(U_e) for U_e in UF for i in range(9)]
Clifford3 = [np.kron(ci, III) for ci in Clif3]
print len(Clifford3)

# [N(Matrix(i)*Matrix([1,1,1])/np.sqrt(3.0),5,chop=True) for i in Clif3]
# BES0
psi1 = (np.kron(III, III) + np.kron(X, X) + np.kron(X * X, X * X))
psi2 = (np.kron(III, III) + np.kron(Z, Z * Z) + np.kron(Z * Z, Z))
psi3 = N(Matrix(psi1) * Matrix(psi2) / 9.0, chop=True)
# print psi3
# get all 216 ent-projs
list_entproj = [N(Matrix(ci) * psi3 * Matrix(ci.T.conj()), 7, chop=True) for ci in Clifford3]


# seperate projectors/////////////////////////////////////////////////////////////////////////////////////////
# D={Z,X,wwXZ，wwwwwXZZ}
z0 = Matrix([1, 0, 0])
z1 = Matrix([0, 1, 0])
z2 = Matrix([0, 0, 1])
# 下述本征函数组中，如eig1中，为三乘三的矩阵，第一行为本征值为1的【本征值，简并度，本征矢量（3*1矩阵）】，
# 第二行本征值为w**2，第三行为w**1.
eig0 = [[(1, 1, [z0]), (N(w1**2), 1, [z1]), (N(w1), 1, [z2])]] + [[(1, 1, [z0]), (N(w1**2), 1, [z2]), (N(w1), 1, [z1])]]
Bas1 = [X, X * X]
Bas2 = [X * Z * Z, X * X * Z, X * Z, X * X * Z * Z]  # ,Z*Z,Z,
eig = [a.eigenvects() for a in Bas1] + eig0 + [a.eigenvects() for a in Bas2]
for i in range(8):
    asd = eig[i][2]
    eig[i][2] = eig[i][1]
    eig[i][1] = asd


# 归一化
def normalize(v):
    norm = sqrt(N(np.conj(v.T) * v, chop=True)[0])
    if (norm) == 0:
        return v
    else:
        return v / S(norm)


project1 = []
project2 = []
list1 = []  # 用于检查pronjector的顺序
list2 = []
for j in range(3):
    for i in range(8):
        eig[i][j][2][0] = normalize(eig[i][j][2][0])
        project1.append(eig[i][j][2][0])
        list1.append([i, j])
for j in range(3):
    for i in [0, 3, 6, 4]:  # 这样定义的projector与tomography中的格式不同，reshape并转置关系
        project2.append(eig[i][j][2][0])
        list2.append([i, j])
# list1,list2


# 获得rou分解到81个基的9*9系数矩阵，约2s /个///////////////////////////////////////////////////////////////////////
Basis_1 = [N(III), N(X), N(X * X), N(Z * Z), N(Z), N(X * Z * Z), N(X * X * Z), N(X * Z), N(X * X * Z * Z)]
Basis_2 = [N(III), N(X), N(X * X), N(Z), N(Z * Z), N(X * Z), N(X * X * Z * Z), N(X * Z * Z), N(X * X * Z)]
Basis = [[Matrix(np.kron(i, j)) for i in Basis_1] for j in Basis_2]

decompose_all = []
for k in range(216):  # 最大为216
    # 注意这里的Basis[j][i]顺序，为的是与基阵Basis一一对应，而且加上共轭项，输出对应基的系数。
    a = Matrix([[N((list_entproj[k] * (Basis[j][i].conjugate().T)).trace(), 6, chop=True)
                 for i in range(9)] for j in range(9)])
    decompose_all.append(a)
rou_n = len(decompose_all)
print 'number of BES %s' % rou_n

# get the 4*2matrix of rous


def substi_IE(m):  # IE表示indice and element，可以完全表示rou，输出indice,k=（-c）mod3
    indice = np.nonzero(m)  # 获取非0元素位置
    index = [1, 3, 5, 7]  # 只提取四个矩阵和数
    b = np.array([indice[1][i] for i in index])
    a = [np.array(m)[indice[0][i]][indice[1][i]] for i in index]
    a = np.where(np.array(a) == (N(1.0, 6, chop=True)), 0, a)  # 这一步放在第一步，防止二次修改
    a = np.where(np.array(a) == (N(w1, 6, chop=True)), 1, a)
    a = np.where(np.array(a) == (N(w1**2, 6, chop=True)), 2, a)
    a = b, (-a.T) % 3  # 提出c，获得k
    a = np.array(a).T
    return a


if __name__ == "__main__":
    # 基分解basis decomposition ,indice,k
    ind_k = [substi_IE(i) for i in decompose_all]
    # [Matrix(i) for i in ind_k]

    data = (III, X, Z, D, Clif3, Clifford3, psi3, list_entproj, eig, project1,
            project2, Basis_1, Basis_2, Basis, decompose_all, rou_n, ind_k)

    ent_proj_decomposition = r'C:\Users\R.He\Documents\testing\IonControl\NCHV\ent_proj_decomposition-%s.pkl' % rou_n
    with open(ent_proj_decomposition, 'w') as f:
        p.dump(data, f)
