
# -*-coding:UTF-8-*-
"""

Created on Fri July 29 2016

@author: heran
"""
#from TomographyQubit import *
#from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
from sympy import *
from sympy.tensor.array import Array, tensorproduct
from numpy import mgrid
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors, ticker, cm
init_printing()

w = exp(2 * pi * 1.0j / 3.0)
III = eye(3)
X = Matrix([[0, 0, 1.0], [1.0, 0, 0], [0, 1.0, 0]])
Z = Matrix([[1.0, 0, 0], [0, w, 0], [0, 0, w**2]])
D = [(w**((2 * i * j) % 6)) * (X**i) * (Z**j) for i in range(3) for j in range(3)]

var('th1 th2 phi1 phi2')
# RZ=Matrix([[cos(a/2.0),-sin(a/2.0),0],[sin(a/2.0),cos(a/2.0),0],[0,0,1]])
# RY=Matrix([[cos(th2/2.0),0,(cos(phi2)+1.0j*sin(phi2))*sin(th2/2.0)],[0,1,0],[(-cos(phi2)+1.0j*sin(phi2))*sin(th2/2.0),0,cos(th2/2.0)]])
# RX=Matrix([[1,0,0],[0,cos(th1/2.0),(-cos(phi1)+1.0j*sin(phi1))*sin(th1/2.0)],[0,(cos(phi1)+1.0j*sin(phi1))*sin(th1/2.0),cos(th1/2.0)]])
RY = Matrix([[cos(th2 / 2.0), 0, 1.0j * (cos(phi2) - 1.0j * sin(phi2)) * sin(th2 / 2.0)],
             [0, 1, 0], [1.0j * (cos(phi2) + 1.0j * sin(phi2)) * sin(th2 / 2.0), 0, cos(th2 / 2.0)]])
RX = Matrix([[1, 0, 0], [0, cos(th1 / 2.0), 1.0j * (cos(phi1) - 1.0j * sin(phi1)) * sin(th1 / 2.0)],
             [0, 1.0j * (cos(phi1) + 1.0j * sin(phi1)) * sin(th1 / 2.0), cos(th1 / 2.0)]])


for i in range(9):
    D[i] = np.array(D[i])


def normalize(v):
    norm = sqrt(N(np.conj(v.T) * v, chop=True)[0])
    if (norm) == 0:
        return v
    else:
        return v / S(norm)


def projectors_of_XYZV(X, Z):
    # give all seperate projectors
    # D={Z,X,wwXZ,wwwwXZZ}
    a = 1 / sqrt(3)
    z0 = Matrix([1, 0, 0])
    z1 = Matrix([0, 1, 0])
    z2 = Matrix([0, 0, 1])
    #[eigenvalue,degeneration,eigenvects]
    eig0 = [(1, 1, [z0]), (w, 1, [z1]), (w**2, 1, [z2])]
    eig1 = (X).eigenvects()
    eig2 = ((w**2) * X * Z).eigenvects()
    eig3 = ((w**4) * X * Z * Z).eigenvects()
    eig = [eig0, eig1, eig2, eig3]
    # change the order of elements as w**0,w**1,w**2...
    for i in range(1, 4):
        exchange = eig[i][2]
        eig[i][2] = eig[i][1]
        eig[i][1] = exchange
    # normalise
    for i in range(4):
        for j in range(3):
            eig[i][j][2][0] = normalize(eig[i][j][2][0])

    # obtain all projectors from all vectors
    proj1 = [[zeros(3) for i in range(3)] for j in range(4)]
    for i in range(4):
        for j in range(3):
            proj1[i][j] = expand(eig[i][j][2][0] * np.conj(eig[i][j][2][0].T))
    return proj1, eig


def all_facet_Aq(proj1):
    # totally 9 q
    arr_a = np.array([1, 0, 1, 2])
    arr_b = np.array([0, 1, 1, 1])

    q = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    for i in range(3):
        for j in range(3):
            q[i][j] = (i * arr_a - j * arr_b) % S(3)

    # get 9 facets
    Aq = [[zeros(3) for i in range(3)]for j in range(3)]
    for a in range(3):
        for b in range(3):
            for i in range(4):
                for j in range(3):
                    if (j) == q[a][b][i]:
                        Aq[a][b] = N(Aq[a][b] + proj1[i][j], 10, chop=True)
            Aq[a][b] = N(Aq[a][b] - eye(3), 10, chop=True)
    return Aq, q


# get wigner function with some unknow elements x,y,c
def Wignerfunc(Aq, rou):
    W = zeros(3)
    for i in range(3):
        for j in range(3):
            W[i, j] = (Aq[i][j] * rou).trace() / 3.0
    W = simplify(W)
    return W


def draw_zone(Aq, rouxyz, x_point, y_point):
    # use lanbdify to subs var elements a,b,to varibles x,y
    area = (rouxyz * rouxyz).trace()
    # area=N(expand(area),5)

    # simulable area, positive wigner function
    sim = zeros(3)
    for i in range(3):
        for j in range(3):
            sim[i, j] = (Aq[i][j] * rouxyz).trace()
    sim = expand(sim)

    # dett=simplify(rouxyz.det())
    fdet = lambdify((a, b), rouxyz.det(), 'numpy')
    ftrace = lambdify((a, b), area, 'numpy')
    print rouxyz.det()
    print area
    # plot setting
    sns.set_style("whitegrid", {'axes.edgecolor': 'black', 'legend.frameon': True})
    sns.set_context("notebook", font_scale=1.4, rc={"lines.linewidth": 3})

    minaxis = -0.5
    maxaxis = 0.5
    resolution = 0.01
    x = np.arange(minaxis, maxaxis, resolution)
    y = np.arange(minaxis, maxaxis, resolution)
    x, y = np.meshgrid(x, y)

    det_aera = fdet(x, y)
    trace_area = ftrace(x, y)
    # simulable area is input by hand, it should change with wigner function slice
    sim00 = -3 * x - 3 * y - 3 * c + 4 / 3.0
    sim02 = y
    sim12 = x

    phasepoint = (x - x_point)**2 + (y - y_point)**2

    fig1 = plt.subplots(figsize=(9, 9))

    sns.plt.contourf(x, y, det_aera, [0, 1], colors='b', alpha=0.4)  # inside area
    sns.plt.contour(x, y, det_aera, [0], colors='r')  # edge

    sns.plt.contour(x, y, trace_area, [1], colors='g', alpha=0.4)  # pure state edge
    sns.plt.contour(x, y, phasepoint, [0.00005], colors='b', alpha=0.7)  # phase point to measure
    # simulable edge
    sns.plt.contour(x, y, sim00, [0], colors='g', alpha=0.4)
    sns.plt.contour(x, y, sim02, [0], colors='g', alpha=0.4)
    sns.plt.contour(x, y, sim12, [0], colors='g', alpha=0.4)

    plt.axis('scaled')
    plt.show()

'''
# fix a point (x,y),but unnecessaryly to be here
x_point = 0.15
y_point = 0.05
print 'x_point=%s,y_point=%s ' % (x_point, y_point)
'''

if __name__ == '__main__':
    # define a common rou
    var('a1 a2 b1 b2 c1 c2 d1 d2')
    rou = Matrix([[a1, b1 - 1.0j * b2, c1 - 1.0j * c2],
                  [b1 + 1.0j * b2, a2, d1 - 1.0j * d2],
                  [c1 + 1.0j * c2, d1 + 1.0j * d2, 1 - a1 - a2]])

    c = 1 / 9.0  # ------an third variable of wigner function
    proj1, eig = projectors_of_XYZV(X, Z)  # ------get all stabilizer projectors
    Aq, q = all_facet_Aq(proj1)  # ------get all 9 facets
    Wigner = Wignerfunc(Aq, rou)  # ------get wigner function of the common rou

    var('a b')
    # ------define a special wigner function, here we have mms=maximum mixed slice
    mms = Matrix([[4 / 9.0 - a - b - c, 1 / 9.0, a], [1 / 9.0, 1 / 9.0, b], [1 / 9.0, 1 / 9.0, c]])
    # ------solve (common wigner function-special wigner function)=0, and get solutions to 6 variables
    re = solve(Wigner - mms, [a1, a2, b1, b2, c1, c2, d1, d2])
    # ------insert the special solutions in rouxyz,and get a special rou corresponding to mms
    rouxyz = Matrix([[re[a1], re[b1] - 1.0j * re[b2], re[c1] - 1.0j * re[c2]],
                     [re[b1] + 1.0j * re[b2], re[a2], re[d1] - 1.0j * re[d2]],
                     [re[c1] + 1.0j * re[c2], re[d1] + 1.0j * re[d2], 1 - re[a1] - re[a2]]])
    # set the remaining 2 variables x,y
    print rouxyz

    '''
    # ------replace a,b with x_point,y_point to get a definite rou, get its 3 eigenvecters to measure
    eigenmixed = (rouxyz.subs('a', x_point).subs('b', y_point)).eigenvects()
    # ------the 3 eigenvalues are the mixed probabilities of the 3 eigenvecters
    mixprobability = [N(eigenmixed[i][0], chop=True) for i in range(3)]

    mixpurestate = [normalize(eigenmixed[i][2][0]) for i in range(3)]  # ------regroup 3 eigenvecters
    pre_angle_list = [pre_angle(mixpurestate[i]) for i in range(3)]  # ------get their preparation pluse lists

    # vectors=[[eig[i][j][2][0] for j in range(3)] for i in range(4)]
    # detect_angle_list=[[detect_angle(vectors[i][j]) for j in range(3)]for i in range(4)]
    checkresult = [angle_check(s, pre_angle_list, mixpurestate, RX, RY) for s in range(3)]  # ------check
    print '-------checkresult---------'
    print checkresult
    # print '-------pre_angle_list---------'
    # print pre_angle_list
    # print pre_angle(mixpurestate[0])
    # print pre_angle(mixpurestate[1])
    # print mixpurestate[1]

    data = (proj1, Aq, q, eigenmixed, mixprobability, mixpurestate, pre_angle_list, rouxyz, rou)
    # dataN=[N(x) for x in data]

    # output all result to a pkl file for future use
    import pickle as p
    import pprint
    wignerf = r'D:\Users\Desktop\HERAN\IonControl\IonControl\NCHV\wignerfunction1\wignerfunctiontopulselistx-%.6fy-%.6f.pkl' % (x_point, y_point)
    with open(wignerf, 'w') as f:
        p.dump(data, f)

    # here,rouxyz should contain variables a and b
    draw_zone(Aq, rouxyz, x_point, y_point)
    '''
    data = [Aq, q]
    import pickle as p
    import pprint
    wignerf = r'negativity_9Aq.pkl'
    with open(wignerf, 'w') as f:
        p.dump(data, f)
    x_point, y_point = [0.0, 0.0]
    draw_zone(Aq, rouxyz, x_point, y_point)
