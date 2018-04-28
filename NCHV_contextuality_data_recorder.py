# -*- coding: utf-8 -*-
"""
Created on Dec 13 2016

@author: jmcui
"""


def GenStablizerList(RabiFre, freq1=200, freq2=207.6372):
    # RabiFre unit MHz
    # freq1 MHz,#freq2 MHz
    a = 0.39182655
    a3 = 1 / 3.0
    a2 = 1 / 2.0
    #[theta1,phai1,theta2,phai2]
    # R2*R1*STATE=[0,0,1].先打入脉冲R1[theta2,phai2]，再打入脉冲R2[theta1,phai1].
    # 所以value应该是[（freq2,theta1,phai1）,（freq2,theta2,phai2）]
    z0 = [0, 0, 1, 1]
    z1 = [1, 0, 0, 0]
    z2 = [0, 0, 0, 0]
    x0 = [a, 0, a2, 1]
    x1 = [a, a3 + 1, a2, a3]
    x2 = [a, -a3 + 1, a2, -a3]
    y0 = [a, 0, a2, -a3]
    y1 = [a, a3 + 1, a2, 1]
    y2 = [a, -a3 + 1, a2, a3]
    v0 = [a, 0, a2, a3]
    v1 = [a, a3 + 1, a2, -a3]
    v2 = [a, -a3 + 1, a2, 1]
    #
    Angle = [[z0, z1, z2], [x0, x1, x2], [y0, y1, y2], [v0, v1, v2]]

    value = [[(freq1 * 2 * pi, Angle[i][j][2] / RabiFre / 2.0, Angle[i][j][3] * pi), (freq2 * 2 *
                                                                                      pi, Angle[i][j][0] / RabiFre / 2.0, Angle[i][j][1] * pi)] for i in range(4) for j in range(3)]
    # colmn and lines have Checked
    return value
