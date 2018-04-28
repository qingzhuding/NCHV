# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2016

@author: jmcui
"""
from TomographyQutritCoolingValiditytest import *
import scipy.linalg as linalg
import pickle as p
from sympy import var
import matplotlib.pyplot as plt


def GenStablizerList(RabiTime1, RabiTime2, freq1, freq2):
    # RabiFre unit MHz
    # freq1 MHz,#freq2 MHz
    a = 0.39182655
    b = 0.5
    c = 1 / 6.0
    # [theta1,phai1,theta2,phai2]
    # R2*R1*STATE=[0,0,1].先打入脉冲R1[theta2,phai2]，再打入脉冲R2[theta1,phai1].所以value应该是[（freq2,theta1,phai1）,（freq2,theta2,phai2）]
    # now
    z0 = [0, 0, 1, b]
    z1 = [1, b, 0, 0]
    z2 = [0, 0, 0, 0]
    x0 = [b, b, a, b]
    x1 = [b, -c, a, 1 + c]
    x2 = [b, 1 + c, a, -c]
    y0 = [b, b, a, -c]
    y1 = [b, -c, a, b]
    y2 = [b, 1 + c, a, 1 + c]
    v0 = [b, b, a, 1 + c]
    v1 = [b, -c, a, -c]
    v2 = [b, 1 + c, a, b]

    Angle = [[z0, z1, z2], [x0, x1, x2], [y0, y1, y2], [v0, v1, v2]]

    value = [[(freq2 * 2 * pi, Angle[i][j][0] * RabiTime2 / 2.0, pi + Angle[i][j][1] * pi),
              (freq1 * 2 * pi, Angle[i][j][2] * RabiTime1 / 2.0, pi + Angle[i][j][3] * pi)]
             for i in range(4) for j in range(3)]
    # colmn and lines have Checked
    return value


def StatePrepareThenTomogra(prepareList, RabiTime1, RabiTime2, freq1, freq2):  # freq2=207.6372
    # RabiTime unit is us, flip 2*pi phase ,drive |0> state to -|0> state
    MeasureLists = GenStablizerList(RabiTime1, RabiTime2, freq1, freq2)
    pulse_list, TimeTagIndex = CombinePrepareMeasureLists(prepareList, MeasureLists)

    waves, Tags = GenMultiWaves(pulse_list)  # 返回所有脉冲的timetag
    t, dt = zip(*Tags)  # t为起始时间，dt为脉冲时长
    TimeTags = [(t[TimeTagIndex[i]], sum(dt[TimeTagIndex[i]:TimeTagIndex[i + 1]]))
                for i in range(len(MeasureLists))]  # [(time,pulse_len) for i in N]
    return waves, TimeTags


def solveX_rou_tomo(r2):  # 根据测量值求出探测的密度矩阵,现在并没有用到，而是直接用下面的函数写了出来.
    # b1,2, c1,2,d1,2是密度矩阵里面的参量

    var('a1 a2 b1 b2 c1 c2 d1 d2')
    w1 = 1 / 3.0
    w2 = 1 / sqrt(3.0)
    v1 = 2 * w1 * (b1 + c1 + d1) + w1
    v2 = -w1 * (b1 + c1 + d1) - w2 * (b2 - c2 + d2) + w1
    # v3=-w1*(b1+c1+d1)+w2*(b2-c2+d2)+w1#与上面两个线性相关
    v4 = -w1 * (b1 + c1 - 2 * d1) - w2 * (b2 + c2) + w1
    v5 = -w1 * (b1 + d1 - 2 * c1) + w2 * (b2 - d2) + w1
    # v6=-w1*(d1+c1-2*b1)+w2*(d2+c2)+w1#与上面两个线性相关
    v7 = -w1 * (b1 + c1) + w2 * (b2 + c2) + 2 * w1 * d1 + w1
    v8 = -w1 * (c1 + d1) + 2 * w1 * b1 - w2 * (c2 + d2) + w1
    # v9=0。。。
    a1 = r2[0]
    a2 = r2[1]

    r3 = np.array([r2[3], r2[4], r2[6], r2[7], r2[9], r2[10]])  # 取所需要的8个测量值
    r1 = np.array([v1, v2, v4, v5, v7, v8])
    solveX = simplify(solve(r1 - r3, b1, b2, c1, c2, d1, d2))
    rou_tomo = np.array([[a1, solveX[b1] - I * solveX[b2], solveX[c1] - I * solveX[c2]],
                         [solveX[b1] + I * solveX[b2], a2, solveX[d1] - I * solveX[d2]],
                         [solveX[c1] + I * solveX[c2], solveX[d1] + I * solveX[d2], 1 - a1 - a2]])
    rou_tomo = np.array(rou_tomo).astype(np.complex)
    return rou_tomo

# 把solveX_rou_tomo函数给直接求出来，输出，不需要求解


def get_rou_tomo(detectedprob):
    '''tomography with wigner function basis
        given detected probalities of 12 unbiased basis
        return rou_tomo
    '''
    v1 = sqrt(3) / 6.0
    a1, a2, a3, m1, m2, m3, m4, m5, m6, m7, m8, m9 = detectedprob
    b1 = 0.5 * (m1 - m4 - m5 + m8)
    b2 = v1 * (-m1 - 2 * m2 - m4 + m5 + 2 * m7 + m8)
    c1 = 0.5 * (m1 - m7 + m5 - m8)
    c2 = v1 * (m1 + 2 * m2 - 2 * m4 - m5 + m7 - m8)
    d1 = 0.5 * (m1 + m4 + m7) - 0.5
    d2 = v1 * (-m1 - 2 * m2 - m4 - 2 * m5 - m7 - 2 * m8 + 3)
    rou_tomo = np.array([[a1, b1 - I * b2, c1 - I * c2], [b1 + I * b2, a2, d1 - I * d2],
                         [c1 + I * c2, d1 + I * d2, 1 - a1 - a2]])
    rou_tomo = np.array(rou_tomo).astype(np.complex)
    return rou_tomo


if __name__ == "__main__":  # -----//------流程图详解见下

    CoolingTime = 4  # ms
    PumpingTime = 400  # us
    DetectingTime = 1000  # us#-----//------设置冷却，泵浦，探测时间，可以根据实际需求进行优化

    x_point, y_point = (-0.11, 0.2)

    wignerf = r'D:wignerfunctiontopulselistx-%.6fy-%.6f.pkl' % (x_point, y_point)
    with open(wignerf, 'r') as f:
        proj1, Aq, q, eigenmixed, mixprobability, mixpurestate, pre_angle_list, rouxyz, rou = p.load(f)
    # print proj1
    # print Aq[0][0]
    rouxyz = rouxyz.subs('a', x_point).subs('b', y_point)  # 得到赋值之后的密度矩阵
    state_numb = len(mixpurestate)  # 其实就是3
    # vector to proj.  mixpurestate是态,mixpurerou是密度矩阵
    mixpurerou = [np.array(mixpurestate[i] * mixpurestate[i].conjugate().T).astype(np.complex)
                  for i in range(state_numb)]
    rouxyz = np.array(rouxyz).astype(np.complex)
    mixprobability = np.array(mixprobability).astype(np.float)
    pre_angle_list = [np.array(i).astype(np.float) for i in pre_angle_list]

    def fidelity_func(detec_rho, rho):
        a = linalg.sqrtm(detec_rho)
        return np.trace(linalg.sqrtm(np.dot(np.dot(a, rho), a)))

# print prepareList
    N = 20
    rho_tomo = np.zeros((state_numb, state_numb), dtype=np.complex)
    rho_threepart = []
    fidelity_list = []
    proba_list = []

    def test_fidelity(prepareList, mixpurerou):
        waves, timeTag = StatePrepareThenTomogra(prepareList, RabiTime1, RabiTime2, freq1, freq2)
        setChannelWave32K(2, waves, amplitude=1.0)
        t0, dt = zip(*timeTag)
        t0 = np.array(t0) + 1.00  # -----//-----compesation delay of AFG, we should triger the AFG in advance
        timeTag = zip(t0, dt)
        InstListData = GenInstList(timeTag, CoolingTime, PumpingTime, DetectingTime)
        k = len(timeTag)  # 输入脉冲数
        Task = SpinTask(InstListData, k * N * 2, Continue=True)
        # -----//-----k组测量，依次进行一次，然后整体重复N次。每次测量包含 cooling 和 detect 两次计数
        # 打印结果，重复运行200次
        th = ThresIntegrator(k * N * 2)  # 阈值判断对象
        Task.StartTask()
        fidelity = []
        for j in range(10):
            temp = th.with_cooling(Task.Read())  # 读取PMT数据,阈值判断，并返回积分数据
            proba = np.average(temp.reshape(N, k), axis=0)  # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
            print j + 1, proba
            detectedstate = get_rou_tomo(proba)
            tra = np.trace(np.dot(detectedstate, detectedstate))  # 求rou平方的迹

            mixpurerou = np.array(mixpurerou).astype(np.complex)
            fidelity0 = fidelity_func(detectedstate, mixpurerou)
            fidelity.append(fidelity0)  # 记录
        return np.array(detectedstate).astype(np.complex), fidelity, proba

        '''
    #-------输入一个任意态，与理论值对比--------
    trystate=mixpurestate[2]

    tryangle=pre_angle(trystate)
    #print tryangle
    mixpurerou=np.array(trystate.dot(trystate.conj().T))
    prepareList=[(freq1*2*pi,tryangle[2]*RabiTime1/2.0,pi*tryangle[3]),
    (freq2*2*pi,tryangle[0]*RabiTime2/2.0,pi*tryangle[1]),(freq1*2*pi,tryangle[4]*RabiTime1/2.0,pi*tryangle[5])] # Note the 2*pi
    detectedstate,fidelity,proba=test_fidelity(prepareList,mixpurerou)
    print'----------trystate----------'
    print mp.matrix(trystate)
    print'--------detectedstate--------'
    print mp.matrix(detectedstate)
    print'---------mixpurerou---------'
    print mp.matrix(mixpurerou)
    print'---------fidelity----------'
    print mp.matrix(fidelity)
    print'---------tryangle----------'
    print mp.matrix(tryangle)  
    print'---------pre_angle_list----------'
    print mp.matrix(pre_angle_list)
    
    import matplotlib.pyplot as plt
    plt.plot(np.real(fidelity))
    plt.show()
    #*****************************************
    
    '''
    # def NCHV(rho,proba):

    for i in range(state_numb):
        # 注意，【信号源角频率，操作时间（不是旋转角度！），相位】
        print pre_angle_list
        prepareList = [(freq1 * 2 * pi, pre_angle_list[i][2] * RabiTime1 / 2.0, pi * pre_angle_list[i][3]),
                       (freq2 * 2 * pi, pre_angle_list[i][0] * RabiTime2 / 2.0, pi * pre_angle_list[i][1]),
                       (freq1 * 2 * pi, pre_angle_list[i][4] * RabiTime1 / 2.0, pi * pre_angle_list[i][5])]
        detectedstate, fidelity, proba = test_fidelity(prepareList, mixpurerou[i])
        proba_list.append(proba)
        rho_threepart.append(detectedstate)
        print '----------mixpurerou%d--------------' % (i + 1)
        print mp.matrix(mixpurerou[i])
        print '----------detectedstate%d--------------' % (i + 1)
        print mp.matrix(detectedstate)
        print '----------fidelity%d--------------' % (i + 1)
        print mp.matrix(fidelity)
        print '\n'
        # print '----------mixpurestate%d--------------'%(i+1)
        # print mp.matrix(mixpurestate[i])
        fidelity_list.append(fidelity)
        rho_tomo += rho_threepart[i] * mixprobability[i]  # 将三个分量概率相加得到最终的密度矩阵

    tra = np.trace(np.dot(rho_tomo, rho_tomo))  # 求rou平方的迹
    Fidelity = fidelity_func(rho_tomo, rouxyz)
    print 'final result'
    print '----------rou_prepare--------------'
    print mp.matrix(rouxyz)
    print '----------rho_tomo--------------'
    print mp.matrix(rho_tomo)
    print '----------Fidelity,tra--------------'
    print 'Fidelity=%0.3f    trace=%0.3f' % (np.real(Fidelity), np.real(tra))
    # print fidelity_list
    '''
    # output all result to a pkl file for future use      
    data=(rouxyz,state_numb,mixpurestate,mixpurerou,mixprobability,rho_threepart,fidelity_list,rho_tomo,Fidelity,tra,N)
    wignerf=r'D:stabilizer_tomo_result_of_x-%fy-%f.pkl'%(x_point,y_point)
    with open(wignerf, 'w') as f:
        p.dump(data,f)
    '''

    plt.plot(np.real(fidelity_list[0]))
    plt.plot(np.real(fidelity_list[1]))
    plt.plot(np.real(fidelity_list[2]))
    plt.show()
    triggerstop()
