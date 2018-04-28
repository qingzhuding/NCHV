# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:15:43 2016

@author: jmcui
"""
from TomographyQubitCoolingValidity import *
import numpy as np
from numpy import pi
from itertools import chain
import scipy.linalg as linalg
from optimal_tomo_function import *
from RamseyProbe import CalibrateRamsey
import datetime


def GenAFG_Waves(pulse_list, SamplingRate=1000):
    '''由波形参数表，创建AFG的波形

    # time unit is us
    # 波形参数表：pulselist
     [(omega,pulselen, phai),(omega,pulselen, phai),(omega,pulselen, phai),...]
    # 返回:波形 和 波形分段位置所标记的timetags,timetag的形式是（该脉冲起始时刻，该脉冲长度）
    '''
    omegaArray, pulse_len, phase = zip(*pulse_list)
    TotalTime = np.sum(pulse_len)
    TotalLen = int(TotalTime * SamplingRate)
    # rescale time
    TotalTime = TotalLen * 1.0 / SamplingRate
    # 给每个点一个时间点 , 0 is included in the time list
    t = np.linspace(0, TotalTime, TotalLen + 1)
    # 将要给每个点一个相位，与所在脉冲有关，在此初始为0
    fai = np.zeros(TotalLen + 1)
    # +omegaArray[0]#给每个点一个频率Hz
    omega = np.zeros(TotalLen + 1)
    # print omegaArray
    # 切片起始起索引，第一个切片的起始索引地址应该为0
    t0_index = 0
    time_tags = np.zeros(len(pulse_len) + 1)
    # 定义n个脉冲的n+1个端点，作为标记点（相变点）
    # 寻找每个脉冲对应切片，并设置相应的相位和频率
    for i in range(len(pulse_len)):
        time_tags[i + 1] = np.sum(pulse_len[0:i + 1])
        # time_tags[i+1]表示第i+1个脉冲结束时刻，等于前i+1个脉冲总时长
        # i最大为N-1，i+1 最大为N
        # 上句话标记了所有脉冲终点时刻的绝对时间，起点为0
        t_index = int(time_tags[i + 1] * SamplingRate)
        # 计算切片重点相对于0点的索引时刻，time_tags[i+1]表示第i+1个脉冲结束时刻，
        # 或者说从0时刻到i脉冲结束时刻的时间
        fai[t0_index:t_index] = phase[i]
        # 切片，并设置这一段的相位
        omega[t0_index:t_index] = omegaArray[i]
        # 设置这一段的角频率
        t0_index = t_index
        # 下一个相位区的起始索引

    return np.cos(omega * t - fai), zip(time_tags[0:len(pulse_len)], pulse_len)
    # 这里给出了所有脉冲的timetag（起始时间，脉冲时长），
    # 然而实际用到的还需要CombinePulseLists中的timetag来分割


def OperationsToPulses(olist, l_para):
    '''由旋转角度生成波形参数表

    # 输入：旋转角度 olist：
    #              [(theta1，fai1, level_index1), (theta2，fai2, level_index2)..]
    #      能级参数 l_para: [(RabiTime1,freq1),(RabiTime2,freq2)...]
    #      level_index 是一个整数，用来选择l_para的参数
    # 输出：波形参数表 ［(omega1,pulselen1, phai1),(omega2,pulselen2, phai2)...］
    '''
    pulse_list = []
    for o in olist:
        if o is not None:
            theta, phai, i = o
            RabiTime, freq = l_para[i]
            pulse = (2 * pi * freq, RabiTime * theta / (2 * pi), phai)
            pulse_list.append(pulse)
    return pulse_list


def CombineOperationLists(prepareLists, MeasureLists):
    ''' 合并两个操作列表, 输出操作－测量元组列表

    注意两个操作列表的大小需要相等
    输入：prepareLists，或者 MeasureLists 形式：
           [prepare1, prepare2, ...]
         其中 prepare1 ＝ [o01,o02,o03,..] 一系列操作元构成一组操作
         因此最终 preLists 的形式：操作元的二维数列
           [[o01,o02,o03,..],[o11,o12,o13,..],]
         其中操作元 o ＝ （theta，fai， i）
    输出：合并结果为 operation_list 形式：
          ［[prepare1，measure1]，[prepare2，measure2]...］
         其中每个列表元素表示一组操作（可以看作一组制备测量过程）
         最终operation_list的形式为：操作元的二维数列
           ［［o1,o2,o3,o4...］,...］

    '''
    n1 = len(prepareLists)
    n2 = len(MeasureLists)

    if n1 != n2:
        raise ValueError("prepare and measure list length is not equal!")
    lists = [prepareLists[i] + MeasureLists[i] for i in range(n1)]
    return lists


def OperationsToWave(olists, l_para):
    '''由操作列表生成任意波形函数表

    olists的形式为：olist 的列表
     ［［o1,o2,o3,o4...］,...］
    每个元组表示一组制备测量，是一系列operation组成的元组，
    可以看作制备测量操作元的组合［o1,o2,o3,o4..］＝ [prepare1，measure1]
    l_para 为能级越迁参数列表，
       [(RabiTime1,Freq1), (RabiTime2,Freq2),....]

    '''
    plists = [OperationsToPulses(olist, l_para) for olist in olists]
    # 将二维lists展成一维list
    plist = []
    for e in plists:
        plist.extend(e)
    # 生成TimeTagIndex 标记测量位置 Ti
    ii = [len(olist) for olist in olists]
    Ti = [0]
    a = 0
    for i in ii:
        a += i
        Ti.append(a)
    # Generate waves and TimeTags
    waves, Tags = GenAFG_Waves(plist)
    t, dt = zip(*Tags)
    N = range(len(ii))
    TimeTags = [(t[Ti[i]], sum(dt[Ti[i]:Ti[i + 1]])) for i in N]
    return waves, TimeTags


def pre_angle(vector):
    ''' get prepare rotation angles from a state vector

    作用顺序为R23*R13*[0,0,1]=vector
    engery level:
    [1]  F=1,m=0           |2>
                      |1>  /
    [2]  F=1,m=1       |  /
                       | /
    [3]  F=0,m=0      |3>

    Transitions:
    0   |1> - |3>
    1   !2> - |3>
    '''
    # 格式重整为横向量的形式
    if np.size(vector.dot(vector.conjugate().T)) != 1:
        vector = np.array([vector[i][0] for i in range(len(vector))])
    # 归一化
    vector = vector / np.sqrt(vector.dot(vector.conjugate().T))
    # 将vector第三分量化为实数，即去除vector的全局相位，
    fai = np.arctan2(np.imag(vector[2]), np.real(vector[2]))
    vector = vector * np.exp(-1.0j * fai)
    # 如果矢量为[1，0，0]，直接给出结果，避免奇异
    if np.abs(np.abs(vector[0]) - 1.0) < 1e-14:
        [a, phi1, b, phi2] = [0, 0, pi, pi / 2.0]
    else:
        b = np.abs(vector[1]) / np.sqrt(1 - np.abs(vector[0]) ** 2)
        if (b - 1.0) > 1e-17:  # 防止float最后一位错误
            b = 1.0
        a = 2 * np.arcsin(b)
        phi1 = pi / 2.0 - np.angle(vector[1])
        b = 2 * np.arcsin(np.abs(vector[0]))
        phi2 = pi / 2.0 - np.angle(vector[0])
    # 制备时，先打R13(b,phi2,0),再打R23(a,phi1,1)到基态
    return [(b, phi2, 0), (a, phi1, 1)]


def m_angle(vector):
    ''' get a measurement rotation angles from a state vector
    '''
    # 测量时，先打R23(a,phi1',1),再打R13(b,phi2',0)到待测基。a,b与制备时相同，相位phi'=phi+pi
    R1, R2 = pre_angle(vector)
    theta1, fai1, i1 = R1
    theta2, fai2, i2 = R2
    return [(theta2, fai2 + pi, i2), (theta1, fai1 + pi, i1)]


def GenStablizerVectors(X):
    ''' 生成X的本征向量

    按照 1, w， w＊w 本征值排列本征向量
    '''
    w, v = np.linalg.eig(X)
    a = np.exp(-2.0j * pi / 3)
    index = [i for iw in [w, w * a, w * a * a]
             for i, e in enumerate(iw) if np.abs(e - 1) < 1e-15]
    vs = [np.array(v[:, i].T)[0] for i in index]
    return vs


def xy_get_rou(a, b):
    '''given phase point(x, y), return density matrix rou
    and plot the point in phase space
    '''
    w = np.exp(2.0j * pi / 3)
    rou = np.matrix([[-b + 4 / 9.0, 0, w**2 * (b - 1 / 9.0)],
                     [0, b + 2 / 9.0, (w - 1) * a - b + w ** 2 / 9.0 + 1 / 3.0],
                     [w * (b - 1 / 9.0), (w**2 - 1) * a - b + w / 9.0 + 1 / 3.0, 1 / 3.0]])
    return rou


def wignerbasis_rou_tomo(detectedprob):
    '''reconstruct density matrix with wigner basis result
    given detected probabilitY LIST of 12 unbiased basis of X,Z,XZ,XZZ
    '''
    v1 = np.sqrt(3) / 6.0
    m1, m2, m3, a1, a2, a3, m6, m4, m5, m8, m9, m7 = detectedprob[0:12]
    b1 = 0.5 * (m1 - m4 - m5 + m8)
    b2 = v1 * (-m1 - 2 * m2 - m4 + m5 + 2 * m7 + m8)
    c1 = 0.5 * (m1 - m7 + m5 - m8)
    c2 = v1 * (m1 + 2 * m2 - 2 * m4 - m5 + m7 - m8)
    d1 = 0.5 * (m1 + m4 + m7) - 0.5
    d2 = v1 * (-m1 - 2 * m2 - m4 - 2 * m5 - m7 - 2 * m8 + 3)
    rou_tomo = np.array([[a1, b1 - 1.0j * b2, c1 - 1.0j * c2], [b1 + 1.0j * b2, a2, d1 - 1.0j * d2],
                         [c1 + 1.0j * c2, d1 + 1.0j * d2, 1 - a1 - a2]])
    return rou_tomo


def Fidelity(rou_tomo, rou):
    '''given rou_tomo and rou, return fidelity
    '''
    def Matrix_sqrtm(M):
        '''get the root of an hermite matrix M
        1, 求本征值，本征态。
        2，本征值求开方
        3，再与原本征态投影算符相乘得到
        '''
        w, v = np.linalg.eig(M)
        M = (np.diag(np.sqrt(w)).dot(v.T)).T.dot(v.conj().T)
        return M

    a = Matrix_sqrtm(rou_tomo)
    b = (a.dot(rou)).dot(a)
    fidelity = np.trace(Matrix_sqrtm(b))
    return fidelity


def rhoTomogra(x_point, y_poin, mlists):
    '''  

    mlists 测量基序列，由多个测量基构成，
        每个测量基由一系列微波操作构成
    '''

    status = 'x_point=%s, y_point=%s' % (x_point, y_point)

    rouxy = xy_get_rou(x_point, y_point)
    w, v11 = np.linalg.eig(rouxy)
    w = np.real(w)
    # rourecord=[]
    # rourecord_optima=[]

    for k in [0, 1, 2]:
        # CalibrateRamsey first
        freq23 = CalibrateRamsey() * 1e-6
        # 配置能级参数
        l_para = [(RabiTime13, freq13), (RabiTime23, freq23)]
        state = np.array(v11[:, k])

        print 'state%s ' % state
        rou = state.dot(state.conj().T)

        # 定义初态 并生成初态制备列表 序列
        so = pre_angle(state)
        plists = [so] * len(mlists)

        prob = np.array([np.abs(v.conj().dot(state))**2 for v in VecList])
        print('theory prob:')
        print(prob)

        olists = CombineOperationLists(plists, mlists)
        waves, TimeTags = OperationsToWave(olists, l_para)
        t1, t2 = zip(*TimeTags)
        t1 = np.array(t1) + 0.03  # compesation delay of AFG, we should triger the AFG in advance
        TimeTags = zip(t1, t2)

        filename = '2017state-%.3f-%.3f' % (x_point, y_point)

        now = datetime.datetime.now()
        timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
        fileinfo = timestamp + ', k=%d' % (k)
        th = ThresDataPool(filename=filename, fileinfo=fileinfo)  # 阈值判断对象
        Integrat = ThresIntegrator(len(Slist) * 3 * 2)  # 积分记录，跟踪保真度
        AFG_cycle = setChannelWave32K(2, waves, amplitude=1.0)  # set AFG ! 注意振幅
        print("AFG cycle was set as %.3f us" % (AFG_cycle / 1000.0))

        InstList = GenInstList(TimeTags, CoolingTime, PumpingTime, DetectingTime)
        Task = SpinTask(InstList, len(olists) * 2, Continue=True)  # set hardware !
        Task.StartTask()

        for i in range(100):
            th.cycle = 0
            while th.cycle < 100:
                data = Task.Read()
                measures = th.with_cooling(data, cooling_threshold=3.0, threshold=1.5)
                temp = Integrat.with_cooling(data, cooling_threshold=3.0, threshold=1.5)  # 返回积分数据
                if np.sum(Integrat.eff_count) < 10:
                    import winsound
                    winsound.Beep(1200, 30)
            rou_tomo = wignerbasis_rou_tomo(temp)  # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
            optimal_rou_tomo = optimal_rho(temp)
            fidelity = Fidelity(rou_tomo, rou)
            optimal_fidelity = Fidelity(optimal_rou_tomo, rou)
            print [status, i + 1, np.sum(Integrat.eff_count), k]
            print temp.reshape(4, 3).T
            print np.sum(temp.reshape(4, 3).T, axis=0)
            print fidelity, optimal_fidelity
        Task.StopTask()
        Task.DAQ_counter.ClearTask()
        # rourecord.append(rou_tomo)
        # rourecord_optima.append(optimal_rou_tomo)
    # print 'total fidelity %s'%Fidelity(np.sum([w[i]*rourecord[i] for i in range(3)]),rouxy)
    # print 'total fidelity optimal %s'%Fidelity(np.sum([w[i]*rourecord_optima[i] for i in range(3)]),rouxy)


def check_point(x, y):
    a = x
    b = y
    dett = 3.0 * a**2 * b - 4 / 3.0 * a**2 + 3.0 * a * b**2 - 7 / 3.0 * a * \
        b + 4 / 9.0 * a - 4 / 3.0 * b**2 + 4 / 9.0 * b - 0.00823045267492819
    if dett < 1.0 + 1e16 and dett > -1e-16:
        return True
    else:
        return False

if __name__ == "__main__":
    from setAFG import *
    global DG5252

    # set pamameters
    '''
    engery level:
    [1]  F=1,m=0           |2>
                      |1>  /
    [2]  F=1,m=1       |  /
                       | /
    [3]  F=0,m=0      |3>

    Transitions:
    0   |1> - |3>
    1   !2> - |3>
    '''
    global freq23
    RabiTime13 = RabiTime1
    RabiTime23 = RabiTime2
    # us,注意，这个时间是将0变成-0态的时间，真正的一个拉比周期是T=RabiTime×2,
    # 用Rabi.py文件或者tomography_Rabi.py程序测量pi脉冲时间，从0翻转到
    freq13 = freq1  # MHz
    freq23 = freq2  # MHz

    # 设置冷却，泵浦，探测时间，可以根据实际需求进行优化
    CoolingTime = 1.5  # ms
    PumpingTime = 150  # us
    DetectingTime = 1000  # us
    #201709 remeasurement 
    task_list = [
        (0.1, -0.05), (0.1, 0.24), (0.1, -0.02), (0.1, 0.15), (0.1, 0.26), (0.1, 0.22),
        (0.1, 0.05), (0.1, 0.02), (0.1111111, 0.1111111), (0.1, 0.30), (0.1, 0.20), (0.1, 0.10),
        (0.1, 0.00), (0.2, 0.0), (0.2, 0.14), (0.2, 0.15), (0.2, 0.13), (0.2, -0.1),
        (0.2, 0.12), (0.2, 0.05), (0.2, -0.05), (0.2, 0.25), (0.2, 0.025), (0.2, -0.01),
        (0.2, 0.2), (0.2, -0.025), (0.2, -0.13), (0.2, 0.1), (0.3, 0.03), (0.3, 0.05),
        (0.3, -0.01), (0.3, -0.05), (0.3, 0.01), (0.3, 0.04), (0.3, -0.025), (0.3, 0.1),
        (0.3, 0.02), (0.05, 0.15), (0.05, 0.20), (0.15, 0.05), (-0.14, 0.25), (0.15, 0.13)]

    #,(0.05,0.25),(-0.075,0.25)
    # check the piont before run tasks
    for i, ta in enumerate(task_list):
        if not check_point(*ta):
            print(i, ta)
            raise ValueError('point out of physical state area.')

    # 生成stablelizer 基失列表
    # X，Z 矩阵，注意使用np.matrix
    X = np.matrix([[0, 0, 1.0], [1.0, 0, 0], [0, 1.0, 0]])
    Z = np.matrix([[1, 0, 0],
                   [0, np.exp(2.0j * pi / 3.0), 0],
                   [0, 0, np.exp(4.0j * pi / 3.0)]])
    Slist = [X, Z, X * Z, X * Z * Z]
    VecList = [GenStablizerVectors(s) for s in Slist]
    VecList = list(chain.from_iterable(VecList))
    # 生成24组测量的 mlists
    mlists = [m_angle(v) for v in VecList]

    np.set_printoptions(4)
    for x_point, y_point in task_list:
        print('---------------------------')
        rhoTomogra(x_point, y_point, mlists)

    print('Misson Complete!')
    print(datetime.datetime.now())

    import winsound
    while(1):
        winsound.Beep(1500, 500)
