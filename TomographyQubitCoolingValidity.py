# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2016

@author: jmcui
"""

from TomographyQubit import *
from pyqtgraph.Qt import QtGui, QtCore
from Config import wordCoolCount, wordCoolCountTriger
from setAFG import RabiTime1, RabiTime2, freq1, freq2

import ctypes
import datetime
import pickle

from numpy import *
# rewrite GenInstList function with cooling counting


def GenInstList(timeTag, CoolingTime, PumpingTime, DetectingTime):  # rewrite GenInstList function with cooling counting
    # Need to change python float to c double, utilize ctypes.c-double().value
    InstListData = []
    InstListData.append((wordCoolCount, OptC.CONTINUE, 0, 1000 * unit.ms)) 
    for ts, pl in timeTag:  # ts: time start, pl :pulse length
        if ts < PumpingTime:
            InstListData.append((wordCoolCount, OptC.CONTINUE, 0, ctypes.c_double(
                CoolingTime * unit.ms).value))  # cooling and coutering
            # InstListData.append((wordCooling,OptC.WAIT,0,ctypes.c_double(1.0*unit.us).value))
            # # wait for line trigger
            InstListData.append((wordPumping, OptC.CONTINUE, 0, ctypes.c_double((PumpingTime - ts) * unit.us).value))
            if ts == 0:
                InstListData.append((wordOperatingTriger, OptC.CONTINUE, 0, ctypes.c_double(pl * unit.us).value))
            else:
                InstListData.append((wordPumpingTriger, OptC.CONTINUE, 0, ctypes.c_double(ts * unit.us).value))
                InstListData.append((wordIdle, OptC.CONTINUE, 0, ctypes.c_double(1.0 * unit.us).value)) # compenstion AOM delay
                InstListData.append((wordOperating, OptC.CONTINUE, 0, ctypes.c_double(pl * unit.us).value))
        else:
            InstListData.append((wordCoolCount, OptC.CONTINUE, 0, ctypes.c_double(
                CoolingTime * unit.ms - (ts - PumpingTime) * unit.us).value))  # cooling and coutering
            # InstListData.append((wordCooling,OptC.WAIT,0,ctypes.c_double(1.0*unit.us).value))
            # # wait for line trigger
            InstListData.append((wordCoolCountTriger, OptC.CONTINUE, 0, ctypes.c_double(
                (ts - PumpingTime + 0.010) * unit.us).value))   # cooling and coutering
            InstListData.append((wordPumping, OptC.CONTINUE, 0, ctypes.c_double((PumpingTime - 0.010) * unit.us).value))
            InstListData.append((wordIdle, OptC.CONTINUE, 0, ctypes.c_double(1.0 * unit.us).value)) # compenstion AOM delay
            InstListData.append((wordOperating, OptC.CONTINUE, 0, ctypes.c_double(pl * unit.us).value))
        InstListData.append((wordDetecting, OptC.CONTINUE, 0, ctypes.c_double(DetectingTime * unit.us).value))
        InstListData.append((wordIdle, OptC.CONTINUE, 0, ctypes.c_double(1.0 * unit.us).value)) # # compenstion AOM delay
    InstListData.append((wordIdle, OptC.BRANCH, 1, ctypes.c_double(1.0 * unit.us).value))  # delay,then next cycle
    return InstListData

# define a thresholder with integrator


class ThresIntegrator:
    cycle = 0

    def __init__(self, size):
        self.cycle = 0  # used for the simple threshold judgement
        self.data = np.zeros(size)

        # using cooling counts to determine whether it's effective
        self.effective_cycles = np.zeros(size / 2, dtype=np.int32) + 1
        self.with_cooling_data = np.zeros(size / 2)

    def simple(self, data, threshold=1):  # simple threhold method, dark state probability
        # dark state is 1, bright state is 0
        self.data += (data < threshold) * 1.0
        self.cycle += 1
        self.probability_simple = self.data / self.cycle
        return self.probability_simple  # return probability of dark state

    # using cooling counts to check the detection validity
    def with_cooling(self, data, cooling_threshold=20.0, threshold=1.5):
        # seperate cooling and detection counts
        cooling_counts, counts = data.reshape(len(data) / 2, 2).T
        # print cooling_counts
        # effective or not in cooling and detection cyles
        eff = cooling_counts > cooling_threshold
        self.eff_count = eff * 1
        # in effevtive cycle and dark state, record
        self.with_cooling_data += (counts < threshold) * 1.0 * eff
        self.probability_with_cooling = self.with_cooling_data / \
            self.effective_cycles  # calculate and save to the probability array
        # increase the array numbers for the effective cycles
        self.effective_cycles += eff * 1
        self.cycle += 1
        # return probability of dark state with cooling check
        return self.probability_with_cooling


class ThresDataPool:
    def __init__(self, filename='temp.dat', fileinfo=None):
        self.cycle = 0
        self.DataPool = []
        self.file = open(filename, 'a')
        cnt_filename=filename+'.cnt'
        self.f_cnt = open(cnt_filename, 'a')
        if fileinfo is None:
            now = datetime.datetime.now()
            timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
            fileinfo = '\n# ' + timestamp + '\n\n'
        else:
            fileinfo = '\n# '+ fileinfo+ '\n\n'
        self.file.write(fileinfo)
        self.f_cnt.write(fileinfo)

    # using cooling counts to check the detection validity
    def with_cooling(self, data, cooling_threshold=20.0, threshold=1.5,save_effective=True):
        # seperate cooling and detection counts
        cooling_counts, counts = data.reshape(len(data) / 2, 2).T
        error = cooling_counts < cooling_threshold  # error cyles
        # error cylces value is negtive, -1 or -2
        if save_effective:
            if np.sum(1*error)>0:
                return None
            else:
                self.cycle += 1
                data=(counts<threshold)*1
                dcnt=counts
                self.DataPool.append(data.tolist())
                self.file.write(", ".join(repr(e) for e in data.tolist()))
                self.f_cnt.write("\t".join(str(e) for e in dcnt.tolist()))
                self.file.write("\n")
                self.f_cnt.write("\n")
        else:
            data = (counts < threshold) * 1 - 2 * error
            dcnt = counts*(1-2* error)-1* error
            self.cycle += 1
            self.DataPool.append(data.tolist())
            self.file.write(", ".join(repr(e) for e in data.tolist()))
            self.f_cnt.write("\t".join(str(e) for e in dcnt.tolist()))
            self.file.write("\n")
            self.f_cnt.write("\n")
        # print(", ".join(repr(e) for e in data.tolist()),, file=self.file)
        # lines = numpy.loadtxt("temp_data",dtype=numpy.int8)
        self.file.flush()
        self.f_cnt.flush()
        return data

    def save(self, filename):
        f = open(filename, 'w')
        pickle.dump(self.DataPool, f, 0)
        f.close()
        self.f_cnt.close()

    def __del__(self):
        self.file.close()


class ThresDataPoolBinary(ThresDataPool):
    def __init__(self, groupsize=24, filename='temp_data', fileinfo=None):
        super(ThresDataPoolBinary, self).__init__(filename, fileinfo)
        self.groupsize = groupsize
        self.effective_cycle = 0

    def with_cooling(self, data, cooling_threshold=20.0, threshold=1.0):
        if len(data) % (2 * self.groupsize):
            print("data size %d is not divisable by 2*groupsize %d\n " %
                  (len(data), 2 * self.groupsize))
        cooling_counts, counts = data.reshape(len(data) / 2, 2).T
        cooling_counts = cooling_counts.reshape(
            len(cooling_counts) / self.groupsize, self.groupsize)
        # this function is under develop


if __name__ == "__main__":  # -----//------流程图详解见下
    global DG5252
    # set pamameters
    # RabiTime=us,注意，这个时间是将0变成-0态的时间，真正的一个拉比周期是T=RabiTime×2

    RabiTime = RabiTime2
    Freq = freq2

    CoolingTime = 4  # ms
    PumpingTime = 400  # us
    DetectingTime = 1000  # us#-----//------设置冷却，泵浦，探测时间，可以根据实际需求进行优化
    print
    '''
    # pulsePrepare=0.0
    # pulsePrepare=RabiTime/2.0
    piTime=RabiTime/2.0  #-----//------用Rabi.py文件或者tomography_Rabi.py程序测量pi脉冲时间，从0翻转到1
    prepare=(RabiTime*1.25, 0.0)#-----//------设置制备态的脉冲（时长，相位）,根据微波旋转矩阵可以求出

    a=2*pi*prepare[0]/RabiTime
    phiadjust=0
    phi1=prepare[1]+phiadjust
    # 这个旋转矩阵是根据段路明组那篇文章中的矩阵来的。
    R1=np.array([[cos(a/2.0),(cos(phi1)+I*sin(phi1))*sin(a/2.0)],[-(cos(phi1)-I*sin(phi1))*sin(a/2.0),cos(a/2.0)]])#-----//------微波旋转矩阵
    initialstate=np.array([[1],[0]])#-----//------初始为0态
    preparestate=np.dot(R1,initialstate)
    rou_prepare=np.dot(preparestate,np.conj(preparestate.T))#-----//------理论上想要制备的态的密度矩阵
    #-----//------下面三个是测量需要用到的旋转操作，将目标态转到0态上，以便进行测量
    zrot=(piTime,0)#Z(pi,0)
    xrot=(piTime/2.0,pi) #X(pi/2.0,0)通过控制 operation 区分 x操作和 z翻转
    yrot=(piTime/2.0,-pi/2.0)#Y(pi/2.0,-pi/2.0)

    pulselist=(prepare,zrot,prepare,yrot)#-----//------生成操作序列，每个元素为一个操作，制备或者旋转，每个操作格式都是（时长，相位），三维的还要有频率

    '''
    # pulsePrepare=0.0
    piTime = RabiTime / 2.0  # pi flip time
    # 根据下面的微波旋转矩阵，有如下的关系
    # Z[0,1]态，生成R（pi,-pi/2.0），测量R（pi，pi/2.0）测量的phi为生成矩阵中的phi+pi
    # X[1,1]态，生成R（pi/2.0,-pi/2.0）,测量R（pi/2.0,pi/2.0）为Z旋转的一半
    # Y[1,i]态，生成R（pi/2.0,pi）,测量R（pi/2.0,0）
    phiadjust = 0
    pulsePrepare = piTime / 2.4
    prepare = (pulsePrepare, pi / 1.5)  # according to initial state

    # xia
    zrot = (piTime, pi / 2.0)
    xrot = (piTime / 2.0, pi / 2.0)  # 通过控制 operation时间 区分 x操作和 z翻转
    yrot = (piTime / 2.0, 0)

    prepare = (piTime * 1.23, -pi / 0.29)

    # 生成初始密度矩阵rou_prepare
    a = 2 * pi * prepare[0] / RabiTime
    phi1 = prepare[1]

    # 这个旋转矩阵是根据Monroe组2007年Daniel Lynn Stick的博士论文16页式2.16来的。
    # R1=np.array([[cos(a/2.0),-I*(cos(phi1)+I*sin(phi1))*sin(a/2.0)],[-I*(cos(phi1)-I*sin(phi1))*sin(a/2.0),cos(a/2.0)]])

    # 下面的矩阵是我们qutrit所用矩阵的二维形式，角度由上面的角度变换而来
    R1 = np.array([[cos(a / 2.0), I * (cos(phi1) + I * sin(phi1)) * sin(a / 2.0)],
                   [I * (cos(phi1) - I * sin(phi1)) * sin(a / 2.0), cos(a / 2.0)]])

    initialstate = np.array([[1], [0]])
    preparestate = np.dot(R1, initialstate)
    print(preparestate)
    rou_prepare = np.dot(preparestate, np.conj(preparestate.T))
    # rou_prepare1=np.array([[0,0],[0,1]])/1.0#理论态，计算机有误差
    # fidelity1=np.trace(np.dot(rou_prepare,rou_prepare1))
    # print rou_prepare
    print(zrot[1])

    pulselist = (prepare, (zrot[0], zrot[1] + pi),
                 prepare, (yrot[0], pi + yrot[1]))
    # pulselist=(prepare,yrot)
    # -----//------利用Genwave函数生成函数值，详细说明见函数内
    WAV, timeTag = GenWave(pulselist, Freq * 2 * pi)
    # -----//------将函数值写入任意波发生器，脉冲幅度在这里设置
    setChannelWave32K(2, WAV, amplitude=1.0)
    # print timeTag
    # plt.plot(WAV)
    # plt.show()
    #-----//------写入函数值之后，最重要的就是确定Gate的时刻和时长，用timeTag来标定每个循环时在何时开放微波，何时关闭。只有微波开启，才会和任意波发生器的微波混合，否则相当于不打微波
    t0, dt = zip(*timeTag)  # -----//------t0时刻列表，dt 脉冲长度列表
    timeTag = ((t0[0], dt[0]), (t0[0], dt[0] + dt[1]), (t0[0], dt[0] + dt[1] / 2.0),
               (t0[2], dt[2] + dt[3]))  # -----//-----z0, z1, x0, y0态制备与旋转脉冲的时间标定

    k = len(timeTag)  # 输入脉冲数
    # timeTag=((start0,pulse length1),(stardt,pulse length2),(start3,pulse
    # length3)...) 对应Z0,Z1,X0,Y0四组测量
    t0, dt = zip(*timeTag)
    # -----//-----compesation delay of AFG, we should triger the AFG in advance
    t0 = np.array(t0) + 1.03
    timeTag = zip(t0, dt)  # -----//-----给t0加上triger触发修正之后重新得到timeTag

    print(timeTag, t0[0])
    N = 100
    # -----//-----将脉冲时序传给控制卡
    InstListData = GenInstList(
        timeTag, CoolingTime, PumpingTime, DetectingTime)
    # -----//-----k组测量，依次进行一次，然后整体重复N次。每次测量包含 cooling 和 detect 两次计数
    Task = SpinTask(InstListData, k * N * 2, Continue=True)
    print(InstListData)
    Task.StartTask()

    # 打印结果，重复运行200次
    th = ThresIntegrator(k * N * 2)  # 阈值判断对象
    import random
    import scipy.linalg as linalg

    def fidelity_func(detec_rho, rho):
        a = linalg.sqrtm(detec_rho)
        fidelity = np.trace(linalg.sqrtm(np.dot(np.dot(a, rho), a)))
        return fidelity
    print('------iujj')
    bbb = preparestate.T[0]
    atate_list = [np.array([1, 0]), np.array([0, 1]), np.array(
        [1, 1]) / sqrt(2), np.array([1, 1.0j]) / sqrt(2)]
    theory_prob = [np.abs(bbb.dot(aa.conjugate()))**2 for aa in atate_list]
    print(theory_prob)
    for j in range(40):
        # 读取PMT数据,阈值判断，并返回积分数据,测试+[random.random()*5 for i in range(k*N*2)
        temp = th.with_cooling(Task.Read())
        # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
        proba = np.average(temp.reshape(N, k), axis=0)
        # print th.effective_cycles
        print(proba)
        print(theory_prob)
        print(proba - theory_prob)
        print('\n')
        # time.sleep(0.05)

        detectestate = (unitymatrix + sigmax * (2 * proba[2] - 1) + sigmay * (
            2 * proba[3] - 1) + sigmaz * (proba[0] - proba[1])) / 2.0
        tra = np.trace(np.dot(detectestate, detectestate))  # 求rou平方的迹
        fidelity0 = fidelity_func(detectestate, rou_prepare)

        # fidelity1=np.trace(np.sqrt(np.dot(np.dot(np.sqrt(preparestate),detectestate),np.sqrt(preparestate))))**2
        print(np.real(fidelity0), np.real(tra))
        print(detectestate)
        print(rou_prepare)

    triggerstop()  # 停止放大器工作，继续cooling
