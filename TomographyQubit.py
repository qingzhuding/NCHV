# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2016

@author: jmcui,heran
"""
import numpy as np
import visa
import os
import sys
import time
import datetime
from struct import pack
import matplotlib.pyplot as plt
from math import *
from cmath import *
from scipy.linalg import expm2, eig, norm
from datetime import datetime
import matplotlib

SITE_ROOT = os.path.dirname(os.path.realpath(__file__))
PARENT_ROOT = os.path.abspath(os.path.join(SITE_ROOT, os.pardir))

sys.path.append(PARENT_ROOT)
import spinpy.OperationCode as OptC
import spinpy.unit as unit
from StateManipulation import SpinTask
from Config import wordCooling, wordPumping, wordOperating, wordDetecting, wordIdle, wordPumpingTriger, wordCoolingTriger, wordOperatingTriger
from setAFG import RabiTime1, RabiTime2, freq1, freq2
I = sqrt(-1.0)
# 泡利矩阵
unitymatrix = np.array([[1, 0], [0, 1]])
sigmax = np.array([[0, 1], [1, 0]])
sigmay = np.array([[0, -I], [I, 0]])
sigmaz = np.array([[1, 0], [0, -1]])


def setChannelWave32K(channel, values, amplitude=0.5):
    global DG5252
    rm = visa.ResourceManager()
    DG5252 = rm.open_resource(u'USB0::0x1AB1::0x0640::DG5T164850185::INSTR')
    # DG5252=rm.open_resource(u'USB0::0x1AB1::0x0640::DG5T182600124::INSTR')
    # 注意使用DG5252播放模式，此时采样率设置需要设置分频数，默认N=0，采样率为1GHz
    DG5252.write(':SOUR%d:BURSt OFF' % channel)
    # set impedance first,then set amplitude!
    DG5252.write(':OUTP%d:LOAD 50' % channel)
    DG5252.write(':SOUR%d:APPL:USER 1000,%f,0,0' %
                 (channel, amplitude))  # 这里频率1000实际上没有作用,[频率，幅度，偏置，相位]
    DG5252.write(':SOUR%d:FUNC:ARB:MODE PLAY' % channel)
    DG5252.write(':SOUR%d:BURS:MODE TRIG' % channel)
    DG5252.write(':SOUR%d:BURS:TRIG:SOUR EXT' % channel)
    DG5252.write(':SOUR%d:BURS:NCYC 1' % channel)
    time.sleep(0.03)
    header = ':SOURCE%d:TRAC:DATA:DAC16 VOLATILE,' % channel

    size = len(values)

    N = ceil(np.real(log(size * 1.0 / 32768) / log(2.0)))
    if N < -1:
        N = -1
    block_size = 32768 * (2**N)  # 32K的2^n 倍
    # print block_size,ceil(log(size*1.0/32768)/log(2.0))
    b = np.zeros(block_size - size, dtype=np.int16)  # 不足区域补零
    values = np.concatenate((values, b))  # lens to interger of 32K
    values = (values - values.min()) * 16383.0 / \
        (values.max() - values.min())  # Normlize to [0,16383]
    values = np.array(values, dtype=np.int16)  # to integter

    Ns = int(ceil(len(values) * 1.0 / 16384))
    for i in range(0, (Ns - 1) * 16384, 16384):
        subvals = values[i:i + 16384]
        waves = pack('<' + 'h' * len(subvals), *subvals)
        DG5252.write_raw(header + 'CON,#5' + '32768' + waves + '\n')
        time.sleep(0.03)

    subvals = values[(Ns - 1) * 16384:]
    waves = pack('<' + 'h' * len(subvals), *subvals)
    sEND = 'END,#532768'  # 532768
    # input the last data segment
    DG5252.write_raw(header + sEND + waves + '\n')
    time.sleep(0.03)
    DG5252.write(':OUTP%d ON' % channel)
    DG5252.write(':SOUR%d:BURSt ON' % channel)
    return block_size


def GenWave(pulse_list, omega=200 * 2 * pi, SamplingRate=1000):
    # time unit is us
    # pulselist form is [(pulselen, phai),(pulselen, phai),(pulselen, phai),......]
    # 返回波形 和 波形分段位置所标记的timetags,timetag的形式是（时刻，此时刻的脉冲长度）
    pulse_len, phase = zip(*pulse_list)
    TotalTime = np.sum(pulse_len)
    TotalLen = int(TotalTime * SamplingRate)
    TotalTime = TotalLen * 1.0 / SamplingRate  # rescale time
    # 0 is include in the time list
    t = np.linspace(0, TotalTime, TotalLen + 1)
    fai = np.zeros(TotalLen + 1)
    # print(t*1000) # check whether time line is correct
    t0_index = 0  # 切片起始起索引，第一个切片的起始索引地址应该为0
    time_tags = np.zeros(len(pulse_len) + 1)  # 定义n个脉冲的n+1个端点，作为标记点（相变点）
    for i in range(len(pulse_len)):  # 寻找每个脉冲对应切片，并设置相应的相位
        # 注意 [i+1] 寻址与[0:i+1]切片的区别，pulse_len的长度为N， i最大为N-1，i+1 最大为N， time_tags[N] 是其最后一个，而pulse_len[0:i+1]，切片0至N-1的N个数据
        time_tags[i + 1] = np.sum(pulse_len[0:i + 1])
        t_index = int(time_tags[i + 1] * SamplingRate)  # 计算切片结束索引
        fai[t0_index:t_index] = phase[i]  # 本次相位区切片，并设置这一段的相位
        t0_index = t_index  # 下一个相位区的起始索引
    return np.cos(omega * t - fai), zip(time_tags[0:len(pulse_len)], pulse_len)


def GenInstList(timeTag, CoolingTime, PumpingTime, DetectingTime):
    InstListData = []
    for ts, pl in timeTag:  # ts: time start, pl :pulse length
        if ts < PumpingTime:
            InstListData.append(
                (wordCooling, OptC.CONTINUE, 0, CoolingTime * unit.ms))
            # wait for line trigger
            InstListData.append((wordCooling, OptC.WAIT, 0, 1.0 * unit.us))
            InstListData.append(
                (wordPumping, OptC.CONTINUE, 0, (PumpingTime - ts) * unit.us))
            if ts == 0:
                InstListData.append(
                    (wordOperatingTriger, OptC.CONTINUE, 0, pl * unit.us))
            else:
                InstListData.append(
                    (wordPumpingTriger, OptC.CONTINUE, 0, ts * unit.us))
                InstListData.append(
                    (wordOperating, OptC.CONTINUE, 0, pl * unit.us))
        else:
            InstListData.append((wordCooling, OptC.CONTINUE, 0,
                                 CoolingTime * unit.ms - (ts - PumpingTime) * unit.us))
            # wait for line trigger
            InstListData.append((wordCooling, OptC.WAIT, 0, 1.0 * unit.us))
            InstListData.append(
                (wordCoolingTriger, OptC.CONTINUE, 0, (ts - PumpingTime + 0.010) * unit.us))
            InstListData.append(
                (wordPumping, OptC.CONTINUE, 0, (PumpingTime - 0.010) * unit.us))
            InstListData.append(
                (wordOperating, OptC.CONTINUE, 0, pl * unit.us))
        InstListData.append((wordDetecting, OptC.CONTINUE,
                             0, DetectingTime * unit.us))
        InstListData.append((wordIdle, OptC.CONTINUE, 0, 1.0 * unit.us))
    # delay,then next cycle
    InstListData.append((wordIdle, OptC.BRANCH, 0, 1.0 * unit.us))
    return InstListData


def triggerstop():  # 停止放大器工作，继续cooling
    ssss = [(wordCooling, OptC.BRANCH, 0, 4.0 * unit.ms)]
    SpinTask(ssss, 4, Continue=False)
    print('MW is stopped')
if __name__ == "__main__":  # -----//------流程图详解见下
    global DG5252
    # set pamameters
    RabiTime = 13.73  # us,注意，这个时间是将0变成-0态的时间，真正的一个拉比周期是T=RabiTime×2

    CoolingTime = 4  # ms
    PumpingTime = 400  # us
    DetectingTime = 600  # us#-----//------设置冷却，泵浦，探测时间，可以根据实际需求进行优化

    # pulsePrepare=0.0
    # pulsePrepare=RabiTime/2.0
    piTime = RabiTime / 2.0  # -----//------用Rabi.py文件或者tomography_Rabi.py程序测量pi脉冲时间，从0翻转到1
    # -----//------设置制备态的脉冲（时长，相位）,根据微波旋转矩阵可以求出
    prepare = (RabiTime / 2.5, pi / 1.4)

    a = 2 * pi * prepare[0] / RabiTime
    phiadjust = 0
    phi1 = prepare[1] + phiadjust
    # 这个旋转矩阵是根据段路明组那篇文章中的矩阵来的。
    R1 = np.array([[cos(a / 2.0), (cos(phi1) + I * sin(phi1)) * sin(a / 2.0)],
                   [-(cos(phi1) - I * sin(phi1)) * sin(a / 2.0), cos(a / 2.0)]])  # -----//------微波旋转矩阵
    initialstate = np.array([[1], [0]])  # -----//------初始为0态
    preparestate = np.dot(R1, initialstate)
    rou_prepare = np.dot(preparestate, np.conj(
        preparestate.T))  # -----//------理论上想要制备的态的密度矩阵
    #-----//------下面三个是测量需要用到的旋转操作，将目标态转到0态上，以便进行测量
    zrot = (piTime, 0)  # Z(pi,0)
    xrot = (piTime / 2.0, pi)  # X(pi/2.0,0)通过控制 operation 区分 x操作和 z翻转
    yrot = (piTime / 2.0, -pi / 2.0)  # Y(pi/2.0,-pi/2.0)

    # -----//------生成操作序列，每个元素为一个操作，制备或者旋转，每个操作格式都是（时长，相位），三维的还要有频率
    pulselist = (prepare, zrot, prepare, yrot)
    # -----//------利用Genwave函数生成函数值，详细说明见函数内
    WAV, timeTag = GenWave(pulselist, 200.0 * 2 * pi)
    # -----//------将函数值写入任意波发生器，脉冲幅度在这里设置
    setChannelWave32K(2, WAV, amplitude=0.5)
    print(timeTag)
    # plt.plot(np.real(WAV))
    # plt.show()
    #-----//------写入函数值之后，最重要的就是确定Gate的时刻和时长，用timeTag来标定每个循环时在何时开放微波，何时关闭。只有微波开启，才会和任意波发生器的微波混合，否则相当于不打微波
    t0, dt = zip(*timeTag)  # -----//------t0时刻列表，dt 脉冲长度列表
    timeTag = ((t0[0], dt[0]), (t0[0], dt[0] + dt[1]), (t0[0], dt[0] + dt[1] / 2.0),
               (t0[2], dt[2] + dt[3]))  # -----//-----z0, z1, x0, y0态制备与旋转脉冲的时间标定
    # timeTag=((start0,pulse length1),(stardt,pulse length2),(start3,pulse
    # length3)...) 对应Z0,Z1,X0,Y0四组测量
    t0, dt = zip(*timeTag)
    # -----//-----compesation delay of AFG, we should triger the AFG 1.290 us in advance
    t0 = np.array(t0) + 1.03
    timeTag = zip(t0, dt)  # -----//-----给t0加上triger触发修正之后重新得到timeTag

    # print timeTag,t0[0]

    N = 100
    # -----//-----将脉冲时序传给控制卡
    InstListData = GenInstList(
        timeTag, CoolingTime, PumpingTime, DetectingTime)
    print(InstListData)
    # -----//-----4组测量，依次进行一次，然后整体重复N次。可以减少因激光功率变化或者离子移动问题造成的误差
    Task = SpinTask(InstListData, 4 * N, Continue=True)
    Task.StartTask()

    # 打印结果，重复运行200次
    k = 4  # 输入脉冲数
    Data = np.zeros(k)  # -----//-----实际PMT采集到的光子数
    plus = np.zeros(k)  # -----//-----阈值判断后得到的比特值
    for j in range(15):
        temp = Task.Read()  # -----//-----读取PMT数据
        for i in range(N):
            for a in range(k):  # -----//-----对读数进行0，1赋值判断
                # print temp[k*i],temp[k*i+1],temp[k*i+2],temp[k*i+3]
                if temp[k * i + a] < 0.9:  # -----//-----阈值设定并赋值0和1
                    plus[a] = 1.0  # -----//-----测到暗态为1
                else:
                    plus[a] = 0.0
            # Data+=temp[k*i:k*i+4]
            Data += plus[0:k]
        proba = 0.01 * Data / (1.0 * (j + 1))  # -----//-----输出概率
        time.sleep(0.05)

        detectestate = (unitymatrix + sigmax * (2 * proba[2] - 1) + sigmay * (
            2 * proba[3] - 1) + sigmaz * (proba[0] - proba[1])) / 2.0
        tra = np.trace(np.dot(detectestate, detectestate))  # 求rou平方的迹
        fidelity0 = np.trace(np.dot(detectestate, rou_prepare))
        # fidelity1=np.trace(np.sqrt(np.dot(np.dot(np.sqrt(preparestate),detectestate),np.sqrt(preparestate))))**2
        print(np.real(fidelity0), np.real(tra))
        print(detectestate)
        print(rou_prepare)

    triggerstop()  # 停止放大器工作，继续cooling
