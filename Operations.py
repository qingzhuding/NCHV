# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 14:15:43 2016

@author: jmcui
"""

import numpy as np
from numpy import pi
import pickle
import datetime
import ctypes
import os
import sys

SITE_ROOT = os.path.dirname(os.path.realpath(__file__))
PARENT_ROOT = os.path.abspath(os.path.join(SITE_ROOT, os.pardir))
sys.path.append(PARENT_ROOT)

import spinpy.OperationCode as OptC
import spinpy.unit as unit
from Config import wordCooling, wordPumping, wordOperating, wordDetecting
from Config import wordIdle, wordPumpingTriger, wordCoolingTriger, wordOperatingTriger
from Config import wordCoolCount, wordCoolCountTriger, wordIdleTriger


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


def GenAFG_TagInst(timeTag, CoolingTime, PumpingTime,
                   DetectingTime, FirstCooling=1000.0 * unit.ms,
                   WaitLineTriger=False, CountInCooling=True):
    ''' GenInstList with line triger

    unit time in timeTag, CoolingTime, PumpingTime, DetectingTime, FirstCooling
    are same with the unit of spincore config
    '''
    # make sure the time data format is float, not long double
    # it can be avoid some problems caused by sympy data type
    CoolingTime = ctypes.c_double(CoolingTime).value
    PumpingTime = ctypes.c_double(PumpingTime).value
    DetectingTime = ctypes.c_double(DetectingTime).value
    FirstCooling = ctypes.c_double(FirstCooling).value
    timeTag = [(ctypes.c_double(ts).value, ctypes.c_double(pl).value)
               for ts, pl in timeTag]

    InstListData = []
    #  first long cooling to make sure start sync with counter
    InstListData.append((wordCooling, OptC.CONTINUE, 0, FirstCooling))
    # whether to wait for the line trigger
    if WaitLineTriger:
        InstListData.append((wordCooling, OptC.WAIT, 0, 1.0 * unit.us))

    if CountInCooling:
        wCool = wordCoolCount
        wCoolTriger = wordCoolCountTriger
    else:
        wCool = wordCooling
        wCoolTriger = wordCoolingTriger

    for ts, pl in timeTag:  # ts: time start, pl :pulse length
        if ts < PumpingTime:
            InstListData.append((wCool, OptC.CONTINUE, 0, CoolingTime))
            if (PumpingTime - ts) > 10 * unit.ns:
                InstListData.append((wordPumping, OptC.CONTINUE, 0, PumpingTime - ts))
            if ts < 10 * unit.ns:
                InstListData.append((wordIdleTriger, OptC.CONTINUE, 0, 1.0 * unit.us))
            else:
                InstListData.append((wordPumpingTriger, OptC.CONTINUE, 0, ts))
                InstListData.append((wordIdle, OptC.CONTINUE, 0, 1.0 * unit.us))
        else:
            tx = CoolingTime + PumpingTime - ts
            if tx < 10 * unit.ns:
                raise ValueError('CoolingTime + PumpingTime is less than lists range time')
            InstListData.append((wCool, OptC.CONTINUE, 0, tx))
            tx = ts - PumpingTime + 10 * unit.ns  # make sure > 10 ns
            InstListData.append((wCoolTriger, OptC.CONTINUE, 0, tx))
            tx = PumpingTime - 10 * unit.ns
            InstListData.append((wordPumping, OptC.CONTINUE, 0, tx))
            InstListData.append((wordIdle, OptC.CONTINUE, 0, 1.0 * unit.us))
        InstListData.append((wordOperating, OptC.CONTINUE, 0, pl))
        InstListData.append((wordDetecting, OptC.CONTINUE, 0, DetectingTime))
        InstListData.append((wordIdle, OptC.CONTINUE, 0, 1.0 * unit.us))
    # delay,then next cycle
    InstListData.append((wordIdle, OptC.BRANCH, 1, 1.0 * unit.us))
    return InstListData


class ThresIntegrator:
    cycle = 0

    def __init__(self, size):
        self.cycle = 0  # used for the simple threshold judgement
        self.data = np.zeros(size)

        # using cooling counts to determine whether it's effective
        self.effective_cycles = np.zeros(size / 2, dtype=np.int32) + 1
        self.with_cooling_data = np.zeros(size / 2)

    def simple(self, data, threshold=1):
        # simple threhold method, dark state probability
        # dark state is 1, bright state is 0
        self.data += (data < threshold) * 1.0
        self.cycle += 1
        self.probability_simple = self.data / self.cycle
        return self.probability_simple  # return probability of dark state

    # using cooling counts to check the detection validity
    def with_cooling(self, data, cooling_threshold=20.0, threshold=1.0):
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
        cnt_filename = filename + '.cnt'
        self.f_cnt = open(cnt_filename, 'a')
        if fileinfo is None:
            now = datetime.datetime.now()
            timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
            fileinfo = '\n# ' + timestamp + '\n\n'
        else:
            fileinfo = '\n# ' + fileinfo + '\n\n'
        self.file.write(fileinfo)
        self.f_cnt.write(fileinfo)

    # using cooling counts to check the detection validity
    def with_cooling(self, data, cooling_threshold=20.0, threshold=1.5, save_effective=True):
        # seperate cooling and detection counts
        cooling_counts, counts = data.reshape(len(data) / 2, 2).T
        error = cooling_counts < cooling_threshold  # error cyles
        # error cylces value is negtive, -1 or -2
        if save_effective:
            if np.sum(1 * error) > 0:
                return None
            else:
                self.cycle += 1
                data = (counts < threshold) * 1
                dcnt = counts
                self.DataPool.append(data.tolist())
                self.file.write(", ".join(repr(e) for e in data.tolist()))
                self.f_cnt.write("\t".join(str(e) for e in dcnt.tolist()))
                self.file.write("\n")
                self.f_cnt.write("\n")
        else:
            data = (counts < threshold) * 1 - 2 * error
            dcnt = counts * (1 - 2 * error) - 1 * error
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

    def __del__(self):
        self.file.close()
        self.f_cnt.close()


if __name__ == "__main__":
    pass
