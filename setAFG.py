# -*- coding: utf-8 -*-

import visa
import time
import numpy as np
from struct import pack

rm = visa.ResourceManager()

DG5252 = rm.open_resource(u'USB0::0x1AB1::0x0640::DG5T164850185::INSTR')
# DG5252=rm.open_resource(u'USB0::0x1AB1::0x0640::DG5T182600124::INSTR')
# 所有的文件都是从此处引用的参数

m = 0  # choose which frequency  m=0 freq1, m=1 freq2
Mode = 'SineWave'
# Mode = 'AFGSineWave' # AFG Wave

RabiTime1 = 20.21  # us

RabiTime2 = 18.67


freq1 = (200.0 * 1e6 + 700) * 1e-6   # MHz, m=0

freq2 = (210.228216e6 - 2300 - 2900 + 10000 - 847 - 7457) * 1e-6  # MHz,m=1


AFG_Amp = 1.0  # 生成任意波的振幅
AFG_N = 8388600  # 1G 采样率， 1000 000 点为1ms的持续时间
channel = 2
# freq2=(189.744927e6)*1e-6  #MHz,m=-1


# set AFG arbitrary wave

def setChannelWave32K(channel, values, amplitude=1.0):
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

    N = np.ceil(np.real(np.log(size * 1.0 / 32768) / np.log(2.0)))
    if N < -1:
        N = -1
    block_size = 32768 * (2**N)  # 32K的2^n 倍
    # print block_size,ceil(log(size*1.0/32768)/log(2.0))
    b = np.zeros(block_size - size, dtype=np.int16)  # 不足区域补零
    values = np.concatenate((values, b))  # lens to interger of 32K
    values = (values - values.min()) * 16383.0 / \
        (values.max() - values.min())  # Normlize to [0,16383]
    values = np.array(values, dtype=np.int16)  # to integter

    Ns = int(np.ceil(len(values) * 1.0 / 16384))
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


if __name__ == '__main__':

    channel = 2
    amplitude = 1.0
    Freqchoice = [freq1 * 1e6, freq2 * 1e6]

    Freq = Freqchoice[m]

    if Mode == 'SineWave':
        DG5252.write(':SOUR%d:BURSt OFF' % channel)
        DG5252.write(':OUTP%d:IMP 50' % channel)
        # 设置波形[<freq>[,<amp>[,<offset>[,<phase>]]]]
        DG5252.write(':SOUR%d:APPL:SIN %f,%f,0,0' % (channel, Freq, amplitude))
        DG5252.write(':SOUR%d:BURS:TRIG:SOUR EXT' % channel)
        DG5252.write(':SOUR%d:BURS:NCYC 1' % channel)
        DG5252.write(':OUTP%d ON' % channel)
        print("Set to %.6f MHz Sine wave" % (Freq * 1.0e-6))
    else:
        t = np.linspace(0, AFG_N * 1e-9, AFG_N)
        w = np.sin(2 * pi * Freq * t)
        setChannelWave32K(channel, w, AFG_Amp)
        print("Set to %.6f MHz AFG Sine wave" % (Freq * 1.0e-6))
        print("Note wave length is %f ms" % (t * 1000))
    print('Frequency m=1 is %.3f Hz' % (freq2 * 1e6))
