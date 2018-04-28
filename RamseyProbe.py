# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 14:15:43 2016

@author: jmcui
"""
from Operations import *
from StateManipulation import SpinTask
from scipy.optimize import curve_fit
from setAFG import *
import pandas as pd
import datetime


def GenRamseyInst(tlist, RabiTime,
                  CoolingTime=1.0 * unit.ms,
                  PumpingTime=100 * unit.us,
                  DetectingTime=1.0 * unit.ms,
                  FirstCooling=1000 * unit.ms):
    ''' Generate ramsey control instruction list

    '''
    Inst = []
    Inst.append((wordCooling, OptC.CONTINUE, 0, FirstCooling))
    for t in tlist:
        tms = int(t / unit.ms)  # ms
        tus = t - tms * unit.ms
        Inst.append((wordCooling, OptC.CONTINUE, 0, CoolingTime))
        Inst.append((wordPumping, OptC.CONTINUE, 0, PumpingTime))
        Inst.append([wordIdleTriger, OptC.CONTINUE, 0, 1.03 * unit.us])  # compensation AOM delay

        Inst.append((wordOperating, OptC.CONTINUE, 0, RabiTime / 4.0))  # pi/2
        if tms > 1.1:  # tms>=2,  0.1 make sure float tail to be  dertermined
            # N *1ms, N>=2 must be,ref spinCore document
            Inst.append((wordIdle, OptC.LONG_DELAY, tms, 1000 * unit.us))
        elif tms > 0.1:
            Inst.append((wordIdle, OptC.CONTINUE, 0, 1000 * unit.us))  # N=1 ,1ms
        if tus > 0.02:                                                       # us total Nms+us
            Inst.append((wordIdle, OptC.CONTINUE, 0, tus))
        Inst.append((wordOperating, OptC.CONTINUE, 0, RabiTime / 4.0))  # Pi/2 Pulse

        Inst.append((wordDetecting, OptC.CONTINUE, 0, DetectingTime))
        Inst.append([wordIdle, OptC.CONTINUE, 0, 1.0 * unit.us])  # compensation AOM delay
    Inst.append((wordIdle, OptC.BRANCH, 1, 1.0 * unit.us))  # delay,then next cycle
    return Inst


def SetAFG_Freq(AFG_Freq=200.0e6, AFG_Amp=1.0,
                CoolingTime=1.0 * unit.ms,  # ms
                PumpingTime=100.0 * unit.us,  # us
                DetectingTime=1000 * unit.us,  # us
                AFG_N=1048570  # 1.04857 ms AFG Pulse Length
                ):
    cycleTime = CoolingTime + PumpingTime + DetectingTime

    if cycleTime < AFG_N * unit.ns:
        raise ValueError('SpinCore cylce time less than AFG pulse Time')

    t_AFG = np.linspace(0, AFG_N * 1e-9, AFG_N)
    w = np.sin(2 * pi * AFG_Freq * t_AFG)
    setChannelWave32K(channel, w, AFG_Amp)


def StartRamsey(AFG_Freq=200.0e6, AFG_Amp=1.0,
                CoolingTime=1.0 * unit.ms,  # ms
                PumpingTime=100.0 * unit.us,  # us
                DetectingTime=1000 * unit.us,  # us
                RabiTime=18.6 * unit.us,
                t0=0.0 * unit.us,  # start us
                t=1000.0 * unit.us,  # stop us
                N=400,  # Points
                AFG_N=1048570  # 1.04857 ms AFG Pulse Length
                ):
    global Task

    t = np.linspace(t0, t, N)
    Inst = GenRamseyInst(t, RabiTime, CoolingTime, PumpingTime, DetectingTime)
    Task = SpinTask(Inst, N, Continue=True)
    SetAFG_Freq(AFG_Freq, AFG_Amp, CoolingTime, PumpingTime, DetectingTime, AFG_N)
    return t


def MeasureRamsey(tlist, LoopNum=10, data=None):
    '''

    '''

    Task.StartTask()
    N = len(tlist)
    if data is None:
        data = np.zeros(N)
        
    for i in range(LoopNum):
        data = data + Task.Read()
        
    return list(fitFringer(tlist, data))+[data]


def StopRamsey():
    global Task
    Task.StopTask()
    Task.DAQ_counter.ClearTask()


def funSin(x, x0, t0, A, B, tau):
    Y = A * np.exp(-(x - x0) / tau) * \
         np.cos(2 * pi * (x - x0) / t0) / 2.0 + B
    return Y


def fitFringer(t, y):
    '''

    '''
    data_fft = abs(np.fft.fft(y))
    pos = max(enumerate(data_fft[1:len(t) / 2]), key=lambda x: x[1])[0]
    xmin = t[0]
    xmax = t[-1]
    xscale = xmax - xmin
    t0 = xscale / (pos + 1)
    A = max(y)
    B = 1
    x0 = 0.0
    tau = 0.5e-3  # 0.5 ms
    op = [x0, t0, A, B, tau]

    popt, pcov = curve_fit(funSin, t, y, op)
    # perr = np.sqrt(np.diag(pcov)

    return (1.0 / popt[1], popt, pcov)


def CalibrateRamsey(filename='RamseyFinger.csv', LoopNum=10):
    store = pd.read_csv(filename, index_col ='index')
    last_freq = store['freq']
    last_freq = last_freq.values[-1]
    #last_freq = 210213284
    
    CoolingTime = 1.0 * unit.ms  # ms
    PumpingTime = 100.0 * unit.us  # us
    DetectingTime = 1000 * unit.us  # us
    RabiTime = 18.6 * unit.us
    t0 = 0.0 * unit.us  # start us
    t = 600.0 * unit.us  # stop us
    N = 200   # Points
    AFG_N = 1048570  # 1.04857 ms AFG Pulse Length
    AFG_Amp = 1.0  # V

    freq = last_freq - 20.0e3  # detune 20 KHz from last freq
    tlist = StartRamsey(freq, AFG_Amp,
                CoolingTime, PumpingTime, DetectingTime, RabiTime,
                t0, t, N, AFG_N)

    tlist = tlist / 1.0e9  # from spincore second to to natual second
    
    df, popt, pcov, data = MeasureRamsey(tlist, LoopNum)
    
    if df < 1.0e3 or df>50e3:
        raise ValueError('Ramsey fit run out')

    current_freq = freq + df
    tau = popt[4]

    ts = pd.to_datetime(datetime.datetime.now())
    point = pd.DataFrame([[ts,current_freq,tau]], columns=['time','freq','tau'])
    #point.to_csv(filename)
    store = store.append(point,ignore_index=True)
    store.index.name='index'
    store.to_csv(filename,columns = ['time','freq','tau'],float_format="%.10e")
    StopRamsey()
    
    return current_freq


if __name__ == '__main__':
    from pyqtgraph.Qt import QtGui, QtCore
    import pyqtgraph as pg
    
    print CalibrateRamsey()

    app = QtGui.QApplication([])
    win = pg.GraphicsWindow(title="Ramsey finger")
    win.resize(1000, 600)
    win.setWindowTitle('Plotting')
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    p6 = win.addPlot(title="Updating plot")
    curve = p6.plot(pen='y')
    curveFit = p6.plot(pen='r')

    AFG_Freq = freq2 * 1e6

    CoolingTime = 1.0 * unit.ms  # ms
    PumpingTime = 100.0 * unit.us  # us
    DetectingTime = 1000 * unit.us  # us
    RabiTime = 18.6 * unit.us
    t0 = 0.0 * unit.us  # start us
    t = 600.0 * unit.us  # stop us
    N = 200   # Points
    AFG_N = 1048570  # 1.04857 ms AFG Pulse Length

    tlist0 = StartRamsey(AFG_Freq, AFG_Amp,
                CoolingTime, PumpingTime, DetectingTime, RabiTime,
                t0, t, N, AFG_N)

    Task.StartTask()
    data = np.zeros(N)

    tlist0 = tlist0 / 1e9

    def update():
        global data
        data = data + Task.Read()
        y = data

        curve.setData(y=y, x=tlist0)
        freq, popt, pcov = fitFringer(tlist0, y)
        # perr = np.sqrt(np.diag(pcov))

        fitY = funSin(tlist0, *popt)
        curveFit.setData(y=fitY, x=tlist0)

        A = popt[2]
        B = popt[3]
        if A < 0:
            fringe = (-A / (2 * B + A))
        else:
            fringe = (A / (2 * B + A))
        print(fringe, freq)

    timer = QtCore.QTimer()
    timer.timeout.connect(update)

    timer.start(200)

    win.show()
    app.exec_()

    print("Exit")
    Task.StopTask()

    ssss = [(wordCooling, OptC.CONTINUE, 0, 4.0 * unit.ms),
            (wordCooling, OptC.BRANCH, 0, 4.0 * unit.ms)
            ]
    Task = SpinTask(ssss, 4, Continue=False)
    Task.StartTask()
