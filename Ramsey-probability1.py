# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2014

@author: jmcui
"""
import numpy as np
import visa,os,sys
from struct import pack
import matplotlib.pyplot as plt
import time
from math import *
from scipy.linalg import expm2,eig,norm
from datetime import datetime
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

#import matplotlib

SITE_ROOT = os.path.dirname(os.path.realpath(__file__))
PARENT_ROOT=os.path.abspath(os.path.join(SITE_ROOT, os.pardir))

sys.path.append(PARENT_ROOT)
import spinpy.OperationCode as OptC
import spinpy.unit as unit
from StateManipulation import SpinTask
import Config
from pyqtgraph.Qt import QtGui, QtCore
from Config import wordCoolCount,wordCoolCountTriger
from setAFG import RabiTime1,RabiTime2,freq1,freq2
from scipy.linalg import expm2,eig,norm
from scipy.optimize import curve_fit,minimize
from TomographyQubitCoolingValidity import *
import random

app = QtGui.QApplication([])
win = pg.GraphicsWindow(title="RABI PROBABILITY")
win.resize(1000,600)
win.setWindowTitle('Plotting')
# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')
curveFit = p6.plot(pen='r')
data = np.random.normal(size=(10,1000))

def mainFun():
    global tlist,Operator,CoolingTime,PumpingTime,DetectingTime,PiPulse
    InstListData=[]
    for t in tlist:
        tms=int(t/1000) #ms
        tus=t-tms*1000
        InstListData.append((wordCoolCount,OptC.CONTINUE,0,ctypes.c_double(CoolingTime*unit.ms).value))
        InstListData.append((wordCooling,OptC.WAIT,0,ctypes.c_double(1.0*unit.us).value))         # wait for line trigger
        InstListData.append((wordPumping,OptC.CONTINUE,0,ctypes.c_double((PumpingTime-1.03)*unit.us).value))
        InstListData.append((wordPumpingTriger,OptC.CONTINUE,0,ctypes.c_double(1.03*unit.us).value))
        InstListData.append((wordOperating,OptC.CONTINUE,0,ctypes.c_double((PiPulse/2.0)*unit.us).value))# pi/2
        if tms>1.1: #tms>=2,0.1 make sure float tail to be  dertermined 
            InstListData.append((wordIdle,OptC.LONG_DELAY,tms,ctypes.c_double(1000*unit.us).value)) # N *1ms, N>=2 must be,ref spinCore document
        elif tms>0.1:
            InstListData.append((wordIdle,OptC.CONTINUE,0,ctypes.c_double(1000*unit.us).value)) # N=1 ,1ms
        if tus>0.02:                                                       # us total Nms+us
            InstListData.append((wordIdle,OptC.CONTINUE,0,ctypes.c_double(tus*unit.us).value))
        InstListData.append((wordOperating,OptC.CONTINUE,0,ctypes.c_double((PiPulse/2.0)*unit.us).value))#Pi/2 Pulse
        InstListData.append((wordDetecting,OptC.CONTINUE,0,ctypes.c_double(DetectingTime*unit.us).value))
        InstListData.append([wordIdle,OptC.CONTINUE,0,2.0*unit.us]) #delay
    InstListData.append((wordIdle,OptC.BRANCH,0,ctypes.c_double(2.0*unit.us).value)) #delay,then next cycle
    return InstListData

    

app = QtGui.QApplication([])
win = pg.GraphicsWindow(title="RAMSEY PROBABILITY")
win.resize(1000,600)
win.setWindowTitle('Plotting')
# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')
curveFit = p6.plot(pen='r')
data = np.random.normal(size=(10,1000))

def funSin(x,x0,t0,A,B):
    Y=A*(1-np.sin(2*pi*(x-x0)/t0))/2.0+B
    return Y
    
def update():
    global Data,N,tlist,X
    global curve, th, ptr, p6
    temp=Task.Read()
    temp=th.with_cooling(temp,cooling_threshold=20.0,threshold=1.0) #读取PMT数据,阈值判断，并返回积分数据
    Data=1-temp # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
    #Data=np.average(temp.reshape(N,RabiPoints),axis=0)
    
    #curve.setData(y=Data,x=X)
    y=Data
    curve.setData(y=y,x=X)
    #print(Operator.data)
    #fit Rammesy
    data_fft=abs(np.fft.fft(y))
    pos=max(enumerate(data_fft[1:len(X)/2]), key = lambda x: x[1])[0]
    xmin=X[0]
    xmax=X[-1]
    xscale=xmax-xmin
    t0= xscale/(pos+1)
    A=max(y)
    B=1
    x0=0
    op=[x0,t0,A,B]
 
    popt, pcov = curve_fit(funSin, X,Data,op)
    perr = np.sqrt(np.diag(pcov))

    fitY=funSin(X,*popt)
    A=popt[2]
    B=popt[3]
    if A<0:
        fringe=(-A/(2*B+A))
    else:
        fringe=(A/(2*B+A))
    print fringe,popt[1], 1000/popt[1]
    curveFit.setData(y=fitY,x=X)
    
timer = QtCore.QTimer()
timer.timeout.connect(update)

if __name__=='__main__':

    CoolingTime=4.0 # ms
    PumpingTime=400.0 #us
    DetectingTime=1000 # u45.9s
    RabiTime=RabiTime2
    fre=freq2
    
    PiPulse=RabiTime*0.5 #us operation time M=0 40.2 M=1 38.1

    N=200   #Points
    t1=600
    #E8257D.write(':FREQ:FIX 9.644818870GHZ')
    Data=np.zeros(N)
    
    th=ThresIntegrator(N*2)
    InstListData=[]
    
    pulselist=((t1,0),(t1,0))

    tlist=np.linspace(0,2*t1,N)
    X=tlist
    InstListData=mainFun()
    #print InstListData
    WAV,timeTag=GenWave(pulselist,fre*2*pi)
    setChannelWave32K(2,WAV,amplitude=1.0) # note 幅度在这里设置
    time.sleep(2)
    Task=SpinTask(InstListData,N*2,Continue=True)
    Task.StartTask()

    timer.start(200)
    win.show()
    app.exec_()
    print "Exit"
    #停止放大器工作，继续cooling
    triggerstop()#停止放大器工作，继续cooling
