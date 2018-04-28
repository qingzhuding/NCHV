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

def funSin(x,x0,t0,A,B):
    Y=A*(1-np.sin(2*pi*(x-x0)/t0))/2.0+B
    return Y
    
def update():
    global Data,N,RabiPoints,X
    global curve, th, ptr, p6
        
    temp=Task.Read()

    temp=th.with_cooling(temp,cooling_threshold=2.0,threshold=1.0) #读取PMT数据,阈值判断，并返回积分数据
    Data=np.average(temp.reshape(N,RabiPoints),axis=0) # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
    
    #Data=np.average(temp.reshape(N,RabiPoints),axis=0)
    
    curve.setData(y=Data,x=X)
    #y=Operator.data
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



    
if __name__ == "__main__":
    #set pamameters.3
    RabiPoints=200
    Data=np.zeros(RabiPoints)
    
    N=1
    
    PumpingTime=400 #us
    DetectingTime=1000 #us
    
    RabiTime=RabiTime2
    fre=freq2
    print RabiTime,fre
    #########
    #pulsePrepare=0.0
    piTime=RabiTime/2.0  # pi flip time
    pulsePrepare=piTime/1.0
    prepare=(pulsePrepare,pi/2.0)#according to initial state
    zrot=(piTime*5,0)
    xrot=(piTime*8,0) # 通过控制 operation 区分 x操作和 z翻转
    yrot=(piTime*6,-pi/2)#此处需要核实相位

    pulselist=(zrot,zrot)
    
    WAV,timeTag=GenWave(pulselist,fre*2*pi)
    setChannelWave32K(2,WAV,amplitude=1.0) # note 幅度在这里设置
    time.sleep(5)
    #print timeTag
    #plt.plot(WAV)
    #plt.show()
    th=ThresIntegrator(RabiPoints*N*2)  # 阈值判断对象
    
    t1,t2=zip(*timeTag)
    CoolingTime=np.sum(t2)*1e-3+1.5 #ms
    #print t2,np.sum(t2)
    X=np.linspace(0,t2[0]+t2[1],RabiPoints)
    timeTag=[(t1[0],tt) for tt in X]
    #print timeTag
    t1,t2=zip(*timeTag)
    t1=np.array(t1)+1.03 # compesation delay of AFG, we should triger the AFG in advance
    timeTag=zip(t1,t2) #  
    #timeTag=((start1,pulse length1),(start2,pulse length2),(start3,pulse length3)...) 对应Z0,Z1,X0,Y0四组测量
    #print timeTag
    
    InstListData=[]
    
    
    InstListData=GenInstList(timeTag,CoolingTime,PumpingTime,DetectingTime)#-----//-----将脉冲时序传给控制卡
    #print InstListData
    Task=SpinTask(InstListData,RabiPoints*N*2,Continue=True)###
    Task.StartTask()

    ##
    timer.start(200)
    win.show()
    app.exec_()
    triggerstop()#停止放大器工作，继续cooling