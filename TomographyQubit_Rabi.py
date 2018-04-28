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
from TomographyQubit import *
from pyqtgraph.Qt import QtGui, QtCore
from Config import wordCoolCount,wordCoolCountTriger
from setAFG import RabiTime1,RabiTime2,freq1,freq2
from scipy.linalg import expm2,eig,norm
from scipy.optimize import curve_fit,minimize

RabiPoints=200

Data=np.zeros(RabiPoints)
def funSin(x,x0,t0,A,B):
    Y=A*(1-np.cos(2*pi*(x-x0)/t0))/2.0+B
    return Y
    
def update():
    global Data,N,RabiPoints,X
    global curve, data, ptr, p6
    temp=Task.Read()
    #for i in range(N):
    #    Data+=temp[RabiPoints*i:RabiPoints*(i+1)]
    Data+=np.sum(temp.reshape(N,RabiPoints),axis=0)
    #print Data
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




app = QtGui.QApplication([])
win = pg.GraphicsWindow(title="Tomagraphy Debug")
win.resize(1000,600)
win.setWindowTitle('Plotting')
# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')
curveFit = p6.plot(pen='r')
data = np.random.normal(size=(10,1000))

timer = QtCore.QTimer()
timer.timeout.connect(update)
    
if __name__ == "__main__":
    #set pamameters
    
    N=1
    
    PumpingTime=150 #us
    DetectingTime=1000 #us
    
    RabiTime=25.0
    fre=freq1
    #########
    #WAV,timeTag=StatePrepareThenTomogra(1.0/RabiTime)
    #print WAV
    #setChannelWave32K(2,WAV)
    #a=np.linspace(0,6*pi,32768*16)
    #setChannelWave32K(2,np.sin(a))
    
    #pulsePrepare=0.0
    piTime=RabiTime/2.0  # pi flip time
    pulsePrepare=piTime/1.0
    prepare=(pulsePrepare,pi/2.0)#according to initial state
    zrot=(piTime*4,0)
    xrot=(piTime*8,0) # 通过控制 operation 区分 x操作和 z翻转
    yrot=(piTime*4,-pi/2)#此处需要核实相位
    #pulselist=(prepare,zrot,prepare,yrot)
    pulselist=(zrot,zrot)
    WAV,timeTag=GenWave(pulselist,fre*2*pi)
    setChannelWave32K(2,WAV,amplitude=1.0) # note 幅度在这里设置

    #plt.plot(WAV)
    #plt.show()
    
    InstListData=[]
    t1,t2=zip(*timeTag)
    CoolingTime=np.sum(t2)*1e-3+1.4 #ms
    print CoolingTime
    X=np.linspace(0,t2[0]+t2[1],RabiPoints)
    timeTag=[(t1[0],tt) for tt in X]
        
    #RabiPoints=4
    #Data=np.zeros(RabiPoints)
    #选择拉比测试还是tomography
    #X=[0,1,2,3]
    #timeTag=((t1[0],t2[0]),(t1[0],t2[0]+t2[1]),(t1[0],t2[0]+t2[1]/2.0),(t1[2],t2[2]+t2[3]))#z0，z1, x0, y0
    
    t1,t2=zip(*timeTag)
    t1=np.array(t1)+1.03 # compesation delay of AFG, we should triger the AFG in advance
    constime=sum(t2)
    timeTag=zip(t1,t2) # 
    #timeTag=((start1,pulse length1),(start2,pulse length2),(start3,pulse length3)...) 对应Z0,Z1,X0,Y0四组测量
    InstListData.append((wordCooling,OptC.CONTINUE,0,100*unit.ms))
    for ts,pl in timeTag: # ts: time start, pl :pulse length
        if ts<PumpingTime:
            InstListData.append((wordCooling,OptC.CONTINUE,0,CoolingTime*unit.ms))
            #InstListData.append((wordCooling,OptC.WAIT,0,1.0*unit.us))         # wait for line trigger
            InstListData.append((wordPumping,OptC.CONTINUE,0,(PumpingTime-ts)*unit.us))
            if ts==0:
                InstListData.append((wordOperatingTriger,OptC.CONTINUE,0,pl*unit.us))
            else :
                InstListData.append((wordPumpingTriger,OptC.CONTINUE,0,ts*unit.us))
                InstListData.append((wordOperating,OptC.CONTINUE,0,pl*unit.us))
        else:
            InstListData.append((wordCooling,OptC.CONTINUE,0,CoolingTime*unit.ms-(ts-PumpingTime)*unit.us))
            #InstListData.append((wordCooling,OptC.WAIT,0,1.0*unit.us))         # wait for line trigger
            InstListData.append((wordCoolingTriger,OptC.CONTINUE,0,(ts-PumpingTime+0.010)*unit.us))
            InstListData.append((wordPumping,OptC.CONTINUE,0,(PumpingTime-0.010)*unit.us))
            InstListData.append((wordOperating,OptC.CONTINUE,0,pl*unit.us))
            #InstListData.append((wordOperating,OptC.CONTINUE,0,constime*unit.us))
            
            
        InstListData.append((wordDetecting,OptC.CONTINUE,0,DetectingTime*unit.us))
        InstListData.append((wordIdle,OptC.CONTINUE,0,1.0*unit.us))
    InstListData.append((wordIdle,OptC.BRANCH,1,1.0*unit.us)) #delay,then next cycle
    
    Task=SpinTask(InstListData,RabiPoints*N,Continue=True)###
    Task.StartTask()
    
    ##
    timer.start(200)
    win.show()
    app.exec_()
    triggerstop()#停止放大器工作，继续cooling