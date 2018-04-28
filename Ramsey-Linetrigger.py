# -*- coding: utf-8 -*-
"""
Created on Wed Dec 31 14:15:43 2014

@author: jmcui
"""
import numpy as np
import visa,os,sys, time
import matplotlib.pyplot as plt
from cmath import *
from scipy.linalg import expm2,eig,norm
from scipy.optimize import curve_fit,minimize
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

SITE_ROOT = os.path.dirname(os.path.realpath(__file__))
PARENT_ROOT=os.path.abspath(os.path.join(SITE_ROOT, os.pardir))

sys.path.append(PARENT_ROOT)
import spinpy.OperationCode as OptC
import spinpy.unit as unit
from StateManipulation import *
from Config import wordCooling,wordPumping,wordOperating,wordDetecting,wordIdle


def mainFun():
    global tlist0,Operator
    CoolingTime=2.0 # ms
    PumpingTime=400.0 #us
    DetectingTime=1000 # u45.9s
    PiPulse=19.6/2.0 #us operation time M=0 41.8 M=1 23.73, M=-1 39.7
    
    t0=0#start us
    t1=1000  #stop us
    N=300   #Points

    #E8257D.write(':FREQ:FIX 9.644818870GHZ')
    InstListData=[[]]
    tlist0=np.linspace(t0,t1,N)
    for t in tlist0:
        tms=int(t/1000) #ms
        tus=t-tms*1000
        InstListData[0].append((wordCooling,OptC.CONTINUE,0,CoolingTime*unit.ms))
        #InstListData[0].append((wordCooling,OptC.WAIT,0,1.0*unit.us))         # wait for line trigger
        InstListData[0].append((wordPumping,OptC.CONTINUE,0,PumpingTime*unit.us))
        InstListData[0].append((wordOperating,OptC.CONTINUE,0,PiPulse/2.0*unit.us))# pi/2
        if tms>1.1: #tms>=2,0.1 make sure float tail to be  dertermined 
            InstListData[0].append((wordIdle,OptC.LONG_DELAY,tms,1000*unit.us)) # N *1ms, N>=2 must be,ref spinCore document
        elif tms>0.1:
            InstListData[0].append((wordIdle,OptC.CONTINUE,0,1000*unit.us)) # N=1 ,1ms
        if tus>0.02:                                                       # us total Nms+us
            InstListData[0].append((wordIdle,OptC.CONTINUE,0,tus*unit.us))
        InstListData[0].append((wordOperating,OptC.CONTINUE,0,PiPulse/2.0*unit.us))#Pi/2 Pulse
        InstListData[0].append((wordDetecting,OptC.CONTINUE,0,DetectingTime*unit.us))
        #InstListData[0].append([wordIdle,OptC.CONTINUE,0,1.0*unit.us]) #delay
    InstListData[0].append((wordIdle,OptC.BRANCH,0,1.0*unit.us)) #delay,then next cycle
    Operator=StateOperator(N,InstListData)
    Operator.start(TimeOut=20.0)


app = QtGui.QApplication([])
win = pg.GraphicsWindow(title="Basic plotting examples")
win.resize(1000,600)
win.setWindowTitle('Plotting')
# Enable antialiasing for prettier plots
pg.setConfigOptions(antialias=True)
p6 = win.addPlot(title="Updating plot")
curve = p6.plot(pen='y')
curveFit = p6.plot(pen='r')
data = np.random.normal(size=(10,1000))
p=0
def funSin(x,x0,t0,A,B):
    Y=A*(1-np.cos(2*pi*(x-x0)/t0))/2.0+B
    return Y
def update():
    global curve, data, ptr, p6,p
    y=Operator.data
    curve.setData(y=y,x=tlist0)
    
    #fit Rammesy
    data_fft=abs(np.fft.fft(y))
    pos=max(enumerate(data_fft[1:len(tlist0)/2]), key = lambda x: x[1])[0]
    xmin=tlist0[0]
    xmax=tlist0[-1]
    xscale=xmax-xmin
    t0= xscale/(pos+1)
    A=max(y)
    B=1
    x0=0
    op=[x0,t0,A,B]
 
    popt, pcov = curve_fit(funSin, tlist0,Operator.data,op)
    perr = np.sqrt(np.diag(pcov))

    fitY=funSin(tlist0,*popt)
    A=popt[2]
    B=popt[3]
    if A<0:
        fringe=(-A/(2*B+A))
    else:
        fringe=(A/(2*B+A))
    #print fringe,popt[1], 1000/popt[1]
    p=p+1
    if (p%10)==0:
        print 1000/popt[1]
    curveFit.setData(y=fitY,x=tlist0)
timer = QtCore.QTimer()
timer.timeout.connect(update)

if __name__=='__main__':
    global Operator
    mainFun()
    timer.start(200)
    win.show()
    app.exec_()
    print "Exit"
    Operator.stop()
    #停止放大器工作，继续cooling
    ssss=[(wordCooling,OptC.BRANCH,0,4.0*unit.ms)]
    SpinTask(ssss,4,Continue=False)
