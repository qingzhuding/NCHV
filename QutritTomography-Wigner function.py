# -*- coding: utf-8 -*-
"""
Created on Fri July 29 2016

@author: jmcui
"""
from TomographyQubit import *

from pyqtgraph.Qt import QtGui, QtCore
from sympy import *
w=exp(2*pi*I/3.0)
 

def GenMultiWaves(pulse_list,SamplingRate=1000): 
    #创建多个频率的波形
    #time unit is us
    #pulselist form is [(omega,pulselen, phai),(omega,pulselen, phai),(omega,pulselen, phai),......]
    # 返回:波形 和 波形分段位置所标记的timetags,timetag的形式是（该脉冲起始时刻，该脉冲长度）
    omegaArray,pulse_len,phase=zip(*pulse_list)
    TotalTime=np.sum(pulse_len)
    TotalLen=int(TotalTime*SamplingRate)
    TotalTime=TotalLen*1.0/SamplingRate #rescale time 
    t=np.linspace(0,TotalTime,TotalLen+1) #给每个点一个时间点 0 is include in the time list
    fai=np.zeros(TotalLen+1) #将要给每个点一个相位，与所在脉冲有关，在此初始为0
    omega=np.zeros(TotalLen+1)#+omegaArray[0]#给每个点一个频率Hz
    #print omegaArray
    t0_index=0 #切片起始起索引，第一个切片的起始索引地址应该为0
    time_tags=np.zeros(len(pulse_len)+1)#定义n个脉冲的n+1个端点，作为标记点（相变点）
    for i in range(len(pulse_len)):#寻找每个脉冲对应切片，并设置相应的相位和频率
        time_tags[i+1]=np.sum(pulse_len[0:i+1])# 注意 [i+1] 寻址与[0:i+1]切片的区别，pulse_len的长度为N， i最大为N-1，i+1 最大为N， time_tags[N] 是其最后一个，而pulse_len[0:i+1]，切片0至N-1的N个数据
        t_index=int(time_tags[i+1]*SamplingRate) #计算切片结束索引时刻，time_tags[i+1]表示第i个脉冲结束时刻，或者说从0时刻到i脉冲结束时刻的时间
        fai[t0_index:t_index]=phase[i] #切片，并设置这一段的相位
        omega[t0_index:t_index]=omegaArray[i] #设置这一段的角频率
        t0_index=t_index #下一个相位区的起始索引
        
    return np.cos(omega*t-fai),zip(time_tags[0:len(pulse_len)],pulse_len)
    #这里给出了所有脉冲的timetag（起始时间，脉冲时长），然而实际用到的还需要CombinePrepareMeasureLists中的timetag来分割
    
def GenStablizerList(RabiFre,freq1=200,freq2=207.6372):
    #RabiFre unit MHz 
    #freq1 MHz
    #freq2 MHz
    a=0.39182655
    a3=1/3.0
    a2=1/2.0
    #[theta1,phai1,theta2,phai2]
    #R2*R1*STATE=[0,0,1].先打入脉冲R1[theta2,phai2]，再打入脉冲R2[theta1,phai1].所以value应该是[（freq2,theta1,phai1）,（freq2,theta2,phai2）]
    z0=[0,0,1,1]
    z1=[1,0,0,0]
    z2=[0,0,0,0]
    x0=[a,0,a2,1]
    x1=[a,a3+1,a2,a3]
    x2=[a,-a3+1,a2,-a3]
    y0=[a,0,a2,-a3]
    y1=[a,a3+1,a2,1]
    y2=[a,-a3+1,a2,a3]
    v0=[a,0,a2,a3]
    v1=[a,a3+1,a2,-a3]
    v2=[a,-a3+1,a2,1]
    #
    Angle=[[z0,z1,z2],[x0,x1,x2],[y0,y1,y2],[v0,v1,v2]]
    
    #pulsetime1=Angle[i][j][2]/RabiFre/2.0  # 注意检查 2，0，1，3 的顺序有没有写错
    #pulsetime2=Angle[i][j][0]/RabiFre/2.0
    #phi2=Angle[i][j][1]*pi
    #phi1=Angle[i][j][3]*pi
    value=[[(freq1*2*pi,Angle[i][j][2]/RabiFre/2.0,Angle[i][j][3]*pi),(freq2*2*pi,Angle[i][j][0]/RabiFre/2.0,Angle[i][j][1]*pi)] for i in range(4) for j in range(3)] # colmn and lines have Checked
    #Generate Pulses as value=[[(Pulse1),(Pulse2)] for i for j]
    return value
    
def CombinePrepareMeasureLists(prepareList,MeasureLists):
    #制备和测量合并 为（制备，测量1），（制备，测量2）（制备，测量3）......的形式
    #请注意 prepareList 应该为 [()]或者[(),()...]的列表形式，不能为单独的（），
    #MeasureLists应该为 [[()..]...] 二维列表形式，代表多个测量基[],每个测量基内有一个或者两个脉冲（）（）
    #同时生成 各个测量开始时间 TimeTag的index
    pulse_list=[]
    TimeTagIndex=np.zeros(len(MeasureLists)+1,dtype=np.int32) # fist index is zero
    ii=0
    for i,m in enumerate(MeasureLists):
        pulse_list.extend(prepareList)
        pulse_list.extend(m)
        ii=ii+len(prepareList)+len(m) # 计算GenMultiWaves中产生的TimeTag中，我们实际每组（准备+测量）对应的下标TimeTagIndex。
        TimeTagIndex[i+1]=ii
        #print i,m#不明白enumerate可以print看看就清楚了
    return pulse_list,TimeTagIndex
    
def GenTomographList(RabiFre,freq1=200,freq2=207.6372):
    #RabiFre, 拉比频率，单位 MHz， 这里的频率全都是f,需要乘以2pi才是角频率
    #无失谐共振频率，最终频率应该为omiga1=12642.8213#MHz，omiga2=omiga1+7.6372#MHz，现在写出的是发生器的频率
    #freq1=200.8213#MHz  千万注意，这里的频率全都是f,需要乘以2pi才是角频率
    #freq2=omiga1+7.6372#MHz  千万注意，这里的频率全都是f,需要乘以2pi才是角频率
    a2=1/2.0
    d1=[0,0,1,-1]
    d2=[0,0,0,0]
    d3=[a2,0,1,1]
    d4=[a2,0,1,-a2]
    d5=[0,0,a2,1]
    d6=[0,0,a2,-a2]
    d7=[a2,0,0,0]
    d8=[a2,-a2,0,0]
    #规整为列表
    tomo_angle=[d1,d2,d3,d4,d5,d6,d7,d8]
    tomo_value=[[(freq1*2*pi,tomo_angle[i][2]/RabiFre/2.0,tomo_angle[i][3]*pi),(freq2*2*pi,tomo_angle[i][0]/RabiFre/2.0,tomo_angle[i][1]*pi)] for i in range(8)] # Please Check This Line
    return tomo_value

def StatePrepareThenTomogra(prepareList,RabiTime,freq1=100.0,freq2=207.6372):#freq2=207.6372
    #RabiTime unit is us, flip 2*pi phase ,drive |0> state to -|0> state
    
    MeasureLists=GenTomographList(1.0/RabiTime,freq1,freq2)
    
    print MeasureLists
    pulse_list,TimeTagIndex=CombinePrepareMeasureLists(prepareList,MeasureLists)
    
    waves,Tags=GenMultiWaves(pulse_list)#返回所有脉冲的timetag
    t,dt=zip(*Tags)#t为起始时间，dt为脉冲时长
    TimeTags=[(t[TimeTagIndex[i]],sum(dt[TimeTagIndex[i]:TimeTagIndex[i+1]])) for i in range(len(MeasureLists))] #[(time,pulse_len) for i in N]
    return waves, TimeTags

    
def detect_and_fidelity(preparedrou):
    # 打印结果，重复运行200次
    k=len(TimeTags)  #输入脉冲数
    Data=np.zeros(k)
    plus=np.zeros(k)
    for j in range(200):
        temp=Task.Read()
        for i in range(N):
            for a in range(k):#对读数进行0，1赋值判断
                #print temp[k*i],temp[k*i+1],temp[k*i+2],temp[k*i+3]
                if temp[k*i+a]<0.9:#阈值设定并赋值0和1
                    plus[a]=1.0 #测到暗态为1
                else:
                    plus[a]=0.0
            #Data+=temp[k*i:k*i+4]
            Data+=plus[0:k]
        proba=0.01*Data/(1.0*(j+1))
        
        detectestate=(unitymatrix+sigmax*(2*proba[2]-1)+sigmay*(2*proba[3]-1)+sigmaz*(proba[0]-proba[1]))/2.0
        tra=np.trace(np.dot(detectestate,detectestate))#求rou平方的迹
        fidelity0=np.trace(np.dot(detectestate,preparestate))
        #fidelity1=np.trace(np.sqrt(np.dot(np.dot(np.sqrt(preparestate),detectestate),np.sqrt(preparestate))))**2
        print np.real(fidelity0),np.real(tra)
        print detectestate
        print preparestate
        #time.sleep(0.05)
    #timer.start(200)
    #print(TimeTags)

if __name__=='__main__':
    RabiTime=0.10 #us
    #Initial state 0
    freq1=200.0
    freq2=207.6372
    CoolingTime=5 #ms
    PumpingTime=200 #us
    DetectingTime=1000 #us
    RabiFre=1.0/RabiTime
    amplitude=1.0
    N=100
    prepareList=[(freq1*2*pi,RabiTime/4.0,0.0)] # Note the 2*pi
    
    waves,TimeTags=StatePrepareThenTomogra(prepareList,RabiTime)
    InstListData=GenInstList(TimeTags,CoolingTime,PumpingTime,DetectingTime)
    #Task=SpinTask(InstListData,len(TimeTags)*N,Continue=True)#配置采集卡
    #Task.StartTask()
    #setChannelWave32K(2,waves,amplitude)
    #detect_and_fidelity(preparestate)
    print(len(TimeTags))
    #get all results from wignerfunctiontopulselist with definite x,y
    import pickle as p,pprint
    wignerf=r'D:wignerfunctiontopulselistx-%fy-%f.pkl'%(0.100000,0.100000)
    with open(wignerf, 'r') as f: 
        proj1,Aq,x_point,y_point,eigenmixed,mixprobability,mixpurestate,pre_angle_list,rouxyz=p.load(f)
    rouxyz=rouxyz.subs('a',x_point).subs('b',y_point)
    print x_point+12
    print pre_angle_list
    print rouxyz
    plt.plot(waves)
    plt.show()
