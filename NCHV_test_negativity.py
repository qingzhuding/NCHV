# -*- coding: utf-8 -*-
"""
Created on 08 17 2016

@author: heran
"""
from TomographyQutritCoolingValiditytest import *
import scipy.linalg as linalg
from sympy import nan
from setAFG import RabiTime1,RabiTime2,freq1,freq2
import pandas as pd
pd.options.display.max_rows =200
pd.options.display.max_columns =200
pd.set_option('max_colwidth',200)
#精度，但是不会设置
pd.set_option('precision',4)
np.set_printoptions(precision=4)


def single_matrix_decomp(M):#------这个是测量角度列表，作用顺序为RY*RX*RYR*vector=[0,0,1]
    M=np.array(M).astype(np.complex)
    eigenvalue,vec=linalg.eig(M)
    vect=[vec[:,i] for i in range(len(eigenvalue))]
    #print vect
    pulse=[pre_angle(vect[i])+np.array([0,1,0,1,0,1]) for i in range(len(eigenvalue))]#同时获得该力学量的三个本征态的三组测量角度列表,加np.array([0,1,0,1,0,1])是为了获得逆矩阵
    return eigenvalue,vect,pulse
def matrix_list_decomp(mlist):
    prob_all=[]
    eigenvecters_all=[]
    pluse_all=[]
    numb=len(mlist)
    for i in range(numb):
        eigvalue,vect,pul=single_matrix_decomp(mlist[i])
        prob_all.append(eigvalue)
        eigenvecters_all.append(vect)
        pluse_all.append(pul)
    return prob_all,eigenvecters_all,pluse_all
    
    
def GenmeasureList(para):
    RabiTime1,RabiTime2,freq1,freq2=para
    mea_angle=[]
    #global Aq_pulse_list
    num=len(Aq_pulse_list)*len(Aq_pulse_list[0])#num=行×列
    #print num
    for a in Aq_pulse_list:
        mea_angle+=a
    measure_value=[[(freq1*2*pi,(mea_angle[i][4])*RabiTime1/2.0,mea_angle[i][5]*pi),
    (freq2*2*pi,(mea_angle[i][0])*RabiTime2/2.0,mea_angle[i][1]*pi),
    (freq1*2*pi,mea_angle[i][2]*RabiTime1/2.0,mea_angle[i][3]*pi)] for i in range(num)] # Please Check This Line
    return measure_value
        
def StatePrepareThenmeasure(prepareList,para):
    RabiTime1,RabiTime2,freq1,freq2=para
    #RabiTime unit is us, flip 2*pi phase ,drive |0> state to -|0> state；pi phase to |1>state
    MeasureLists=GenmeasureList(para)
    pulse_list,TimeTagIndex=CombinePrepareMeasureLists(prepareList,MeasureLists)
    waves,Tags=GenMultiWaves(pulse_list)#返回所有脉冲的timetag
    t,dt=zip(*Tags)#t为起始时间，dt为脉冲时长
    TimeTags=[(t[TimeTagIndex[i]],sum(dt[TimeTagIndex[i]:TimeTagIndex[i+1]])) for i in range(len(MeasureLists))] #[(time,pulse_len) for i in N]
    return waves,TimeTags

def chektimetag(timeTag):

    print'--------timeTag--------'
    print timeTag
    
    #验证timetag正确
    def timetagcheck(stateangle,mea_angle_list):
        a=[]
        b=0
        for mea_angle in mea_angle_list:
            b=((stateangle[0]+mea_angle[0])*RabiTime2+(stateangle[2]+mea_angle[2])*RabiTime1+(stateangle[4]+mea_angle[4])*RabiTime1)*0.5
            a.append(b)
        return a
    
    print [i[1] for i in timeTag]
    listofaq=[]
    for i in Aq_pulse_list:
        listofaq+=i
    print listofaq
    print timetagcheck(pre_angle_list[0],listofaq)
    
    
def measure(prepareList,para):
    global turn
    RabiTime1,RabiTime2,freq1,freq2=para #参数设置
    waves,timeTag=StatePrepareThenmeasure(prepareList,para)#产生波形和TimeTag
    k=len(timeTag)  #输入脉冲数
    t0,dt=zip(*timeTag)
    t0=np.array(t0)+1.030#-----//-----compesation delay of AFG, we should triger the AFG in advance
    timeTag=zip(t0,dt)
    plusetimesum=np.sum([i[1] for i in timeTag])#脉冲总时长
    print'--------timeTag--------'
    print mp.matrix(np.array(timeTag).T)
    print '%d pulse time sum:%s us\n'%(k,plusetimesum)
    print '-----measure data--------------'
    
    #如果总脉冲时间小于cooling时间，会发生时序错误，下式保证冷却时间不小于4ms
    CoolingTime=int((plusetimesum)/1000.0)+1.5  #至少需要4ms的冷却,先于setChannelWave32K定义
    setChannelWave32K(2,waves,amplitude=1.0)#写入波形
    time.sleep(2)
    InstListData=GenInstList(timeTag,CoolingTime,PumpingTime,DetectingTime)
    Task=SpinTask(InstListData,k*N*2,Continue=True)#-----//-----k组测量，依次进行一次，然后整体重复N次。每次测量包含 cooling 和 detect 两次计数
    #chektimetag(timeTag)
    
    # 打印结果，重复运行200次
    th=ThresIntegrator(k*N*2)  #阈值判断对象
    Task.StartTask()
    
    every_prob=[]#记录所有的探测概率
    for j in range(100):#********************Repeat numbers
        t0=time.time()
        data1=Task.Read(TimeOut=30.0)
	#data1=simulationtask(k*N*2)
        #print data1
        print mp.mpf(time.time()-t0)#读取循环100次的时间，约为4.5秒
        temp=th.with_cooling(data1,cooling_threshold=20.0,threshold=1.0) #读取PMT数据,阈值判断，并返回积分数据
        prob=np.average(temp.reshape(N,k),axis=0) # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
        #print 'prob:%s '%prob.reshape(k/state_numb,state_numb)
        difference=np.array(prob)-np.array(probability_theory)
        print 'turn%s/3 the %s time, effective count: %s '%(turn+1,j+1,np.sum(th.eff_count))
        print 'max error: NO.%s  %s , abs average:%s  '%(np.argmax(np.abs(difference)),np.max(np.abs(difference)),np.average(np.abs(difference)))
        #print 'prob     %s'%prob.reshape(k/state_numb,state_numb)
        #print 'theory    %s'%np.array(probability_theory).reshape(k/state_numb,state_numb)
        #print 'difference%s'%(np.array(prob)-np.array(probability_theory)).reshape(k/state_numb,state_numb)
        sumprob=np.sum(prob.reshape(k/state_numb,state_numb),axis=1)
        
        print 'sum prob  %s    fangcha:  %s   average:  %s'%(sumprob,np.std(sumprob),np.average(sumprob))
        #pd.set_option('expand_frame_repr', False)#是否换行
        df2=pd.DataFrame({ 'difference' :difference,'prob' : prob,'theory' : probability_theory},
        columns=['prob', 'theory', 'difference'])
        print df2.T
        #print df2
        every_prob.append(prob)#存储每一步的概率值,以备后用
        
    N_temp=temp.reshape(N,k)
    #print len(temp)
    #temp.reshape(N,k),只需要对着里面的N个元素（每个元素是k个值）分别求getmeanvalue（）,就得到N个值，再用np.std即可
    return prob,every_prob,N_temp
    
def getmeanvalue(Aq_prob,rho3prob,values):
    values=np.array(values)
    #print values
    rho3prob=np.array(rho3prob)
    Aq_prob=np.array(Aq_prob)
    r=values.T*rho3prob
    v=r[:,0]+r[:,1]+r[:,2]
    w=np.array(v).reshape(len(Aq_prob),3)
    u=np.sum(w*Aq_prob,axis=1)
    return np.real(u)
    
#check whether the measure angle are right,if all [0,0,1],right
def mea_angle_check(mea_angle,vect):
    th1,phi1,th2,phi2,th3,phi3=mea_angle*pi
    RYR=np.array([[cos(th3/2.0),0,1.0j*(cos(phi3)-1.0j*sin(phi3))*sin(th3/2.0)],[0,1,0],[1.0j*(cos(phi3)+1.0j*sin(phi3))*sin(th3/2.0),0,cos(th3/2.0)]])
    RYS=np.array([[cos(th2/2.0),0,1.0j*(cos(phi2)-1.0j*sin(phi2))*sin(th2/2.0)],[0,1,0],[1.0j*(cos(phi2)+1.0j*sin(phi2))*sin(th2/2.0),0,cos(th2/2.0)]])
    RXS=np.array([[1,0,0],[0,cos(th1/2.0),1.0j*(cos(phi1)-1.0j*sin(phi1))*sin(th1/2.0)],[0,1.0j*(cos(phi1)+1.0j*sin(phi1))*sin(th1/2.0),cos(th1/2.0)]])
    zer0_check=RYS.dot(RXS.dot(RYR.dot(vect)))
    return zer0_check

if __name__ == "__main__":#-----//------流程图详解见下
    
    #---------------------------set pamameters------------------------------------
    PumpingTime=400 #us
    DetectingTime=1000 #us#-----//------设置冷却，泵浦，探测时间，可以根据实际需求进行优化
    N=100
    print 'N=%s'%N
    x_point,y_point=[0.1,-0.04]
    x_point,y_point=[0.15,0.19]
    print '-----pamameters--------'
    print 'x_point=%s,y_point=%s '%(x_point,y_point)
    print 'RabiTime1,RabiTime2,freq1,freq2%s\n'%[RabiTime1,RabiTime2,freq1,freq2]
    
    import pickle as p,pprint
    wignerf=r'D:\Users\Desktop\HERAN\IonControl\IonControl\NCHV\wignerfunction1\wignerfunctiontopulselistx-%.6fy-%.6f.pkl'%(x_point,y_point)
    with open(wignerf, 'r') as f: 
        proj1,Aq,q,eigenmixed,mixprobability,mixpurestate,pre_angle_list,rouxyz,rou=p.load(f)
    #print proj1
    rouxyz=rouxyz.subs('a',x_point).subs('b',y_point)#得到赋值之后的密度矩阵
    rouxyz=np.array(rouxyz).astype(np.complex)
    state_numb=len(mixpurestate)#其实就是3
    #vector to proj.  mixpurestate是态,mixpurerou是密度矩阵
    q=[np.array(a).astype(np.complex) for a in q]
    mixpurestate=[np.array(a).astype(np.complex) for a in mixpurestate]
    mixpurerou=[np.array(mixpurestate[i]*mixpurestate[i].conjugate().T).astype(np.complex) for i in range(state_numb)]
    mixprobability=np.array(mixprobability).astype(np.float)
    pre_angle_list=[np.array(i).astype(np.float) for i in pre_angle_list]

    Aq_measure_list=[np.array(Aq[i][j]).astype(np.complex) for i in range(3) for j in range(3)]#所有需要测量的矩阵的一维列表
    Aq_measure_list=[np.array(Aq[0][0]).astype(np.complex) for i in range(3) for j in range(3)]
    mea_num=state_numb*len(Aq_measure_list)#脉冲数
    Aq_prob_list,Aq_eigenvecter_list,Aq_pulse_list=matrix_list_decomp(Aq_measure_list)#得到每个矩阵对应的三个本征态及其组合概率，测量角度列表
    theory_result=[np.trace(np.dot(rouxyz,a)) for a in Aq_measure_list]
    print '------3 eigenstates------'
    print pd.Series(mixpurestate)
    print '------3 eigenvalues------'
    print pd.Series(mixprobability)    
    print '----Aq_measure_list----------'
    print pd.Series(Aq_measure_list)
    print '----------Aq_prob_list-------------------------'
    print pd.Series(Aq_prob_list)
    print '----------Aq_eigenvecter_list-------------------------'
    print pd.Series(Aq_eigenvecter_list)
    print '----------Aq_pulse_list-------------------------'
    print pd.Series(Aq_pulse_list)

    result=[]#记录27个测量的三组最终结果
    para=[RabiTime1,RabiTime2,freq1,freq2]#打包参数
    all_prob=[[],[],[]]#好像还需要累积所有100次的每次的结果。。
    N_independent_results=[]#记录N次独立测量结果，用于求errorbar
    
    theory_prob_of_each_eigenvects=[[],[],[]]
    asd=[]
    c=[]

    for i in range(len(Aq_measure_list)):
        #print i
        #print '----------------------------------'
        for j in range(state_numb):
            a=Aq_pulse_list[i][j]
            b=Aq_eigenvecter_list[i][j]
            #print a,b
            theory_prob_of_each_eigenvects[0].append((np.abs(b.dot((mixpurestate[0].T[0]).conjugate())))**2)
            theory_prob_of_each_eigenvects[1].append((np.abs(b.dot((mixpurestate[1].T[0]).conjugate())))**2)
            theory_prob_of_each_eigenvects[2].append((np.abs(b.dot((mixpurestate[2].T[0]).conjugate())))**2)
            asd.append((np.abs(b.dot(mixpurestate[0].conjugate())))**2)
            print mea_angle_check(a,b)
        #print Aq_eigenvecter_list[i][2]
    #print np.array(theory_prob_of_each_eigenvects)


    for i in range(state_numb):
        turn=i
        probability_theory=theory_prob_of_each_eigenvects[i]
        print '------input state %d------'%(i+1)
        print 'state %d=%s'%(i+1,mixpurestate[i].T[0])
        print 'angle %d=%s\n'%(i+1,pre_angle_list[i])
        prepareList=[(freq1*2*pi,pre_angle_list[i][2]*RabiTime1/2.0,pi*pre_angle_list[i][3]),
        (freq2*2*pi,pre_angle_list[i][0]*RabiTime2/2.0,pi*pre_angle_list[i][1]),(freq1*2*pi,pre_angle_list[i][4]*RabiTime1/2.0,pi*pre_angle_list[i][5])] # Note the 2*pi
        prob,i_prob,N_temp=measure(prepareList,para)
        
        N_independent_results.append(N_temp)#将N组k个独立测量值记录下来
        all_prob[i].append(i_prob)
        result.append(prob)
        
    #求9个期望值的error bar
    ww=[]
    for i in range(N):
        for j in range(state_numb):
            ww.append(N_independent_results[j][i])#将记录下来的结果重新排列为三次分量测量的27个值放在一起，共计N组
    mm=np.array(ww).reshape(N,mea_num*state_numb)#分成N组，分别对应N次独立测量结果
    qq=[list(m.reshape(state_numb,mea_num)) for m in mm]#将每组重排列为三个array，每个array对应三个态矢下的27个测量结果
    yy=[getmeanvalue(Aq_prob_list,mixprobability,opp) for opp in qq]#计算所有的独立期望值，N组，每组有9个，分别对应9个测量矩阵
    samples=np.array(yy).T#将yy中的结果转置为采样列表，每一个array包含每个测量矩阵的N次采样结果
    errorbar=[np.std(s) for s in samples]#得到errorbar
    average2=[np.average(s)for s in samples]#获得N次测量均值，应与之前算出的结果相同
    
    print 'errorbar'
    print errorbar
    print 'average2'
    print average2#should be the same as value_of_Aq
    
    value_of_Aq=getmeanvalue(Aq_prob_list,mixprobability,result)
    print '------9 Aq values of x=%f, y=%f--------------'%(x_point,y_point)
    print value_of_Aq
    print '----------theory_result---------------'
    print np.real(theory_result)
    min=np.argmin(average2)
    ABmin=[x_point,y_point,average2[min],errorbar[min],theory_result[min]]
    print '-----------ABmin-------------'
    print ABmin
    timerecord=time.strftime("date%y-%m-%d time%H-%M", time.localtime())
    data=(ABmin,theory_prob_of_each_eigenvects,average2,theory_result,errorbar,N_independent_results,Aq_prob_list,Aq_eigenvecter_list,mixpurestate,mixpurerou,mixprobability,timerecord)
    import pickle as p,pprint
    wignerf=r'D:\Users\Desktop\HERAN\IonControl\IonControl\NCHV\NCHVnegativity1\NCHVnegativityofx-%.6fy-%.6f.pkl'%(x_point,y_point)
    with open(wignerf, 'w') as f:
        p.dump(data,f)
    triggerstop()
    plt.errorbar(range(len(Aq_measure_list)),average2,yerr=errorbar,fmt='o',color='r')
    plt.plot(range(len(Aq_measure_list)),theory_result,'D')
    plt.show()
    

    '''
    print '----------Aq_prob_list-------------------------'
    print Aq_prob_list
    print '----------Aq_eigenvecter_list-------------------------'
    print Aq_eigenvecter_list
    print '----------Aq_pulse_list-------------------------'
    print Aq_pulse_list
    '''
    