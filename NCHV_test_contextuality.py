# -*- coding: utf-8 -*-
"""
Created on 08 17 2016

@author: heran
"""
from TomographyQutritCoolingValiditytest import *
from NCHV_test_negativity import *
import scipy.linalg as linalg
from sympy import nan
from wignerfunctiontopulselist import x_point,y_point
from setAFG import RabiTime1,RabiTime2,freq1,freq2

# -*- coding: utf-8 -*-
"""
Created on 08 17 2016

@author: heran
"""
from TomographyQutritCoolingValiditytest import *
import scipy.linalg as linalg
from sympy import nan
from wignerfunctiontopulselist import x_point,y_point

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
    return waves, TimeTags
def setTask(timeTag):
    t0,dt=zip(*timeTag)
    t0=np.array(t0)+1.03 #-----//-----compesation delay of AFG, we should triger the AFG in advance
    timeTag=zip(t0,dt)
    InstListData=GenInstList(timeTag,CoolingTime,PumpingTime,DetectingTime)
    k=len(timeTag)  #输入脉冲数
    print'--------timeTag--------'
    print timeTag
    Task=SpinTask(InstListData,k*N*2,Continue=True)#-----//-----k组测量，依次进行一次，然后整体重复N次。每次测量包含 cooling 和 detect 两次计数
    return Task
def measure(prepareList,para):
    RabiTime1,RabiTime2,freq1,freq2=para #参数设置
    waves,timeTag=StatePrepareThenmeasure(prepareList,para)#产生波形和TimeTag
    setChannelWave32K(2,waves,amplitude=1.0)#写入波形
    Task=setTask(timeTag)#设置时序
    k=len(timeTag)  #输入脉冲数
    #print 'k:   %d'%k
    every_prob=[]#记录所有的探测概率
    # 打印结果，重复运行200次
    th=ThresIntegrator(k*N*2)  #阈值判断对象
    Task.StartTask()
    eff=[]#record which count is valid
    count=[]
    cooling_threshold=2
    counting_threshold=1
    for j in range(40):#********************Repeat numbers
        t0=time.time()
        data=Task.Read()
        #记录所有探测结果与是否有效
        cooling_counts_once,count_once=data.reshape(len(data)/2,2).T # seperate cooling and detection counts
        eff_once=(cooling_counts_once>cooling_threshold)*1.0
        eff.append(eff_once)
        count.append(count_once)
        print mp.mpf(time.time()-t0)#读取循环100次的时间，约为4.5秒
        temp=th.with_cooling(data,cooling_threshold,counting_threshold) #读取PMT数据,阈值判断，并返回积分数据
        prob=np.average(temp.reshape(N,k),axis=0) # 重整为N行k列数据，axis=0对N行求平均值，获得k个平均数
        print j+1,np.sum(th.effective_cycles)
        #print prob.reshape(3,9)
        every_prob.append(prob)#存储每一步的概率值,以备后用
        N_temp=temp.reshape(N,k)
    #print len(temp)
    #temp.reshape(N,k),只需要对着里面的N个元素（每个元素是k个值）分别求getmeanvalue（）,就得到N个值，再用np.std即可
    return prob,every_prob,N_temp,eff,count
    
def getmeanvalue(Aq_prob,rho3prob,values):
    values=np.array(values)
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
    CoolingTime=4 #ms
    PumpingTime=400 #us
    DetectingTime=1000 #us#-----//------设置冷却，泵浦，探测时间，可以根据实际需求进行优化
    N=10
    
    x_point,y_point=[-0.11,0.2]
    import pickle as p,pprint
    wignerf=r'D:\Users\Desktop\HERAN\IonControl\IonControl\NCHV\wignerfunction1\wignerfunctiontopulselistx-%.6fy-%.6f.pkl'%(x_point,y_point)
    with open(wignerf, 'r') as f: 
        proj1,Aq,q,eigenmixed,mixprobability,mixpurestate,pre_angle_list,rouxyz,rou=p.load(f)
    #print proj1
    x_point=np.array(x_point).astype(np.float)
    y_point=np.array(y_point).astype(np.float)
    rouxyz=rouxyz.subs('a',x_point).subs('b',y_point)#得到赋值之后的密度矩阵
    rouxyz=np.array(rouxyz).astype(np.complex)
    state_numb=len(mixpurestate)#其实就是3
    #vector to proj.  mixpurestate是态,mixpurerou是密度矩阵
    q=[np.array(a).astype(np.complex) for a in q]
    mixpurestate=[np.array(a).astype(np.complex) for a in mixpurestate]
    mixpurerou=[np.array(mixpurestate[i]*mixpurestate[i].conjugate().T).astype(np.complex) for i in range(state_numb)]
    mixprobability=np.array(mixprobability).astype(np.float)
    pre_angle_list=[np.array(i).astype(np.float) for i in pre_angle_list]
    

    Aq_measure_list=[state_numb**3*np.eye(state_numb)-np.array(Aq[i][j]).astype(np.complex) for i in range(3) for j in range(3)]#所有需要测量的矩阵的一维列表
    #Aq_measure_list.append(np.eye(state_numb))#添加在两个qutrit的态相同的情况下，对直积项I的测量
    Aq_prob_list,Aq_eigenvecter_list,Aq_pulse_list=matrix_list_decomp(Aq_measure_list)#得到每个矩阵对应的三个本征态及其组合概率，测量角度列表
    theory_result=[np.trace(np.dot(rouxyz,a)) for a in Aq_measure_list]
    mea_num=state_numb*len(Aq_measure_list)
    result=[]#记录27个测量的三组最终结果
    para=[RabiTime1,RabiTime2,freq1,freq2]#打包参数
    all_prob=[[],[],[]]#好像还需要累积所有100次的每次的结果。。
    all_eff=[]
    all_count=[]
    N_independent_results=[]#记录N次独立测量结果，用于求errorbar
    for i in range(len(Aq_measure_list)):
        #print i
        #print '----------------------------------'
        for j in range(state_numb):
            a=Aq_pulse_list[i][j]
            b=Aq_eigenvecter_list[i][j]
            #print mea_angle_check(a,b)
        #print Aq_eigenvecter_list[i][2]
        
    for i in range(state_numb):
        prepareList=[(freq1*2*pi,pre_angle_list[i][2]*RabiTime1/2.0,pi*pre_angle_list[i][3]),
        (freq2*2*pi,pre_angle_list[i][0]*RabiTime2/2.0,pi*pre_angle_list[i][1]),(freq1*2*pi,pre_angle_list[i][4]*RabiTime1/2.0,pi*pre_angle_list[i][5]),] # Note the 2*pi
        prob,i_prob,N_temp,eff_of_i,count_of_i=measure(prepareList,para)
        print len(eff_of_i[0])
        N_independent_results.append(N_temp)#将N组k个独立测量值记录下来
        all_prob[i].append(i_prob)
        all_eff.append(eff_of_i)
        all_count.append(count_of_i)
        result.append(prob)
    #写两个函数，一个是将all_eff,all_count中的最后一项与前9项的值分开。得到的结果代入第二个函数
    #第二个函数，将输入的qutrit1，和qutrit2的结果处理得到最终的结果。
    #def get2result(all_eff,all_count):
        
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
    print average2
    value_of_Aq=getmeanvalue(Aq_prob_list,mixprobability,result)
    print '------9 Aq values of x=%f, y=%f--------------'%(x_point,y_point)
    print value_of_Aq
    print '----------theory_result---------------'
    print np.real(theory_result)
    triggerstop()
    plt.errorbar(range(len(Aq_measure_list)),average2,yerr=errorbar,fmt='-*',color='r')
    plt.show()
    

    '''
    print '----------Aq_prob_list-------------------------'
    print Aq_prob_list
    print '----------Aq_eigenvecter_list-------------------------'
    print Aq_eigenvecter_list
    print '----------Aq_pulse_list-------------------------'
    print Aq_pulse_list
    '''
    
 