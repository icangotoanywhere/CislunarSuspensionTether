# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 21:58:05 2024

@author: Www
"""

# -*- coding: utf-8 -*-

"""
Created on Mon May 10 23:36:13 2021
@author: 汪巍巍（wang weiwei)
"""
import time
import datetime
print('启动时间',datetime.datetime.now())
start=time.time()


# %%这是一个关于月球天梯动力学的程序
#%%%开始，加载各种东西
import numpy as np
import math
#import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib
matplotlib.rc("font", family='Microsoft YaHei')
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
from scipy.integrate import solve_ivp
from scipy.optimize import root
from scipy.fft import fft, fftfreq
#import bisect
#%%%分段数--------------------
#N=300#质点总数量,段数，月表为0
#name='右偏0.07-椭圆-粘性-10-无for-轿厢爬升'
name='例6.1'
#--------------------------
#这里将Setting文件的参数读取出来并赋值
f=open('Setting-算'+name+'.txt',encoding='utf-8')
f1=f.read().split('\n')
keys = {}
for line in f1[1:]:
    if line:
        key, value = line.split('=')
        keys[key] = value
f.close()
N=int(keys['N'])#质点总数量,段数，月表为0
c=float(keys['c']) #阻尼系数
x_end=float(keys['x_end'])  #天梯端点所处的位置，单位m,放到前面参数表去了
theta_0=float(keys['theta_0'])  #初始摆角
v_0=float(keys['v_0']) #有没有初速度的开关
epsilon=float(keys['epsilon']) #月球公转偏心率
print(f1)
# N=345#质点总数量,段数，月表为0,放到前面参数表去了
num =np.arange(0,N+1,1) #质点序号数组
print (num)
#%%%基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m**3*s**(-2)
mu_m = 4.9028E12 #月球引力常数，单位m**3*s**(-2)
omega = 2*math.pi/(27.32*24*3600) #月球公转角速度，单位S**(-1)
#epsilon=0.0549  #月球公转偏心率
a=3.84403E8 #月心距离地心的距离，半长轴，单位m
h=np.sqrt(a*mu_e*(1-epsilon**2))   #月心角动量

#--这几个几何参数仅用于计算绳子的参数--
x_e=-4800000 #地心，单位m
y_e=0
x_m = 3.84403E8 #月心距离地心的距离，单位m
y_m=0
x_ms = 3.7826E8+4403000 #月面距离地心的距离，单位m
# paras = [mu_e,mu_m,omega,x_e,x_m,x_ms]
# print(paras)
#x_end= 0.5E8+4000000 #天梯端点所处的位置，单位m
l_step=(x_ms-x_end)/N #每段天梯的长度，单位m


#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m**(-3)
S = 6000000000 #材料抗拉强度，单位N*m**(-2)，即Pa
#c=1e-10#阻尼系数
E=9E11 #材料的弹性模量
k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms,x_end-l_step,-l_step)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)

#%%%以下是研究变截面的情况-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式
    return 0.0002*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
A=np.ones(N+1)
for i in num:
    A[i]=Area(xi[i])  #所有节点的截面积

#%%%每个点的质量----------
def mass2(x): #定义某一段的质量
    return sp.integrate.quad(Area,x,x+l_step)[0]*rho #使用scipy的积分
#for i in xi[1:total]: print(mass2(i))
mass=np.zeros(N+1)
for i in num[1:N+1]:
    mass[i]=mass2(xi[i]) #python切片的规则很奇怪，要注意
#print("质量阵是这样的",mass)
#%%%每段的胡克系数-------------
def Areafan(x): #定义x处的横截面积Area的倒数，用于求胡克系数
    return Area(x)**(-1)
for i in num[1:N+1]:
    #k[i]=E*Area((xi[i-1]+xi[i])/2)/l_step #相邻两端中点值，k[0]是总平均，不过用不上
    k[i] = E/(sp.integrate.quad(Areafan,xi[i],xi[i-1])[0])#这是用积分的方法求胡克系数

#----轿厢的参数-----
vc=150  #轿厢沿绳的速度 单位 m/s
#mc=1000000 #轿厢的质量 单位 kg
mc=float(keys['mc']) #轿厢的质量 单位 kg

N1=np.argmax(A)
N2=int(np.floor(3*N/4))
N3=N
x00=x_ms
x10=xi[N1]
x20=xi[N2]
x30=x_end
M1=sp.integrate.quad(Area,x10,x00)[0]*rho
M2=sp.integrate.quad(Area,x20,x10)[0]*rho
M3=sp.integrate.quad(Area,x30,x20)[0]*rho
m1=3500000
m2=1200000
m3=900000
k1= E/(sp.integrate.quad(Areafan,x10,x00)[0])
k2= E/(sp.integrate.quad(Areafan,x20,x10)[0])
k3= E/(sp.integrate.quad(Areafan,x30,x20)[0])
l1=x00-x10
l2=x10-x20
l3=x20-x30

#前方高能预警，本程序的核心科技即将来临

#----确定地月参数的位置，真近点角和时间的关系---
#---开普勒方程，EE(t）函数
def EE(t):
    def kepler(E):  #开普勒方程，用于椭圆轨道
        return E-epsilon*np.sin(E)-omega*t
    ang_E = root(kepler, omega*t,tol=1E-24).x[0]
    return ang_E  #直接用scipy求解开普勒方程



#%%以下是定义欧拉-拉格朗日方程，即微分方程组------------------
def EL(t, q):
    Q=np.empty(shape=(2*3))#

    ang_E=EE(t)

    xm=x_m*(1-epsilon*np.cos(ang_E))
    xe=x_e*(1-epsilon*np.cos(ang_E+np.pi))
    x0=xm-1740000

    x1=q[0]
    x2=q[1]
    x3=q[2]
    v1=q[3]
    v2=q[4]
    v3=q[5]

    dl1=np.maximum(0,x0-x1-l1)
    dl2=np.maximum(0,x1-x2-l2)
    dl3=np.maximum(0,x2-x3-l3)

    ro1=M1/(x0-x1)
    ro2=M2/(x1-x2)
    ro3=M3/(x2-x3)

    cg = np.mat([[m1 + M1 / 3 + M2 / 3,M2 / 6,0],[M2 / 6,m2 + M2 / 3 + M3 / 3,M3 / 6],[0,M3 / 6,m3 + M3 / 3]])

    F1=ro1*mu_e*1/(-x1+xe)-ro1*mu_m*1/(xm-x1)+ro2*mu_e*1/(x1-xe)+ro2*mu_m*1/(xm-x1)+m1*mu_m*1/(xm-x1)**2-m1*mu_e*1/(x1-xe)**2+k1*dl1-k2*dl2+((2*M1+2*M2+6*m1)*x1+x0*M1+x2*M2)*omega**2/6
    F2=ro2*mu_e*1/(-x2+xe)-ro2*mu_m*1/(xm-x2)+ro3*mu_e*1/(x2-xe)+ro3*mu_m*1/(xm-x2)+m2*mu_m*1/(xm-x2)**2-m2*mu_e*1/(x2-xe)**2+k2*dl2-k3*dl3+((2*M2+2*M3+6*m2)*x2+x1*M2+M3*x3)*omega**2/6
    F3=ro3*mu_e*1/(-x3+xe)-ro3*mu_m*1/(xm-x3)+m3*mu_m*1/(xm-x3)**2-m3*mu_e*1/(x3-xe)**2+k3*dl3+omega**2*((2*M3+6*m3)*x3+M3*x2)/6

    #F1=M1*mu_e*np.log((x0-xe)/(x1-xe))/(x0-x1)**2+M1*mu_m*np.log((xm-x1)/(xm-x0))/(x0-x1)**2+M1*mu_e/(x0-x1)/(-x1+xe)-M1*mu_m/(x0-x1)/(xm-x1)+M2*mu_e*np.log((x2-xe)/(x1-xe))/(x1-x2)**2+M2*mu_m*np.log((-xm+x1)/(-xm+x2))/(x1-x2)**2+M2*mu_e/(x1-x2)/(x1-xe)+M2*mu_m/(x1-x2)/(xm-x1)+m1*mu_m/(xm-x1)**2-m1*mu_e/(x1-xe)**2+k1*dl1-k2*dl2
    #F2=-M2*mu_e*np.log((x2-xe)/(x1-xe))/(x1-x2)**2-M2*mu_m*np.log((-xm+x1)/(-xm+x2))/(x1-x2)**2+M2*mu_e/(x1-x2)/(-x2+xe)-M2*mu_m/(x1-x2)/(xm-x2)+M3*mu_e*np.log((x3-xe)/(x2-xe))/(x2-x3)**2+M3*mu_m*np.log((-xm+x2)/(-xm+x3))/(x2-x3)**2+M3*mu_e/(x2-x3)/(x2-xe)+M3*mu_m/(x2-x3)/(xm-x2)+m2*mu_m/(xm-x2)**2-m2*mu_e/(x2-xe)**2+k2*dl2-k3*dl3
    #F3=-M3*mu_e*np.log((x3-xe)/(x2-xe))/(x2-x3)**2-M3*mu_m*np.log((-xm+x2)/(-xm+x3))/(x2-x3)**2+M3*mu_e/(x2-x3)/(-x3+xe)-M3*mu_m/(x2-x3)/(xm-x3)+m3*mu_m/(xm-x3)**2-m3*mu_e/(x3-xe)**2+k3*dl3


    b = np.array([F1, F2, F3])  # b代表常数列
    #print(b)
    #print(dl1,dl2,dl3)
    jie = np.linalg.solve(cg, b)


    a1=jie[0]
    a2=jie[1]
    a3=jie[2]


    Q[0]=v1
    Q[1]=v2
    Q[2]=v3
    Q[3]=a1
    Q[4]=a2
    Q[5]=a3


    return Q


#%%定义完毕-下面给出初始条件------------------------------------------

q0=[x10,x20,x30,0,0,0]
#print("初始条件是",q0)
#print(EL(1,q0))
#%%龙格库塔求解-------------------------------------------------
tianshu=60 #仿真天数
t_span=[0, tianshu*24*3600]
t_eval=range(0,int(t_span[1]),6)#600代表10分钟，每1分钟输出一次结果
sol = solve_ivp(EL, t_span, q0, method='DOP853', t_eval=t_eval, rtol=1e-3,atol=1e-6,first_step=6)
print("输出步数为",len(sol.t))
#print(sol.y)

x_result=sol.y[0:3]
v_result=sol.y[3:]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

matplotlib.rc("font", family='Microsoft YaHei')

#%%绘制截图-圆形----------
fig,ax = plt.subplots(figsize=(9,4.5))

ax.plot(sol.t,x_result[0],'k',linewidth=0.8)
ax.plot(sol.t,x_result[1],'k',linewidth=0.8)
ax.plot(sol.t,x_result[2],'k',linewidth=0.8)

ax.text(0, x10*1.05, '$\mathdefault{x_1}$',fontsize=15, ha="center", va="center")
ax.text(0, x20*1.08, '$\mathdefault{x_2}$',fontsize=15, ha="center", va="center")
ax.text(0, x30*1.1, '$\mathdefault{x_3}$',fontsize=15, ha="center", va="center")

tline=np.arange(0, tianshu+1, 1)*24*3600
ymoon=np.zeros(tianshu+1)
yearth=np.zeros(tianshu+1)
for i in np.arange(0, tianshu+1, 1):
    ymoon[i]=x_m*(1-epsilon*np.cos(EE(tline[i])))
    yearth[i]=x_e*(1-epsilon*np.cos(EE(tline[i])+np.pi))

ax.plot(tline,ymoon,'grey',linewidth=2)
ax.plot(tline,yearth,'b',linewidth=4)

ax.text(0, x_m*1, 'Moon',fontsize=12, ha="center", va="center")
ax.text(0, 0.9E7, 'Earth',fontsize=12, ha="center", va="center")


def day(x):
   return x/24/3600
secax = ax.secondary_xaxis('top', functions=(day, day))
secax.set_xlabel('time (day)')
secax.set_xticks(np.arange(0, tianshu+2, 2))
ax.set_xlabel('time (s)')
ax.set_ylabel('x (m)')
#ax.set_title('Earth-Moon Planer',y=-0.12)
plt.show
plt.savefig('坐标时程'+'-'+name+'.svg')


#%%静态的数据图，典型点的傅里叶变换-------------------
fig, ax = plt.subplots(1,1,figsize=(9,4.5))
ax=plt.subplot(1,1,1)
Nt=len(sol.t)
Nn=100
x1_fft = 2.0/Nt*np.abs(fft(x_result[0])[1:Nn])
x2_fft = 2.0/Nt*np.abs(fft(x_result[1])[1:Nn])
x3_fft = 2.0/Nt*np.abs(fft(x_result[2])[1:Nn])
freq = fftfreq(Nt, sol.t[1] - sol.t[0])[1:Nn]
wid=0.00000008
#plt.yscale('log')
ax.bar(freq, x1_fft, width=wid, label= '$\mathdefault{x_1}$')
ax.bar(freq+wid, x2_fft, width=wid, label= '$\mathdefault{x_2}$')
ax.bar(freq+2*wid, x3_fft, width=wid, label= '$\mathdefault{x_3}$')

plt.axvline(x=4.4e-7, color="black", linestyle="--")
plt.axvline(x=6.65e-6, color="black", linestyle="--")
plt.axvline(x=8.25e-6, color="black", linestyle="--")
plt.axvline(x=1e-5, color="black", linestyle="--")
day1=1/6.65e-6/24/3600
day2=1/8.25e-6/24/3600
day3=1/10e-6/24/3600

#ax.set_xlim(0.00000001,0.00001)
#ax.set_ylim(0.0000001,fft_up)
ax.set_xlabel('$frequency (Hz)$') # the horizontal axis represents the time
ax.set_ylabel('Amplitude (m)')
#plt.axvline(x=1.75e-5, color="black", linestyle="--")
#plt.axvline(x=2.237e-5, color="black", linestyle="--")
def invert(x):
    # 1/x with special treatment of x == 0
    x = np.array(x).astype(float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 1 / x[~near_zero]
    return x/24/3600

secax = ax.secondary_xaxis('top', functions=(invert, invert))
secax.set_xlabel('Period (day)')
secax.set_xticks([0.6,1,day3,day2,day1,2,27.3])

plt.legend(ncol=1,labelspacing=0.1,columnspacing=0.4)
plt.tight_layout()
plt.show
plt.savefig('..\图像结果\算'+name+'-'+'傅里叶.svg',bbox_inches='tight')





end=time.time()
print('Running time: %s Seconds'%(end-start))
# end
