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
import os
import numpy as np
import math
#import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib
matplotlib.rc("font", family='Microsoft YaHei')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
from scipy.optimize import root


#%%%读取算例的参数并赋值--------------------

name='例8'
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
#--------------------------
num =np.arange(0,N+1,1) #质点序号数组
#print (num)
#%%%基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m^3*s^(-2)
mu_m = 4.9028E12 #月球引力常数，单位m^3*s^(-2)
omega = 2*math.pi/(27.32*24*3600) #月球公转角速度，单位S^(-1)
#epsilon=0.0549  #月球公转偏心率
a=3.84403E8 #月心距离地心的距离，半长轴，单位m
h=np.sqrt(a*mu_e*(1-epsilon**2))   #月心角动量

#--这几个几何参数仅用于计算绳子的参数--
x_e=0 #地心，单位m
y_e=0
x_m = 3.84403E8 #月心距离地心的距离，单位m
y_m=0
x_ms = 3.7826E8+4403000 #月面距离地心的距离，单位m
# paras = [mu_e,mu_m,omega,x_e,x_m,x_ms]
# print(paras)
#x_end= 0.5E8+4000000 #天梯端点所处的位置，单位m,放到前面参数表去了
ls=(x_ms-x_end)/N #每段天梯的长度，单位m


#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
#c=1e-10#阻尼系数
E=9E11 #材料的弹性模量
k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms,x_end-ls,-ls)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)

#%%%以下是研究变截面的情况-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式
    return 0.01*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
#for x in xi:
#    print(x,"处的截面积为",Area(x))

#%%%每个点的质量----------
def mass2(x): #定义某一段的质量
    return sp.integrate.quad(Area,x,x+ls)[0]*rho #使用scipy的积分
#for i in xi[1:total]: print(mass2(i))
mass=np.zeros(N+1)
for i in num[1:N+1]:
    mass[i]=mass2(xi[i]) #python切片的规则很奇怪，要注意
#print("质量阵是这样的",mass)
#%%%每段的胡克系数-------------
def Areafan(x): #定义x处的横截面积Area的倒数，用于求胡克系数
    return Area(x)**(-1)
for i in num[1:N+1]:
    #k[i]=E*Area((xi[i-1]+xi[i])/2)/ls #相邻两端中点值，k[0]是总平均，不过用不上
    k[i] = E/(sp.integrate.quad(Areafan,xi[i],xi[i-1])[0])#这是用积分的方法求胡克系数



#前方高能预警，本程序的核心科技即将来临
#%%以下是定义欧拉-拉格朗日方程，即微分方程组------------------
def EL(t, q):
    Q=np.empty(shape=(4*N+4))#4N+4维向量

    #----确定地月参数的位置---

    def kepler(E):  #开普勒方程，用于椭圆轨道
        return E-epsilon*np.sin(E)-omega*t
    ang_E = root(kepler, omega*t,tol=1E-24).x[0]  #直接用scipy求解开普勒方程
    theta=np.arctan(np.sqrt((1+epsilon)/(1-epsilon))*np.tan(ang_E/2))*2
    if (theta < 0):
        theta=2*np.pi+theta
    r=a*(1-epsilon**2)/(1+epsilon*np.cos(theta)) #月心模长
    theta_dot=mu_e*(1+epsilon*np.cos(theta))/(h*r)
    #---月心位置
    xm=r*np.cos(theta) #这几个变量是t的函数
    ym=r*np.sin(theta) #这几个变量是t的函数
    xe=0 #这几个变量可以是t的函数
    ye=0 #这几个变量可以是t的函数


    #----一维拆成四维----
    x=q[0:4*N+1:4] #每行第0列是x
    y=q[1:4*N+2:4] #每行第1列是y
    vx=q[2:4*N+3:4] #每行第2列是vx
    vy=q[3:4*N+4:4] #每行第3列是vy


    #----算出各段几何参数----
    le=np.sqrt((x-xe)**2+(y-ye)**2)  #到地心距离
    lm=np.sqrt((x-xm)**2+(y-ym)**2)  #到月心距离
    ll=np.append(0,np.sqrt(np.diff(x)**2+np.diff(y)**2))  #杆长
    lr=np.append(0,np.maximum(0,1-ls/ll[1:]))  #伸长率，柔性体现在这里，绳子没有推力，只能拉力-体现在max函数中

    #----求结果Q-----
    Q[0:4*N+1:4]=vx #降阶，2阶方程将为1阶,dx/dt=v(x)
    Q[1:4*N+2:4]=vy #降阶，2阶方程将为1阶,dy/dt=y(x)

    #----根（月球）的运动速度和加速度------
    Q[0]=-mu_e/h*np.sin(theta)*x_ms/a         #有根树系统的根（即运动已知的点）
    Q[1]=mu_e/h*(epsilon+np.cos(theta))*x_ms/a         #前四个数值确定，说明这是有根树系统的根（即运动已知的点）
    Q[2]=-mu_e/h*np.cos(theta)*theta_dot*x_ms/a    #有根树系统的根（即运动已知的点）
    Q[3]=-mu_e/h*np.sin(theta)*theta_dot*x_ms/a      #如果月球是静止的，速度为0，这里就都等于0

    #----求加速度，这里去掉了for循环---
    Q[4+2:4*N+2:4]=-mu_m*(x[1:N]-xm)/(lm[1:N]**3)\
                 -mu_e*(x[1:N]-xe)/(le[1:N]**3)\
                 -k[1:N]*np.diff(x[0:N])*lr[1:N]/mass[1:N]\
                 -k[2:]*(-np.diff(x[1:]))*lr[2:]/mass[1:N]\
                 -c*(vx[1:N]-y[1:N]*omega)
    Q[4+3:4*N+3:4]=-mu_m*(y[1:N]-ym)/(lm[1:N]**3)\
                 -mu_e*(y[1:N]-ye)/(le[1:N]**3)\
                 -k[1:N]*np.diff(y[0:N])*lr[1:N]/mass[1:N]\
                 -k[2:]*(-np.diff(y[1:]))*lr[2:]/mass[1:N]\
                 -c*(vy[1:N]+x[1:N]*omega)
    #特别要注意的是，最后一组自由度，弹性势能项只和前一组有关，与后一组无关了
    Q[4*N+2]=-mu_m*(x[N]-xm)/(lm[N]**3) \
             -mu_e*(x[N]-xe)/(le[N]**3) \
             -k[N]*(x[N]-x[N-1])*lr[N]/mass[N] \
             -c*(vx[N]-y[N]*omega)
    Q[4*N+3]=-mu_m*(y[N]-ym)/(lm[N]**3) \
             -mu_e*(y[N]-ye)/(le[N]**3) \
             -k[N]*(y[N]-y[N-1])*lr[N]/mass[N] \
             -c*(vy[N]+x[N]*omega)
    return Q


#%%定义完毕-下面给出初始条件------------------------------------------
#theta_0=0.07  #初始摆角
x=(xi-xi[0])*math.cos(theta_0)+x_ms*(1-epsilon)
y=(xi[0]-xi)*math.sin(theta_0)

#vx=np.zeros(N+1)
#vy=np.zeros(N+1)
vx=-mu_e/h*np.sin(theta_0)*y/a
vy=mu_e/h*(epsilon+np.cos(theta_0))*x/a
vy[N-2:]=vy[N-2:]*200
vx=vx*v_0
vy=vy*v_0
q0=np.transpose([x,y,vx,vy]).reshape(-1)
#print("初始条件是",q0)
#print(EL(1,q0))
#%%龙格库塔求解-------------------------------------------------
#sol=sp.integrate.RK45(EL, 0, q0, 10, max_step=0.1, rtol=1e-10, atol=1e-10)
t_span=[0, 27.3*24*3600]
t_eval=range(0,int(t_span[1]),600)
sol = solve_ivp(EL, t_span, q0, method='RK45', t_eval=t_eval, rtol=1e-3,atol=1e-6,first_step=5)
print("输出步数为",len(sol.t))
#print(sol.y)
print("Q向量的行数为",sol.y.shape[0])  #shape[0]的意思是行数

x_result=sol.y[0:N*4+1:4]
y_result=sol.y[1:N*4+2:4]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#%%输出结果数据--
np.savetxt('..\数值结果\算'+name+'-'+'outt.txt', sol.t, fmt="%d")
np.savetxt('..\数值结果\算'+name+'-'+'outx.txt', x_result)
np.savetxt('..\数值结果\算'+name+'-'+'outy.txt', y_result)

#%%绘制截图-圆形----------
fig,ax = plt.subplots(figsize=(10,10))
ax.set_xlim(-x_m*1.2,x_m*1.2)
ax.set_ylim(-x_m*1.2,x_m*1.2)

img = plt.imread('earth.png')#画地球
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球

L=range(0,x_result.shape[1] ,6*3)  #6代表1小时
for i in L:
   ax.plot(y_result[:,i],x_result[:,i],'k',linewidth=0.8)

imgm = plt.imread('moon3.png')#读月球
for i in L:
   ax.imshow(imgm,extent=(y_result[0,i]-4e6,y_result[0,i]+4e6,x_result[0,i]-4e6,x_result[0,i]+4e6),aspect='auto')#画月亮

#----
ha=["left","right"]
va=["bottom","top"]
for i in range(0,x_result.shape[1] ,6*24):  #6代表1小时
    ax.text(y_result[0,i],x_result[0,i],'t='+str(int(i/144))+'day',\
            ha=ha[bool(i > 14*144)],va=va[7*144 < i < 21*144],\
            size=10) #布尔值用来使文字出现在外围
#----
ax.set_xlabel('km')
ax.set_title('Earth-Moon Planer',y=-0.12)
plt.show

end=time.time()
print('Running time: %s Seconds'%(end-start))

