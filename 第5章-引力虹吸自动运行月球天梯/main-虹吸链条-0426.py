# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 21:13:05 2024

@author: 汪巍巍
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
#import bisect
#%%%分段数--------------------

name='例1'

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
c_mass_p=float(keys['c_mass_p']) #货物占绳和空筐质量的比例
speed_initial=float(keys['speed_initial']) #初速度，单位m/s
print(f1)
# N=345#质点总数量,段数，月表为0,放到前面参数表去了
num =np.arange(0,N+1,1) #质点序号数组
print (num)
#%%%基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m^3*s^(-2)
mu_m = 4.9028E12 #月球引力常数，单位m^3*s^(-2)
omega = -2*math.pi/(27.32*24*3600) #月球公转角速度，单位S^(-1)
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
#x_end= 0.5E8+4000000 #天梯端点所处的位置，单位m
ls=(x_ms-x_end)/N #每段天梯的长度，单位m


#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
#c=1e-10#阻尼系数
E=8E11 #材料的弹性模量
k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms-ls/2,x_end,-ls)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)

#%%%截面积-----------------
Area=0.0001 #平方米，1厘米粗的
#%%%每一段加上筐的质量----------
mass_b=ls*Area*rho
mass_p=c_mass_p*mass_b
#mass=ls*Area*rho
#%%%每段的胡克系数-------------
k=E*Area/ls

M=2*N



#前方高能预警，本程序的核心科技即将来临
#%%以下是定义欧拉-拉格朗日方程，即微分方程组------------------
def YunDong(t, q):
    Q=np.zeros(shape=4*(M+1))#定义Q
    #---地心和月心位置
    xm=x_m
    ym=y_m
    xe=0 #这几个变量可以是t的函数
    ye=0 #这几个变量可以是t的函数
    #----一维拆成四维----
    x=q[0:M] #
    y=q[M:2*M] #
    vx=q[2*M:3*M] #
    vy=q[3*M:4*M] #


    xce=q[-4]
    yce=q[-3]
    vxce=q[-2]
    vyce=q[-1]



    xcm=x_ms
    ycm=0

    #----算出各段几何参数----
    le=np.sqrt((x-xe)**2+(y-ye)**2)  #到地心距离
    lm=np.sqrt((x-xm)**2+(y-ym)**2)  #到月心距离
    difx=np.diff(np.append(x,x[0]))
    dify=np.diff(np.append(y,y[0]))
    #difrx=-np.diff(np.append(x[-1],x))
    #difry=-np.diff(np.append(y[-1],y))
    ll=np.sqrt(difx**2+dify**2) #杆长
    lr=np.maximum(0,1-ls/ll)  #伸长率，柔性体现在这里，绳子没有推力，只能拉力-体现在max函数中



    ############################################################
    #接下来这部分是找哪几个点和滑轮有接触，根据x的位置判断
    dist_ce=np.sqrt((x-xce)**2+(y-yce)**2)
    mininum_ce=np.argmin(dist_ce) #距离动滑轮最小近的点的序号
    nceb=mininum_ce-1 #最近点的前一个点序号
    ncea=mininum_ce+1 #最近点的后一个点序号
    if ncea==M :
        ncea=0
    if dist_ce[ncea]<dist_ce[nceb]:
        neb=mininum_ce
        nea=ncea
    else:
        neb=nceb
        nea=mininum_ce
    neaa=nea+1
    if neaa==M:
        neaa=0
    if neb==-1:
        neb=M-1

    dist_cm=np.sqrt((xcm-x)**2+(ycm-y)**2)
    mininum_cm=np.argmin(dist_cm) #距离动滑轮最小近的点的序号
    ncmb=mininum_cm-1 #最近点的前一个点序号
    ncma=mininum_cm+1 #最近点的后一个点序号
    if ncma==M :
        ncma=0
    if dist_cm[ncma]<dist_cm[ncmb]:
        nmb=mininum_cm
        nma=ncma
    else:
        nmb=ncmb
        nma=mininum_cm
    nmaa=nma+1
    if nmaa==M:
        nmaa=0
    if nmb==-1:
        nmb=M-1
    #############################################################

    #print(nmb,nma,neb,nea)

    #np.clip(x,xce,xcm,out=x)
    #if x[]
    #vx[nea]=vx[neaa]
    #vy[nea]=vy[neaa]
    #vx[nea]=(-vx[neb]+vx[neaa])/2
    #vy[nea]=(vy[neb]+vy[neaa])/2
    vx[nea]=np.sqrt(vx[nea]**2+vy[nea]**2)*(x[neaa]-xce)/dist_ce[neaa]
    vy[nea]=np.sqrt(vx[nea]**2+vy[nea]**2)*(y[neaa]-yce)/dist_ce[neaa]

    #vx[nma]=vx[nmaa]
    #vy[nma]=vy[nmaa]
    #vx[nma]=(-vx[nmb]+vx[nmaa])/2
    #vy[nma]=(vy[nmb]+vy[nmaa])/2
    vx[nma]=np.sqrt(vx[nma]**2+vy[nma]**2)*(x[nmaa]-xcm)/dist_cm[nmaa]
    vy[nma]=np.sqrt(vx[nma]**2+vy[nma]**2)*(y[nmaa]-ycm)/dist_cm[nmaa]

    lce=dist_ce[nea]+dist_ce[neb]
    dle=np.maximum(lce-ls,0)
    lcm=dist_cm[nma]+dist_cm[nmb]
    dlm=np.maximum(lcm-ls,0)

    mass=np.ones(shape=M)*mass_b

    if nmaa<neb:
        mass[nmaa:neb]=mass_b+mass_p
        mass[neb]=mass_b+mass_p
    elif nmaa>neb:
        mass[nmaa:]=mass_b+mass_p
        mass[:neb]=mass_b+mass_p
        mass[neb]=mass_b+mass_p



    massce=mass_b*N/3


    ax=-mu_m*(x-xm)/(lm**3)\
       -mu_e*(x-xe)/(le**3)\
       +k*difx*lr/mass\
       +k*(-np.roll(difx,1))*(np.roll(lr,1))/mass\
       -c*vx/mass\
       +(omega**2)*x\
       +2*omega*vy
    ay=-mu_m*(y-ym)/(lm**3)\
       -mu_e*(y-ye)/(le**3)\
       +k*dify*lr/mass\
       +k*(-np.roll(dify,1))*(np.roll(lr,1))/mass\
       -c*vy/mass\
       +(omega**2)*y\
       -2*omega*vx

   ########################################################################
    ax[nea]=-mu_m*(x[nea]-xm)/(lm[nea]**3)\
            -mu_e*(x[nea]-xe)/(le[nea]**3)\
            +k*lr[nea]*(x[neaa]-x[nea])/mass[nea]\
            +k*dle*(xce-x[nea])/dist_ce[nea]/mass[nea]\
            -c*vx[nea]/mass[nea]\
            +x[nea]*(omega**2)\
            +2*omega*vy[nea]
    ax[neb]=-mu_m*(x[neb]-xm)/(lm[neb]**3)\
            -mu_e*(x[neb]-xe)/(le[neb]**3)\
            +k*lr[neb-1]*(x[neb-1]-x[neb])/mass[neb]\
            +k*dle*(xce-x[neb])/dist_ce[neb]/mass[neb]\
            -c*vx[neb]/mass[neb]\
            +x[neb]*(omega**2)\
            +2*omega*vy[neb]
    ay[nea]=-mu_m*(y[nea]-ym)/(lm[nea]**3)\
            -mu_e*(y[nea]-ye)/(le[nea]**3)\
            +k*lr[nea]*(y[neaa]-y[nea])/mass[nea]\
            +k*dle*(yce-y[nea])/dist_ce[nea]/mass[nea]\
            -c*vy[nea]/mass[nea]\
            +y[nea]*(omega**2)\
            -2*omega*vx[nea]
    ay[neb]=-mu_m*(y[neb]-ym)/(lm[neb]**3)\
            -mu_e*(y[neb]-ye)/(le[neb]**3)\
            +k*lr[neb-1]*(y[neb-1]-y[neb])/mass[neb]\
            +k*dle*(yce-y[neb])/dist_ce[neb]/mass[neb]\
            -c*vy[neb]/mass[neb]\
            +y[neb]*(omega**2)\
            -2*omega*vx[neb]


    ax[nma]=-mu_e*(x[nma]-xe)/(le[nma]**3)\
            -mu_m*(x[nma]-xm)/(lm[nma]**3)\
            +k*lr[nma]*(x[nmaa]-x[nma])/mass[nma]\
            +k*dlm*(xcm-x[nma])/dist_cm[nma]/mass[nma]\
            -c*vx[nma]/mass[nma]\
            +x[nma]*(omega**2)\
            +2*omega*vy[nma]
    ax[nmb]=-mu_e*(x[nmb]-xe)/(le[nmb]**3)\
            -mu_m*(x[nmb]-xm)/(lm[nmb]**3)\
            +k*lr[nmb-1]*(x[nmb-1]-x[nmb])/mass[nmb]\
            +k*dlm*(xcm-x[nmb])/dist_cm[nmb]/mass[nmb]\
            -c*vx[nmb]/mass[nmb]\
            +x[nmb]*(omega**2)\
            +2*omega*vy[nmb]
    ay[nma]=-mu_e*(y[nma]-ye)/(le[nma]**3)\
            -mu_m*(y[nma]-ym)/(lm[nma]**3)\
            +k*lr[nma]*(y[nmaa]-y[nma])/mass[nma]\
            +k*dlm*(ycm-y[nma])/dist_cm[nma]/mass[nma]\
            -c*vy[nma]/mass[nma]\
            +y[nma]*(omega**2)\
            -2*omega*vx[nma]
    ay[nmb]=-mu_e*(y[nmb]-ye)/(le[nmb]**3)\
            -mu_m*(y[nmb]-ym)/(lm[nmb]**3)\
            +k*lr[nmb-1]*(y[nmb-1]-y[nmb])/mass[nmb]\
            +k*dlm*(ycm-y[nmb])/dist_cm[nmb]/mass[nmb]\
            -c*vy[nmb]/mass[nmb]\
            +y[nmb]*(omega**2)\
            -2*omega*vx[nmb]

    axce=-mu_m*(xce-xm)/(np.sqrt((xce-xm)**2+(yce-ym)**2)**3)\
         -mu_e*(xce-xe)/(np.sqrt((xce-xe)**2+(yce-ye)**2)**3)\
         -k*(lce-ls)*(xce-x[nea])/dist_ce[nea]/massce\
         -k*(lce-ls)*(xce-x[neb])/dist_ce[neb]/massce\
         -c*vxce/massce\
         +xce*(omega**2)\
         +2*omega*vyce

    ayce=-mu_m*(yce-ym)/(np.sqrt((xce-xm)**2+(yce-ym)**2)**3)\
         -mu_e*(yce-ye)/(np.sqrt((xce-xe)**2+(yce-ye)**2)**3)\
         -k*(lce-ls)*(yce-y[nea])/dist_ce[nea]/massce\
         -k*(lce-ls)*(yce-y[neb])/dist_ce[neb]/massce\
         -c*vyce/massce\
         +yce*(omega**2)\
         -2*omega*vxce

#----求结果Q-----
    Q[0:M]=vx #降阶，2阶方程将为1阶,dx/dt=v(x)
    Q[M:2*M]=vy #降阶，2阶方程将为1阶,dy/dt=y(x)
    Q[2*M:3*M]=ax #
    Q[3*M:4*M]=ay #
    Q[-4]=vxce
    Q[-3]=vyce
    Q[-2]=axce
    Q[-1]=ayce

    return Q


#%%定义完毕-下面给出初始条件------------------------------------------
#theta_0=0.07  #初始摆角
pian=ls/200
x=np.append(xi-pian,xi[::-1]+pian)
y=x*0
#vx=np.zeros(N+1)
#vy=np.zeros(N+1)
vx=np.zeros(2*N)
vx[:N]=-1*speed_initial
vx[N:]=speed_initial
vy=np.zeros(2*N)
ce=[x_end,0,0,0]

q0=np.concatenate((x,y,vx,vy,ce),axis=0)
#print("初始条件是",q0)
#print(EL(1,q0))
#%%龙格库塔求解-------------------------------------------------
#sol=sp.integrate.RK45(EL, 0, q0, 10, max_step=0.1, rtol=1e-10, atol=1e-10)
t_span=[0, 27.3*24*3600]
t_eval=range(0,int(t_span[1]),600)
sol = solve_ivp(YunDong, t_span, q0, method='DOP853', t_eval=t_eval, rtol=1e-3,atol=1e-6,max_step=6)
print("输出步数为",len(sol.t))
#print(sol.y)
print("Q向量的行数为",sol.y.shape[0])  #shape[0]的意思是行数

x_result=sol.y[0:M]
y_result=sol.y[M:2*M]
vx_result=sol.y[2*M:3*M]
vy_result=sol.y[3*M:4*M]
ce_result=sol.y[4*M:]

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#%%输出结果数据--
np.savetxt('..\数值结果\算'+name+'-'+'outt.txt', sol.t, fmt="%d")
np.savetxt('..\数值结果\算'+name+'-'+'outx.txt', x_result)
np.savetxt('..\数值结果\算'+name+'-'+'outy.txt', y_result)
np.savetxt('..\数值结果\算'+name+'-'+'outvx.txt', vx_result)
np.savetxt('..\数值结果\算'+name+'-'+'outvy.txt', vy_result)
np.savetxt('..\数值结果\算'+name+'-'+'outce.txt', ce_result)

#%%绘制截图-展开2旋转----------
fig,ax = plt.subplots(figsize=(8,3))
#ax.set_xlim(-x_m*0.1,x_m*5)
ax.set_ylim(-x_m*0.05,x_m*1.2)

imgm = plt.imread('moon3.png')#读月球
imge = plt.imread('earth.png')#画地球
L3=range(0,x_result.shape[1],6*4)  #
Ld=range(0,x_result.shape[1],60*24) #1天

pianyi=np.arange(x_result.shape[1])*2e7*len(L3)/x_result.shape[1]
#---坐标转换---

u = y_result+pianyi
v = x_result
#uc = yc_result*np.cos(ang)-xc_result*np.sin(ang)+x_ms*ang
#vc = yc_result*np.sin(ang)+xc_result*np.cos(ang)


#--画线画月亮--

for i in L3:
   ax.plot(u[:,i],v[:,i],'k',linewidth=0.8)
   ax.plot([u[-1,i],u[0,i]],[v[-1,i],v[0,i]],color='g', linewidth=2)
   ax.scatter(u[:,i],v[:,i],c='k',s=3)
   ax.scatter(u[0,i],v[0,i],c='g',s=3)
   ax.scatter([pianyi[i],sol.y[-3,i]+pianyi[i]],[x_ms,sol.y[-4,i]],c='b',s=6)
   ax.imshow(imgm,extent=(pianyi[i]-4e6,pianyi[i]+4e6,x_m-4e6,x_m+4e6),aspect='equal')#画月亮
   ax.imshow(imge,extent=(pianyi[i]-6e6,pianyi[i]+6e6,-6e6,6e6),aspect='equal')#画地球
   #ax.scatter(uc[i],vc[i],c='k')
#--设时间标--
xticklables=[]
uu=[]
for i in Ld :
   xticklables.append('day'+str(int(i/144)))
   uu.append(u[0,i])
def ziji(x):
    return x
ax.set_xticks(uu,xticklables,minor=False,rotation=45)
secax = ax.secondary_xaxis('top', functions=(ziji, ziji))
secax.set_xlabel('u (m)')
secax.set_xticks(np.arange(0, 3, 0.1)*1e9)
#ax.set_xlabel('time')
ax.set_ylabel('v (m)')
plt.tight_layout()
plt.show
plt.savefig('..\图像结果\算'+name+'-'+'展开.svg',bbox_inches='tight')

end=time.time()
print('Running time: %s Seconds'%(end-start))


