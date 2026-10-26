

# -*- coding: utf-8 -*-
"""
Created on 2021-08-17 17:56
@author: 汪巍巍（wang weiwei)
"""
# 这是一个关于月球天梯动力学的程序
# import os
# os.environ['CUDA_VISIBLE_DEVICES'] = '1'
import numpy as np
import math
#import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
#------------------------
N=int(5) # 质点总数量,段数，月表为0
#--------------------------
num =np.arange(0,N+1,1) #质点序号数组
#for i  in num: print(num[i])
print ("num=",num)
#-----基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m^3*s^(-2)
mu_m = 4.9028E12 #月球引力常数，单位m^3*s^(-2)
omega = 2*math.pi/(27.32*24*3600) #月球公转角速度，单位S^(-1)
x_e = -4800000 #地心距离地月公共质心的距离，单位m
y_e=0
x_m = 3.80E8 #月心距离地月公共给质心的距离，单位m
y_m=0
x_ms = 3.7826E8 #月面距离地月公共质心的距离，单位m
# paras = [mu_e,mu_m,omega,x_e,x_m,x_ms]
# print(paras) 
x_end= 0.4E8 #天梯端点所处的位置，单位m
l_step=(x_ms-x_end)/N #每段天梯的长度，单位m
#-----------材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
c=5e-10 #阻尼系数
K=3E10 #材料的弹性模量
KK=np.ones(N+1)*K
k=np.ones(N+1)*K
#k=np.zeros(N+1)
#----------------------
xi=np.arange(x_end,x_ms+l_step,l_step)  #各节点的位置
#print("xi 这个数组是这样的 ", xi)
#xi这个数组，我们这里是从绳端开始的，后面都要倒过来，要注意

#---以下是研究变截面的情况-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式
    return 0.01*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
#for i in xi: print(i,"处的截面积为",Area(i))

#-----每个点的质量----------
def mass2(x): #定义某一段的质量
    return sp.integrate.quad(Area,x-l_step,x)[0]*rho #使用scipy的积分
#for i in xi[1:total]: print(mass2(i))
mass=[0]
for i in num[1:N+1]: 
    mass.append(mass2(xi[-i])) #python切片的规则很奇怪，要注意

#print("质量阵是这样的",mass[-N])#因为我们是正着积分，所以这个数组要倒过来
#------每段的胡克系数-------------
for i in num: k[i]=K*Area((xi[-i-1]+xi[-i])/2)


#-----五力平衡非线性方程组-------
def F(xy0): #定义五力平衡非线性代数方程组，()中是未知量.(这个自变量不包括x0,y0)
    xy01=[x_ms,y_m]
    xy=np.append(np.append(xy01[0],xy0[0:N]),np.append(xy01[1],xy0[N:]))#（这里还要把x0,y0贴回来）
    FX=np.zeros(N+1)
    FY=np.zeros(N+1)
    x=np.zeros(N+1)
    y=np.zeros(N+1)
    de=np.zeros(N+1)
    de2=np.zeros(N+1)
    dm=np.zeros(N+1)
    dm2=np.zeros(N+1)
    l_new=np.zeros(N+1)
    lr=np.zeros(N+1)
    for i in np.arange(N+1):
      x[i]=xy[i]#一维拆成二维
      y[i]=xy[i+N+1]#一维拆成二维
      de2[i]=(x[i]-x_e)**2+(y[i]-y_e)**2#到地球距离平方
      de[i]=np.sqrt(de2[i])#到地球距离
      dm2[i]=(x[i]-x_m)**2+(y[i]-y_m)**2#到月球距离平方
      dm[i]=np.sqrt(dm2[i])#到月球距离
      l_new[i]=np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)#新杆长
      lr[i]=max(0,(l_new[i]-l_step)/l_step)#刚度关
      #lr[i]=(l_new[i]-l_step)/l_step#刚度开
    for i in np.arange(1,N):
        FX[i]=-mu_e*mass[i]/de2[i]*(x[i]-x_e)/de[i]\
              -mu_m*mass[i]/dm2[i]*(x[i]-x_m)/dm[i]\
              +mass[i]*x[i]*omega**2\
              -k[i]*lr[i]*(x[i]-x[i-1])/l_new[i]\
              -k[i+1]*lr[i+1]*(x[i]-x[i+1])/l_new[i+1]
        FY[i]=-mu_e*mass[i]/de2[i]*(y[i]-y_e)/de[i]\
              -mu_m*mass[i]/dm2[i]*(y[i]-y_m)/dm[i]\
              +mass[i]*y[i]*omega**2\
              -k[i]*lr[i]*(y[i]-y[i-1])/l_new[i]\
              -k[i+1]*lr[i+1]*(y[i]-y[i+1])/l_new[i+1]
    FX[N]=-mu_e*mass[N]/de2[N]*(x[N]-x_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(x[N]-x_m)/dm[N]\
          +mass[N]*x[N]*omega**2\
          -k[N]*(x[N]-x[N-1])*lr[N]/l_new[N]
    FY[N]=-mu_e*mass[N]/de2[N]*(y[N]-y_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(y[N]-y_m)/dm[N]\
          +mass[N]*y[N]*omega**2\
          -k[N]*(y[N]-y[N-1])*lr[N]/l_new[N]
    FXY=np.append(FX[1:],FY[1:])
                     
    print("x方向的合力",FX)  
    print("y方向的合力",FY)
    return FXY  #输出的是数组

#------定义迭代方程----------------
def PHI(xy0): #定义五力平衡非线性代数方程组，()中是未知量.(这个自变量不包括x0,y0)
    xy01=[x_ms,y_m]
    xy=np.append(np.append(xy01[0],xy0[0:N]),np.append(xy01[1],xy0[N:]))#（这里还要把x0,y0贴回来）
    FX=np.zeros(N+1)
    FY=np.zeros(N+1)
    x=np.zeros(N+1)
    y=np.zeros(N+1)
    de=np.zeros(N+1)
    de2=np.zeros(N+1)
    dm=np.zeros(N+1)
    dm2=np.zeros(N+1)
    l_new=np.zeros(N+1)
    lr=np.zeros(N+1)
    for i in np.arange(N+1):
      x[i]=xy[i]#一维拆成二维
      y[i]=xy[i+N+1]#一维拆成二维
      de2[i]=(x[i]-x_e)**2+(y[i]-y_e)**2#到地球距离平方
      de[i]=np.sqrt(de2[i])#到地球距离
      dm2[i]=(x[i]-x_m)**2+(y[i]-y_m)**2#到月球距离平方
      dm[i]=np.sqrt(dm2[i])#到月球距离
      l_new[i]=np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2)#新杆长
      lr[i]=max(0,(l_new[i]-l_step)/l_step)#刚度关
      #lr[i]=(l_new[i]-l_step)/l_step#刚度开
    for i in np.arange(1,N):
        FX[i]=-mu_e*mass[i]/de2[i]*(x[i]-x_e)/de[i]\
              -mu_m*mass[i]/dm2[i]*(x[i]-x_m)/dm[i]\
              +mass[i]*x[i]*omega**2\
              -k[i]*lr[i]*(x[i]-x[i-1])/l_new[i]\
              -k[i+1]*lr[i+1]*(x[i]-x[i+1])/l_new[i+1]
        FY[i]=-mu_e*mass[i]/de2[i]*(y[i]-y_e)/de[i]\
              -mu_m*mass[i]/dm2[i]*(y[i]-y_m)/dm[i]\
              +mass[i]*y[i]*omega**2\
              -k[i]*lr[i]*(y[i]-y[i-1])/l_new[i]\
              -k[i+1]*lr[i+1]*(y[i]-y[i+1])/l_new[i+1]
    FX[N]=-mu_e*mass[N]/de2[N]*(x[N]-x_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(x[N]-x_m)/dm[N]\
          +mass[N]*x[N]*omega**2\
          -k[N]*(x[N]-x[N-1])*lr[i]/l_new[N]
    FY[N]=-mu_e*mass[N]/de2[N]*(y[N]-y_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(y[N]-y_m)/dm[N]\
          +mass[N]*y[N]*omega**2\
          -k[N]*(y[N]-y[N-1])*lr[i]/l_new[N]
    for i in np.arange(1,N+1):      
          x[i]=(FX[i]+mu_e*mass[i]/de2[i]*(x[i]-x_e)/de[i])/(mu_e*mass[i]/de2[i]/de[i])+x_e
          y[i]=(FY[i]+mu_e*mass[i]/de2[i]*(y[i]-y_e)/de[i])/(mu_e*mass[i]/de2[i]/de[i])+y_e
          
    FXY=np.append(x[1:],y[1:])
    print(FXY)
    return FXY  #输出的是数组





#-------给出初值------------
theta0=0.1
x0=((xi[-1:-N-2:-1]-xi[N])*math.cos(theta0)+xi[N])[1:]
y0=((xi[N]-xi[-1:-N-2:-1])*math.sin(theta0))[1:]




print(np.arange(N+1))

PHI(PHI(np.append(x0,y0)))
print('每段长',l_step)
print('x初始值',x0)
print('y初始值',y0)
print('EA数组',k)
print('质量阵',mass)

'''
import scipy.optimize as SO
from scipy.optimize import root
from scipy.optimize import least_squares as lts
from scipy.optimize import minimize
resl=root(F,np.append(x0,y0),method='lm')
#resl=SO.excitingmixing(F,np.append(x0,y0))
print('解是',resl)
print('每段长',l_step)
print('x初始值',x0)
print('y初始值',y0)
print('EA数组',k)
print('质量阵',mass)
'''








