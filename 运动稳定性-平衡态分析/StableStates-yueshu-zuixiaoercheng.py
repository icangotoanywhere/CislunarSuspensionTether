# -*- coding: utf-8 -*-
"""
Created on 2021-08-17 17:56
@author: 汪巍巍（wang weiwei)
"""
# 这是一个关于月球天梯平衡态的程序
# import os
# os.environ['CUDA_VISIBLE_DEVICES'] = '1'#试图使用GPU
import numpy as np
import math
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
#---以上为加载各种必备的工具包------
#---程序正式开始--------
N=3 # 绳子总段数,节点，月表为0，终点为N
#--------------------------
num =np.arange(0,N+1,1) #质点序号数组 
#for i  in num: print(num[i])
print ("num=",num)
#-----基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m^3*s^(-2)
mu_m = 4.9028E12 #月球引力常数，单位m^3*s^(-2)
omega = 2*math.pi/(27.32*24*3600) #月球公转角速度，单位S^(-1)
x_e = -4800000 #地心距离地月公共质心的距离（地心x坐标），单位m
y_e=0 #地心的y坐标
x_m = 3.80E8 #月心距离地月公共给质心的距离（月心x坐标），单位m
y_m=0 #月心y坐标
x_ms = 3.7826E8 #月面距离地月公共质心的距离（月心减去月球半径），单位m
x_end= 0.4E8 #天梯端点所处的位置，单位m
l_s=(x_ms-x_end)/N #每段天梯的长度，单位m
#-----------材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
c=5e-10 #阻尼系数,这个程序中没用到
E=3E10 #材料的弹性模量,单位N*m^(-2)
k=np.ones(N+1)*E #胡克系数初值，都设为1
#----------------------
xi=np.arange(x_end,x_ms+l_s,l_s)  #各节点的位置
#print("xi 这个数组是这样的 ", xi)
#xi这个数组，我们这里是从绳端开始的，后面都要倒过来，要注意

#---以下是研究变截面的情况，各处的截面积-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式（详见论文中的方程7）
    return 0.01*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
#for i in xi: print(i,"处的截面积为",Area(i))

#-- ---每个点的质量----  ------
def mass2(x): #定义某一段的质量
    return sp.integrate.quad(Area,x-l_s,x)[0]*rho #每段截面积积分乘密度
#for i in xi[1:total]: print(mass2(i))
mass=[0]
for i in num[1:N+1]: 
    mass.append(mass2(xi[-i])) #python切片的规则很奇怪，要注意
#print("质量阵是这样的",mass[-N])#因为我们是正着积分，所以这个数组要倒过来

#------每段的胡克系数-------------
for i in num: k[i]=E*Area((xi[-i-1]+xi[-i])/2)
#注意，是倒着来的，取得是中点的值
#如果要更精确，应该每段截面积积分除以段长
#所以段取短一些，这里更精确

#-----五力平衡非线性方程组-------
def F(xy0): #定义五力平衡非线性代数方程组，()中是未知量.(这个自变量不包括x0,y0)
    xy=np.append(np.append(x_ms,xy0[0:N]),np.append(y_m,xy0[N:]))
    #这是把第0个点贴回来，就是月面那个点
    FX=np.zeros(N+1) #初始化x方向合力
    FY=np.zeros(N+1) #初始化y方向合力
    x=np.zeros(N+1) #初始化
    y=np.zeros(N+1) #初始化
    de=np.zeros(N+1) #初始化到地球距离
    de2=np.zeros(N+1) #初始化到地球距离平方
    dm=np.zeros(N+1) #初始化到月球距离
    dm2=np.zeros(N+1) #初始化到月球距离平方
    l_new=np.zeros(N+1) #初始化各段新长度
    lr=np.zeros(N+1) #初始化各段伸长率
    for i in np.arange(N+1): #从0到N
      x[i]=xy[i]#一维拆成二维，前半段是x坐标，第0个是月面
      y[i]=xy[i+N+1]#一维拆成二维，后半段是y坐标，第0个是月面
      de2[i]=(x[i]-x_e)**2+(y[i]-y_e)**2#到地心距离平方
      de[i]=np.sqrt(de2[i])#到地心距离
      dm2[i]=(x[i]-x_m)**2+(y[i]-y_m)**2#到月心距离平方
      dm[i]=np.sqrt(dm2[i])#到月心距离
      l_new[i]=np.sqrt((x[i]-x[i-1])**2+(y[i]-y[i-1])**2) #新杆长
      lr[i]=max(0,(l_new[i]-l_s)/l_s)#刚度关
      #lr[i]=(l_new[i]-l_step)/l_s#刚度开
    for i in np.arange(1,N): #从1到N
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
    #最后一个点的写法不同，如下
    FX[N]=-mu_e*mass[N]/de2[N]*(x[N]-x_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(x[N]-x_m)/dm[N]\
          +mass[N]*x[N]*omega**2\
          -k[N]*lr[N]*(x[N]-x[N-1])/l_new[N]
    FY[N]=-mu_e*mass[N]/de2[N]*(y[N]-y_e)/de[N]\
          -mu_m*mass[N]/dm2[N]*(y[N]-y_m)/dm[N]\
          +mass[N]*y[N]*omega**2\
          -k[N]*lr[N]*(y[N]-y[N-1])/l_new[N]
    FXY=np.append(FX[1:],FY[1:]) #贴成一个向量，从第1个点开始，不管0点
    lstsqr=np.linalg.norm(FXY) #FXY的2-范数                 
    print("x方向的合力",FX)  
    print("y方向的合力",FY)
    print('合力数组',FXY)
    print('合力数组模',lstsqr)
    return lstsqr  #输出的是FXY的2-范数，一个标量
    # return FXY #输出的是FXY，一个向量

#-------给出初值------------
theta0=0.0 #初始给一个偏角，弧度
x0=((xi[-1:-N-2:-1]-xi[N])*math.cos(theta0)+xi[N])[1:] #初始各点的x坐标
y0=((xi[N]-xi[-1:-N-2:-1])*math.sin(theta0))[1:]  #初始各点的y坐标

#----走一个--------
F(np.append(x0,y0))

'''
#cons是约束条件
cons = ({'type': 'ineq', 'fun': lambda x: (x[0]-x_ms)**2+(x[0+N]-y_m)**2-l_s**2-0.1},
        {'type': 'ineq', 'fun': lambda x: (x[1]-x[1-1])**2+(x[1+N]-x[1+N-1])**2-l_s**2-0.1}
        ,{'type': 'ineq', 'fun': lambda x: (x[2]-x[2-1])**2+(x[2+N]-x[2+N-1])**2-l_s**2-0.1})
'''

import scipy.optimize as SO
from scipy.optimize import root
from scipy.optimize import least_squares as lts
from scipy.optimize import minimize

'''
resl=minimize(F,np.append(x0,y0),method='COBYLA',constraints=cons)
#resl=SO.excitingmixing(F,np.append(x0,y0))
print('解是',resl)
'''

print('每段长',l_s)
print('x初始值',x0)
print('y初始值',y0)
print('EA数组',k)
print('质量阵',mass)








