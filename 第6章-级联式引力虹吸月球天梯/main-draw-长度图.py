# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 00:17:20 2022

@author: Www
"""

# %%这是一个关于虹吸链条的程序
#%%%开始，加载各种东西
# import os
# os.environ['CUDA_VISIBLE_DEVICES'] = '1'
import numpy as np
import math
#import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib
matplotlib.rc("font", family='Microsoft YaHei')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy.fft as fft
import sys
import heapq



#%%%基本几何参数--------------------
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
x_end= 0.4E8 #天梯端点所处的位置，单位m，前面已经有了

#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
#c=1e-9#阻尼系数，画图用不着
E=9E11 #材料的弹性模量，画图用不着
#k=np.ones(N+1) #胡克系数初始化
xi=np.arange(50000000,130000000,1000)  #各节点的x轴位置
xii=np.array([5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,12.0])*10000000

yii=[1645,1422,1256,1122,992,897,802,706,610,508,403,278,113]
#print("xi 这个数组是这样的 ", xi)
def v(x):
    return np.sqrt((2*mu_e/(x-x_e)-2*mu_m/(x-x_m)-(x*omega)**2)/2-(2*mu_e/(x_ms-x_e)-2*mu_m/(x_ms-x_m)-(x_ms*omega)**2)/2)


#%%绘制速度变化曲线----------
fig, ax = plt.subplots(figsize=(8,3))

yi=[]
for i in xi:
    yi.append(v(i))


ax.plot(xi,yi,color='r', linewidth=2,label='Eq.(5-12)')
ax.scatter(xii,yii,label='Numerical Simulation')
ax.legend()
#ax.set_xticks(np.arange(0,t_result[-1]+10+864000,864000),labels=uu2)
#ax.text(15*86400, 150,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p), fontsize=10,ha="left", va="center",  color="black")
ax.text(15*86400, 50, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km', fontsize=10,ha="left", va="center",  color="black")
ax.set_ylabel('speed (m/s)')
ax.set_xlabel('$\mathdefault{x_{end}}$ (m)')

plt.show
plt.savefig('..\图像结果\长度与速度.svg',bbox_inches='tight')
