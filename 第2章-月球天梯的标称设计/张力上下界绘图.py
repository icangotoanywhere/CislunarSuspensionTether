# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 13:10:49 2022

@author: Www
"""
#把Maple里画的月球天梯张力上下限的绘图移植到python中
import numpy as np
import math
#import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib
matplotlib.rc("font", family='Microsoft YaHei')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import patheffects

n=10000

C_e = -4800000
C_m = 380000000
mu_e = 3.986E14
mu_m = 4.9028E12
omega = 2.6616996246359261859816393057264E-6
rhoxA=0.001 #意思是rho乘以A

x = np.linspace(15000000, 321000000, n)
xaxe=np.zeros(n)
yy=np.linspace(-2500,18000, 10)


#%%---设置拉力函数
def T(x):
    return (-mu_e/(x - C_e) - 1/2*x**2*omega**2 + mu_m/(x - C_m))*rhoxA

Tmax=T(321919886.3)-T(x)
Tmin=T(377500000)-T(x)



fig,ax = plt.subplots(figsize=(11,7))
ax.set_ylim(-2500,18000)
ax.set_xlim(-15000000,390000000)

xtick=np.array([0,50,100,150,200,250,300,350])*1000000
xticklabels=np.array([0,50,100,150,200,250,300,350]).astype(int).astype(str)
ax.set_xticks(xtick)
ax.set_xticklabels(xticklabels)


#%%---地球和月球
ax.fill_betweenx(yy,C_e-6500000, C_e+6500000, facecolor='blue')
ax.fill_betweenx(yy,C_m-3000000, C_m+3000000, facecolor='grey')



ax.plot(x, Tmax, color='black',linewidth=2.5,label='Tmax, at L1')
ax.plot(x, Tmin, color='black',linewidth=2.5,linestyle=(0, (2, 1, 1, 1)),label='Tmin at $x_{ms} $')

ax.fill_between(x, Tmax, Tmin, where=Tmax>=Tmin, facecolor='green')
ax.vlines(321000000,-2500,18000,linestyle='--',color='black',linewidth=1.5)# L1点
ax.hlines(0,-15000000,390000000,linestyle='solid',color='black',linewidth=1.5)

ax.vlines(110842000,-2500,16000,linestyle='--',color='black',linewidth=1)# L1点
ax.vlines(40000000,-2500,12000,linestyle='--',color='black',linewidth=1)# L1点
ax.vlines(90000000,-2500,12000,linestyle='--',color='black',linewidth=1)# L1点
ax.vlines(30000000,-2500,16000,linestyle='--',color='black',linewidth=1)# L1点



ax.text(C_e,18200,'Earth',horizontalalignment ='center')
ax.text(C_m,18200,'Moon',horizontalalignment ='center')
ax.text(321000000,18200,'Lagrange Point 1',horizontalalignment ='center')

#ax.text(380000000,-2600,r"($10^6$ m)",horizontalalignment ='center',verticalalignment ='top')
#ax.text(-33000000,15000,r"(N)",horizontalalignment ='right',verticalalignment ='top')

ax.text(0,-3800,"Earth-Moon CMC",horizontalalignment ='center',verticalalignment ='top',rotation=0)


ax.text(C_e+6800000,1000,'Earth Surface',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"),rotation=90.0)
ax.text(C_m-3300000,1000,'Lunar Surface',horizontalalignment ='right',bbox=dict(pad=0, ec="none", fc="white"),rotation=90.0)




ax.set(xlabel='End-Point Position of the Tethers, $ x_{end} $  ($10^6$ m)', ylabel='black')
ax.set_ylabel('Inner Tension Force  (N)')


ax.hlines(11000,40000000,90000000,linestyle='solid',color='black',linewidth=1.5)
ax.arrow(80000000,11000, 10000000, 0,width=1,length_includes_head=True,head_width=400,head_length=5000000,fc='black',ec='none')
ax.arrow(50000000,11000, -10000000, 0,width=1,length_includes_head=True,head_width=400,head_length=5000000,fc='black',ec='none')
ax.text(65000000,11200,'Proper $ x_{end} $',horizontalalignment ='center',verticalalignment ='bottom',bbox=dict(pad=0, ec="none", fc="white"))


ax.hlines(14000,30000000,110000000,linestyle='solid',color='black',linewidth=1.5)
ax.arrow(100000000,14000, 10000000, 0,width=1,length_includes_head=True,head_width=400,head_length=5000000,fc='black',ec='none')
ax.arrow(40000000,14000, -10000000, 0,width=1,length_includes_head=True,head_width=400,head_length=5000000,fc='black',ec='none')
ax.text(75000000,14200,'Valid $ x_{end} $',horizontalalignment ='center',verticalalignment ='bottom',bbox=dict(pad=0, ec="none", fc="white"))



ax.grid()
ax.legend(loc=10)


plt.show()

plt.savefig('最大最小张力.svg', format='svg')