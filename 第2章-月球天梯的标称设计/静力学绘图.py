# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 22:45:11 2022

@author: Www
"""
#把Maple里画的不同长度等截面月球天梯的绘图移植到python中
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
Ln = np.array([20.0,30,50,90,110,150,200,250])*1000000
#打了小数点后面的0，是为了让数字成为浮点数，而不是整数

C_e = -4800000
C_m = 380000000
mu_e = 3.986E14
mu_m = 4.9028E12
omega = 2.6616996246359261859816393057264E-6
rhoxA=0.001 #意思是rho乘以A


fig,ax = plt.subplots(figsize=(11,7))
ax.set_ylim(-2000,15000)

x = np.linspace(0, 377000000, n)
xn=np.searchsorted(x,Ln) #端点在x数组中的位置
xaxe=np.zeros(n)

yy=np.linspace(-2000,15000, 10)
#%%---设置拉力函数
T = (-mu_e/(x-C_e)-x**2*omega**2/2+mu_m/(x-C_m))*rhoxA
Tl= (-mu_e/(Ln-C_e)-Ln**2*omega**2/2+mu_m/(Ln-C_m))*rhoxA

#下面这个y纯粹是画图用的中间变了
y=np.zeros(shape=(np.shape(T)[0],np.shape(Tl)[0]))
#%%---画线
for i in np.arange(np.shape(Ln)[0]): 
    y[:,i]=T-Tl[i]
ax.plot(x[xn[0]:], y[:,0][xn[0]:],color='grey',linestyle='solid')
ax.plot(x[xn[1]:], y[:,1][xn[1]:],color='grey')
ax.plot(x[xn[2]:], y[:,2][xn[2]:],color='grey')
ax.plot(x[xn[3]:], y[:,3][xn[3]:],color='grey')
ax.plot(x[xn[5]:], y[:,5][xn[5]:],color='grey')
ax.plot(x[xn[6]:], y[:,6][xn[6]:],color='grey')
ax.plot(x[xn[7]:], y[:,7][xn[7]:],color='grey')
ax.plot(x[xn[4]:], y[:,4][xn[4]:],linewidth=3,color='red')
#%%---横坐标标签
xtick=np.append(np.array([0,20.0,30,50,90,110,150,200,250])*1000000,321000000)
xticklabels=np.append(np.array([0,20.0,30,50,90,110,150,200,250]).astype(int).astype(str),'L1')

ax.set_xticks(xtick)
ax.set_xticklabels(xticklabels)
#%%---地球和月球
ax.fill_betweenx(yy,C_e-6500000, C_e+6500000, facecolor='blue')
ax.fill_betweenx(yy,C_m-3000000, C_m+3000000, facecolor='grey')

#%%---画网格线
ax.grid()
ax.plot(x,xaxe,linestyle='solid',color='black',linewidth=1.5)
ax.vlines(321000000,-2000,15000,linestyle='--',color='black',linewidth=1.5)# L1点

ax.text(C_e,15200,'Earth',horizontalalignment ='center')
ax.text(C_m,15200,'Moon',horizontalalignment ='center')
ax.text(321000000,15200,'Lagrange Point 1',horizontalalignment ='center')

ax.text(380000000,-2100,r"($10^6$ m)",horizontalalignment ='center',verticalalignment ='top')
ax.text(-33000000,15000,r"(N)",horizontalalignment ='right',verticalalignment ='top')

ax.text(0,-2600,"Earth-Moon CMC",horizontalalignment ='right',verticalalignment ='top',rotation=0)


ax.text(C_e+6500000,11000,'Earth Surface',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(C_m-3000000,7200,'Lunar Surface',horizontalalignment ='right',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(100000000,11500,'Constant Section Tethers with Various Length,\nEquivalently Various End Point \n （with Linear Density 1 g/m）',
        horizontalalignment ='left',verticalalignment ='top',
        bbox=dict(pad=0, ec="none", fc="none"))

ax.text(260000000,14200,'L20',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(260000000,9550,'L30',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(260000000,5400,'L50',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(260000000,2400,'L90',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(260000000,1600,'L110',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"),color='red')
ax.text(260000000,800,'L150',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(260000000,260,'L200',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
ax.text(300000000,0,'L250',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))

ax.set_xlabel('x')
ax.set_ylabel('Inner Tension Force')

plt.savefig('等截面变长度图.svg', format='svg')

