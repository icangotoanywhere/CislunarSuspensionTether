# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 22:57:37 2022

@author: Www
"""

#把Maple里画各类材料的绘图移植到python中
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


fig,ax = plt.subplots(3,1,figsize=(9,12))
ax[0].set_ylim(-250,200)
ax[1].set_ylim(-6,3)
ax[2].set_ylim(-30,15)





x = np.linspace(10000000, 377000000, n)
xn=np.searchsorted(x,Ln) #端点在x数组中的位置
xaxe=np.zeros(n)

yy=np.linspace(-2000,15000, 10)


#%%---截面比的对数
def A (x,s,rho):
    return np.exp((-mu_e/(x - C_e) - x**2*omega**2/2 + mu_m/(x - C_m))*rho/s)

rho=np.array([7.85,7.85,2.7,4.5,4.5,1.12,1.75,2.4,1.4,1,1,1,1.2,1.2,0.6,0.5])*1E3
s=np.array([0.3,0.1,0.35,1,0.5,0.1,2.5,3.43,3.6,0.5,60,6,0.4,0.2,0.17,0.12])*1E9
r=np.zeros(shape=(n,16))
for i in np.arange(16):
   r[:,i]=np.log10(A(x,s[i],rho[i])/A(50000000,s[i],rho[i]))






#%%---横坐标标签和网格
xtick=np.array([0,50,100,150,200,250,300,350])*1000000
xticklabels=np.array([0,50,100,150,200,250,300,350]).astype(int).astype(str)
ax[0].set_xticks(xtick)
ax[0].set_xticklabels(xticklabels)
ax[1].sharex(ax[0])
ax[2].sharex(ax[0])
ax[0].get_yticklabels()[0].set_color('red')
ax[0].get_yticklabels()[-1].set_color('red')
ax[1].get_yticklabels()[0].set_color('red')
ax[1].get_yticklabels()[-1].set_color('red')
ax[2].get_yticklabels()[0].set_color('red')
ax[2].get_yticklabels()[-1].set_color('red')

ax[0].grid()
ax[1].grid()
ax[2].grid()


ax[0].text(C_e,220,'Earth',horizontalalignment ='center')
ax[0].text(C_m,220,'Moon',horizontalalignment ='center')
ax[0].text(321000000,220,'Lagrange Point 1',horizontalalignment ='center')
#ax[0].text(380000000,-2100,r"($10^6$ m)",horizontalalignment ='center',verticalalignment ='top')
#ax[0].text(-33000000,15000,r"(N)",horizontalalignment ='right',verticalalignment ='top')
#ax[0].text(0,-2600,"Earth-Moon CMC",horizontalalignment ='right',verticalalignment ='top',rotation=0)
#ax[0].text(C_e+6500000,11000,'Earth Surface',horizontalalignment ='left',bbox=dict(pad=0, ec="none", fc="white"))
#ax[0].text(C_m-3000000,7200,'Lunar Surface',horizontalalignment ='right',bbox=dict(pad=0, ec="none", fc="white"))


ax[0].set_xlabel('x  ($10^6$ m)')
ax[0].set_ylabel(r"$\lg{\frac{A(x)}{A(50*10^6 m)}}$")
ax[1].set_xlabel('x  ($10^6$ m)')
ax[1].set_ylabel(r"$\lg{\frac{A(x)}{A(50*10^6 m)}}$")
ax[2].set_xlabel('x  ($10^6$ m)')
ax[2].set_ylabel(r"$\lg{\frac{A(x)}{A(50*10^6 m)}}$")

ax[0].text(15000000,120,'Metals',size=12,weight='black',ha ='left',va ='bottom',bbox=dict(pad=0, ec="none", fc="white"))
ax[1].text(15000000,2,'Composite Materials',size=12,weight='black',ha ='left',va ='bottom',bbox=dict(pad=0, ec="none", fc="white"))
ax[2].text(15000000,10,'Woody Materials',size=12,weight='black',ha ='left',va ='bottom',bbox=dict(pad=0, ec="none", fc="white"))







#%%---地球和月球
ax[0].fill_betweenx(yy,C_e-6500000, C_e+6500000, facecolor='blue')
ax[0].fill_betweenx(yy,C_m-3000000, C_m+3000000, facecolor='grey')
ax[0].plot(x,xaxe,linestyle='solid',color='black',linewidth=1.5)
ax[0].vlines(321000000,-2000,15000,linestyle='--',color='black',linewidth=1.5)# L1点

ax[1].fill_betweenx(yy,C_e-6500000, C_e+6500000, facecolor='blue')
ax[1].fill_betweenx(yy,C_m-3000000, C_m+3000000, facecolor='grey')
ax[1].plot(x,xaxe,linestyle='solid',color='black',linewidth=1.5)
ax[1].vlines(321000000,-2000,15000,linestyle='--',color='black',linewidth=1.5)# L1点

ax[2].fill_betweenx(yy,C_e-6500000, C_e+6500000, facecolor='blue')
ax[2].fill_betweenx(yy,C_m-3000000, C_m+3000000, facecolor='grey')
ax[2].plot(x,xaxe,linestyle='solid',color='black',linewidth=1.5)
ax[2].vlines(321000000,-2000,15000,linestyle='--',color='black',linewidth=1.5)# L1点


#%%---那些个线
ax[0].plot(x,r[:,0],label='1) High Strength Steel')
ax[0].plot(x,r[:,1],label='2) High Strength Steel')
ax[0].plot(x,r[:,2],label='3) Aluminum Alloy')
ax[0].plot(x,r[:,3],label='4) Titanium Alloy')
ax[0].plot(x,r[:,4],label='5) Titanium Alloy')
ax[1].plot(x,r[:,6],label='7) Carbon Fiber')
ax[1].plot(x,r[:,7],label='8) Fiberglass')
ax[1].plot(x,r[:,8],label='9) Kevlar')
ax[1].plot(x,r[:,10],label='11) Carbon Nanotube')
ax[1].plot(x,r[:,11],label='12) Carbon Nanotube')
ax[2].plot(x,r[:,12],label='13) Bamboo (wet)')
ax[2].plot(x,r[:,13],label='14) Bamboo (wet)')
ax[2].plot(x,r[:,14],label='15) Bamboo (dry)')
ax[2].plot(x,r[:,15],label='16) Pinewood')


ax[0].legend(loc=8)
ax[1].legend(loc=8)
ax[2].legend(loc=8)



plt.show()
plt.savefig('各类材料.svg', format='svg')