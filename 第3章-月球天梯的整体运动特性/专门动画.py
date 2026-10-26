# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 00:17:20 2022

@author: Www
"""

# %%这是一个关于月球天梯动力学的程序
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

#%%%分段数，文件名--------------------
#N=340#质点总数量,段数，月表为0
name='3400-圆-粘性-10-4万公里-无偏角-随月启动'

#%%读数据----
x_result=np.loadtxt('outx-'+name+'.txt')
y_result=np.loadtxt('outy-'+name+'.txt')
#yc_result=np.loadtxt('outyc-'+name+'.txt')
#xc_result=np.loadtxt('outxc-'+name+'.txt')
t_result=np.loadtxt('outt-'+name+'.txt')

N=np.shape(x_result)[0]-1#质点总数量,段数，月表为0
num =np.arange(0,N+1,1) #质点序号数组
print (N)
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
x_end= 0.4E8 #天梯端点所处的位置，单位m
l_step=(x_ms-x_end)/N #每段天梯的长度，单位m

#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
c=1e-9#阻尼系数
E=9E11 #材料的弹性模量
k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms,x_end-l_step,-l_step)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)

#%%%以下是研究变截面的情况-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式
    return 0.01*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
#for x in xi:
#    print(x,"处的截面积为",Area(x))

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

#%%绘制截图-圆形----------
fig,ax = plt.subplots(figsize=(7,7))
ax.set_xlim(-x_m*1.2,x_m*1.2)
ax.set_ylim(-x_m*1.2,x_m*1.2)

img = plt.imread('earth.png')#画地球
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球

L=range(0,x_result.shape[1] ,6*3)  #6代表1小时
for i in L:
   ax.plot(y_result[:,i],x_result[:,i],'k',linewidth=0.8)
   #ax.scatter(yc_result[i],xc_result[i],c='k')

imgm = plt.imread('moon3.png')#读月球
for i in L:
   ax.imshow(imgm,extent=(y_result[0,i]-4e6,y_result[0,i]+4e6,x_result[0,i]-4e6,x_result[0,i]+4e6),aspect='auto')#画月亮

#----
ha=["left","right"]
va=["bottom","top"]
for i in range(0,x_result.shape[1] ,6*24):  #6代表1小时
    ax.text(y_result[0,i],x_result[0,i],str(int(i/144))+'day',\
            ha=ha[bool(i > 14*144)],va=va[7*144 < i < 21*144],\
            size=10) #布尔值用来使文字出现在外围
#----
#ax.set_xlabel('m')
#ax.set_title('Earth-Moon Planer',y=-0.12)
plt.show
plt.savefig('例1圆饼.svg',bbox_inches='tight')
#%%绘制截图-展开1映射----------
'''
fig,ax = plt.subplots(figsize=(20.4,5.2))
ax.set_xlim(-x_m*0.1,x_m*5)
ax.set_ylim(-x_m*0.05,x_m*1.1)

img = plt.imread('earth.png')#画地球
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球


#---坐标转换---
v=np.sqrt(y_result**2+x_result**2)
ang=np.angle(x_result+y_result*1j)
ang2=np.where(ang < -0.1,ang+2*math.pi,ang)
u=x_ms*ang2  #展开后的y


imgm = plt.imread('moon3.png')#读月球
#--画线画月亮--
L3=range(0,x_result.shape[1]-200 ,6*3)  #6代表1小时，这是3小时
Ld=range(0,x_result.shape[1]-100,6*24) #1天
for i in L3:
   ax.plot(u[:,i],v[:,i],'k',linewidth=0.8)
   ax.imshow(imgm,extent=(u[0,i]-4e6,u[0,i]+4e6,v[0,i]-4e6,v[0,i]+4e6),aspect='auto')#画月亮
   ax.scatter(yc_result[i],xc_result[i],c='k')
#--设时间标--
xticklables=[]
uu=[]
for i in Ld :
   xticklables.append('t='+str(int(i/144))+'day')
   uu.append(u[0,i])
ax.set_xticks(uu,xticklables,minor=False,rotation=45)

#----
#ax.set_xlabel('time')
ax.set_title('Extented to time-paner',y=1)
plt.show
'''
#%%绘制截图-展开2旋转----------
fig,ax = plt.subplots(figsize=(14,5))
ax.set_xlim(-x_m*0.1,x_m*5)
ax.set_ylim(-x_m*0.05,x_m*1.1)
imge = plt.imread('earth.png')#画地球
#---坐标转换---
ang0=np.angle(x_result[0,:]+y_result[0,:]*1j)
ang=np.where(ang0 < 0,ang0+2*math.pi,ang0)
u = y_result*np.cos(ang)-x_result*np.sin(ang)+x_ms*ang
v = y_result*np.sin(ang)+x_result*np.cos(ang)
#uc = yc_result*np.cos(ang)-xc_result*np.sin(ang)+x_ms*ang
#vc = yc_result*np.sin(ang)+xc_result*np.cos(ang)
imgm = plt.imread('moon3.png')#读月球
#--画线画月亮--
L3=range(0,x_result.shape[1]-200 ,6*3)  #6代表1小时，这是3小时
Ld=range(0,x_result.shape[1]-100,6*24) #1天
for i in L3:
   ax.plot(u[:,i],v[:,i],'k',linewidth=0.8)
   ax.imshow(imgm,extent=(u[0,i]-4e6,u[0,i]+4e6,v[0,i]-4e6,v[0,i]+4e6),aspect='equal')#画月亮
   ax.imshow(imge,extent=(u[0,i]-6e6,u[0,i]+6e6,-6e6,6e6),aspect='equal')#画地球
   #ax.scatter(uc[i],vc[i],c='k')
#--设时间标--
xticklables=[]
uu=[]
for i in Ld :
   xticklables.append(str(int(i/144))+'day')
   uu.append(u[0,i])
def ziji(x):
    return x
ax.set_xticks(uu,xticklables,minor=False,rotation=45)
secax = ax.secondary_xaxis('top', functions=(ziji, ziji))
secax.set_xlabel('u (m)')
secax.set_xticks(np.arange(0, 3, 0.1)*1e9)
#ax.set_xlabel('time')
ax.set_ylabel('v (m)')
plt.show
plt.savefig('例1展开.svg',bbox_inches='tight')
#%%以下就要画动画啦，不过这种方法涉及到imshow，非常慢-----
'''
fig, ax = plt.subplots(figsize=(10,10))
plt.xlim(-x_m*1.2, x_m*1.2)  #由于我使用的坐标系的问题，要旋转90度来画图
plt.ylim(-x_m*1.2, x_m*1.2)
img = plt.imread('earth.png')#画地球
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球
imgm = plt.imread('moon3.png')#读月球

artists=[]
for i in range(0,len(sol.t)):
    line, = ax.plot(y_result[:,0], x_result[:,0],'k',linewidth=0.8) #画出初始图形
    line.set_ydata(x_result[:,i])  # update the data
    line.set_xdata(y_result[:,i])
    lunar=ax.scatter(y_result[0,i],x_result[0,i],c='k')
    #ims=ax.imshow(imgm,extent=(y_result[0,i]-5e6,y_result[0,i]+5e6,x_result[0,i]-5e6,x_result[0,i]+5e6),aspect='auto')
    #artists.append([line,ims])
    artists.append([lunar,line])

ani = animation.ArtistAnimation(fig=fig, artists=artists, repeat=False, interval=1)

plt.show()
'''
#%%静态的数据图，展开转换后各点的u值和v值-------------------
fig, ax = plt.subplots(2,1,sharex='all',figsize=(12,7))
j=5 #选取的数量
biaoqian=["the root","the upper-middle point","the middle point","the lower-middle point","the end"]
ii=0
s=np.arange(j) #选取点的序号
for i in range(j):
    s[i]=np.array_split(range(N+1),j)[i][0]#平均分成j段，取各段第一个点
ax1=plt.subplot(2,1,1)
for i in s:
    plt.plot(t_result, u[i]-x_ms*ang, label=biaoqian[ii])
    #plt.xlabel('$time (s)$') # the horizontal axis represents the time
    plt.ylabel('$u-u_0$ (m)')
    ii=ii+1
plt.legend(fontsize=10,ncol=j,bbox_to_anchor=(0, 0), loc=2, borderaxespad=1)
def day(x):
   return x/24/3600
secax = ax1.secondary_xaxis('top', functions=(day, day))
secax.set_xlabel('time (day)')
secax.set_xticks(np.arange(0, 28, 1))
ax2=plt.subplot(2,1,2)
ii=0
for i in s:
    plt.plot(t_result, v[i], label=biaoqian[ii])
    plt.ylabel('v (m)')
    plt.xlabel('$time (s)$') # the horizontal axis represents the time
    ii=ii+1
plt.show()
plt.savefig('例1uv.svg',bbox_inches='tight')
#%%五个点五条线的标签

#%%静态的数据图，典型点的应力-------------------
fig, ax = plt.subplots(1,1,figsize=(10,5))
s=np.arange(j) #选取点的序号
biaoqian=["the root","the upper-middle point","the middle point","the lower-middle point","the end"]
ii=0
l_step=np.sqrt((x_result[10][0]-x_result[11][0])**2+(y_result[10][0]-y_result[11][0])**2)
for i in range(j):
    s[i]=np.array_split(range(N+1),j)[i][0]#平均分成j段，取各段第一个点
for i in s+5:
    stress=(np.sqrt((x_result[i]-x_result[i+1])**2+(y_result[i]-y_result[i+1])**2)-l_step)/(Area((xi[i]+xi[i])/2))
    plt.plot(t_result, stress, label=biaoqian[ii])
    ii=ii+1
plt.axhline(y=0.0, color="black", linestyle="--")
plt.xlabel('$time (s)$') # the horizontal axis represents the time
plt.ylabel('stress')
def day(x):
    return x/24/3600
secax = ax.secondary_xaxis('top', functions=(day, day))
secax.set_xlabel('time (day)')

plt.legend()
plt.show()
plt.savefig('例1应力.svg',bbox_inches='tight')
#%%静态的数据图，典型点的傅里叶变换-------------------
fig, ax = plt.subplots(1,1,figsize=(9,5))
ax=plt.subplot(1,1,1)
s=np.arange(j) #选取点的序号
for i in range(j):
    s[i]=np.array_split(range(N+1),j)[i][0]#平均分成j段，取各段第一个点
ii=0
for i in s:
    us=u[i]-u[0]#相对u位移
    #stress=(np.sqrt((x_result[i]-x_result[i+1])**2+(y_result[i]-y_result[i+1])**2)-l_step)/(Area((xi[i]+xi[i])/2))
    ftn=np.int64(np.power(2,np.ceil(np.log2(t_result.size))))#取采样数的下一个2的n次方
    freq = np.fft.fftfreq(ftn, d=600)#频率
    absfft=np.abs(2/ftn*fft.fft(us,ftn))#fft的幅值
    ax.plot(freq, absfft,label=biaoqian[ii])
    ax.set_xlim(0.0000005,0.000015)
    ii=ii+1
ax.set_xlabel('$frequency (Hz)$') # the horizontal axis represents the time
ax.set_ylabel('Amplitude (m)')
plt.axvline(x=2.848e-6, color="black", linestyle="--")
plt.axvline(x=1.06e-5, color="black", linestyle="--")
def invert(x):
    # 1/x with special treatment of x == 0
    x = np.array(x).astype(float)
    near_zero = np.isclose(x, 0)
    x[near_zero] = np.inf
    x[~near_zero] = 1 / x[~near_zero]
    return x/24/3600

secax = ax.secondary_xaxis('top', functions=(invert, invert))
secax.set_xlabel('Period (day)')
secax.set_xticks([0.5, 1, 1.5, 2,3,4,5,10])
plt.legend()
plt.show()
plt.savefig('例1傅里叶.svg',bbox_inches='tight')
#%%以下就要画动画啦，动画就是先画个图，然后在图的基础上更新-----
fig, ax = plt.subplots(figsize=(10,10))
plt.xlim(-x_m*1.2, x_m*1.2)  #由于我使用的坐标系的问题，要旋转90度来画图
plt.ylim(-x_m*1.2, x_m*1.2)
img = plt.imread('earth.png')#画地球
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球

def animate(i):
    line, = ax.plot(y_result[:,0], x_result[:,0],'k',linewidth=0.8) #画出初始图形
    line.set_ydata(x_result[:,i*3])  # update the data
    line.set_xdata(y_result[:,i*3])
    lunar=ax.scatter(y_result[0,i*3],x_result[0,i*3],c='k')
    #climber=ax.scatter(yc_result[i*3],xc_result[i*3],c='r')
    return lunar,line#,climber
# Init only required for blitting to give a clean slate.
def init():
    line, = ax.plot(y_result[:,0], x_result[:,0],'k',linewidth=0.8) #画出初始图形
    line.set_ydata(x_result[:,0])
    line.set_xdata(y_result[:,0])
    lunar=ax.scatter(y_result[0,0],x_result[0,0],c='k')
   # climber=ax.scatter(yc_result[0],xc_result[0],c='r')
    return lunar,line#,climber
# call the animator.  blit=True means only re-draw the parts that have changed.
# blit=True dose not work on Mac, set blit=False
# interval= update frequency

ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              interval=1,
                              frames=int(len(t_result)/3),
                              init_func=init,
                              blit=True)
#ani.save('动画'+str(N)+'-'+name+'.mp4',writer='ffmpeg')
plt.show()


