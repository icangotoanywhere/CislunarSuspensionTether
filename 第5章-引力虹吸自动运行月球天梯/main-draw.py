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
#%%%分段数--------------------

#name=sys.argv[1] #这里是批量运行传入名称参数必须的语法
name='例1' #这里是批量运行传入名称参数必须的语法
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




#%%读数据----
#xc_result=np.loadtxt('outxc-'+name+'.txt')
t_result=np.loadtxt('..\数值结果\算'+name+'-'+'outt.txt')
x_result=np.loadtxt('..\数值结果\算'+name+'-'+'outx.txt')
y_result=np.loadtxt('..\数值结果\算'+name+'-'+'outy.txt')
vx_result=np.loadtxt('..\数值结果\算'+name+'-'+'outvx.txt')
vy_result=np.loadtxt('..\数值结果\算'+name+'-'+'outvy.txt')
ce_result=np.loadtxt('..\数值结果\算'+name+'-'+'outce.txt')




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
#x_end= 0.4E8 #天梯端点所处的位置，单位m，前面已经有了
ls=2*(x_ms-x_end)/N #每段天梯的长度，单位m

#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
#c=1e-9#阻尼系数，画图用不着
E=9E11 #材料的弹性模量，画图用不着
#k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms,x_end-ls,-ls)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)




#%%%截面积-----------------
Area=0.0001 #平方米，1厘米粗的
#%%%每一段加上筐的质量----------
mass_b=ls*Area*rho
mass_p=c_mass_p*mass_b
#mass=ls*Area*rho
#%%%每段的胡克系数-------------
k=E*Area/ls

'''
#%%绘制截图-展开2旋转----------
fig,ax = plt.subplots(figsize=(8,3))
#ax.set_xlim(-x_m*0.1,x_m*5)
ax.set_ylim(-x_m*0.05,x_m*1.2)

imgm = plt.imread('moon3.png')#读月球
imge = plt.imread('earth.png')#画地球
L3=range(0,x_result.shape[1],6*12)  #
Ld=range(0,x_result.shape[1],30*24) #每5天

u = y_result
v = x_result
#uc = yc_result*np.cos(ang)-xc_result*np.sin(ang)+x_ms*ang
#vc = yc_result*np.sin(ang)+xc_result*np.cos(ang)

for i in L3:
   #ax.plot(u[:,i],v[:,i],'k',linewidth=0.1)
   ax.plot([u[-1,i]+i*600*omega*x_ms,u[0,i]+i*600*omega*x_ms],[v[-1,i],v[0,i]],color='r', linewidth=2)
   ax.scatter(u[:,i]+i*600*omega*x_ms,v[:,i],c='k',s=0.1)
   #ax.scatter(u[0,i],v[0,i],c='g',s=1)
   ax.scatter([i*600*omega*x_ms,ce_result[-3,i]+i*600*omega*x_ms],[x_ms,ce_result[-4,i]],c='black',s=6)
   ax.imshow(imgm,extent=(i*600*omega*x_ms-4e6,i*600*omega*x_ms+4e6,x_m-4e6,x_m+4e6),aspect='equal')#画月亮
   ax.imshow(imge,extent=(i*600*omega*x_ms-6e6,i*600*omega*x_ms+6e6,-6e6,6e6),aspect='equal')#画地球
   #ax.scatter(uc[i],vc[i],c='k')
#--设时间标--
xticklables=[]
uu=[] #标签位置
for i in Ld :
   xticklables.append('day'+str(int(i/144)))
   uu.append(i*600*omega*x_ms)
def ziji(x):
    return x
ax.set_xticks(uu,xticklables,minor=False,rotation=45)
#secax = ax.secondary_xaxis('top', functions=(ziji, ziji))
#secax.set_xlabel('u (m)')
#secax.set_xticks(np.arange(0, 3, 0.1)*1e9)
#ax.set_xlabel('time')
ax.set_ylabel('v (m)')
plt.tight_layout()
plt.show
plt.savefig('..\图像结果\算'+name+'-'+'展开.svg',bbox_inches='tight')

#end=time.time()
#print('Running time: %s Seconds'%(end-start))

'''
#%%绘制速度变化曲线----------
fig, ax = plt.subplots(figsize=(8,3))

v_result=np.sqrt(vx_result**2+vy_result**2)
v_mean=np.mean(v_result, axis=0)

ax.plot(t_result,v_mean,color='black', linewidth=2)

uu2=[]
for i in np.arange(0,t_result[-1]+10+864000,864000) :
   uu2.append('day'+str(int(i/86400)))

ax.set_xticks(np.arange(0,t_result[-1]+10+864000,864000),labels=uu2)
ax.text(15*86400, 150,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p), fontsize=10,ha="left", va="center",  color="black")
ax.text(15*86400, 50, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km', fontsize=10,ha="left", va="center",  color="black")
ax.set_ylabel('speed (m/s)')

plt.show
plt.savefig('..\图像结果\算'+name+'-'+'速度.svg',bbox_inches='tight')

#%%绘制应力变化曲线----------
fig, ax = plt.subplots(figsize=(8,3))

f_result=(np.sqrt((np.diff(np.append(x_result,x_result[0:1,:],axis=0),axis=0))**2+(np.diff(np.append(y_result,y_result[0:1,:],axis=0),axis=0))**2)-ls)*k
f_mean=np.mean(f_result, axis=0)

ax.plot(t_result,f_mean,color='black', linewidth=2)

uu2=[]
for i in np.arange(0,t_result[-1]+10+864000,864000) :
   uu2.append('day'+str(int(i/86400)))

ax.set_xticks(np.arange(0,t_result[-1]+10+864000,864000),labels=uu2)
ax.set_ylabel('mean inner force(N)')

ax.text(15*86400, 1.4e6,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p), fontsize=10,ha="left", va="center",  color="black")
ax.text(15*86400, 1.2e6, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km', fontsize=10,ha="left", va="center",  color="black")

plt.show
plt.savefig('..\图像结果\算'+name+'-'+'内力.svg',bbox_inches='tight')


#%%绘制最大应力变化曲线----------
fig, ax = plt.subplots(figsize=(8,3))

f_result=(np.sqrt((np.diff(np.append(x_result,x_result[0:1,:],axis=0),axis=0))**2+(np.diff(np.append(y_result,y_result[0:1,:],axis=0),axis=0))**2)-ls)*k
f_bigest=np.mean(np.sort(f_result, axis=0)[-10:-1,:],axis=0)

ax.plot(t_result,f_bigest,color='black', linewidth=2)

uu2=[]
for i in np.arange(0,t_result[-1]+10+864000,864000) :
   uu2.append('day'+str(int(i/86400)))

ax.set_xticks(np.arange(0,t_result[-1]+10+864000,864000),labels=uu2)
ax.set_ylabel('bigest inner force(N)')

ax.text(15*86400, 1.6e6,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p), fontsize=10,ha="left", va="center",  color="black")
ax.text(15*86400, 1.2e6, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km', fontsize=10,ha="left", va="center",  color="black")

plt.show
plt.savefig('..\图像结果\算'+name+'-'+'最大内力.svg',bbox_inches='tight')

#%%绘制最大的开度变化曲线----------
fig, ax = plt.subplots(figsize=(8,3))

y_result_sort=np.sort(y_result,axis=0)
zuida=np.mean(y_result_sort[0:10,:], axis=0)
zuixiao=np.mean(y_result_sort[-10:-1,:], axis=0)

kai_result=-zuida+zuixiao

ax.plot(t_result,kai_result,color='black', linewidth=2)

uu2=[]
for i in np.arange(0,t_result[-1]+10+864000,864000) :
   uu2.append('day'+str(int(i/86400)))

ax.set_xticks(np.arange(0,t_result[-1]+10+864000,864000),labels=uu2)
ax.set_ylabel('open distance(m)')

ax.text(15*86400, 2.5e6,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p), fontsize=10,ha="left", va="center",  color="black")
ax.text(15*86400, 1.2e6, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km', fontsize=10,ha="left", va="center",  color="black")

plt.show
plt.savefig('..\图像结果\算'+name+'-'+'开度.svg',bbox_inches='tight')


'''
#%%绘制某一时刻的图像----------
fig, ax = plt.subplots(figsize=(3,8))
plt.xlim(-x_m*0.2, x_m*0.2)  #由于我使用的坐标系的问题，要旋转90度来画图
plt.ylim(-x_m*0.2, x_m*1.2)

img = plt.imread('earth.png')#画地球

day=93

ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='equal')#画地球
ax.imshow(imgm,extent=(-4e6,4e6,x_m-4e6,x_m+4e6),aspect='equal')#画月亮
ax.plot(y_result[:,day*144], x_result[:,day*144],'k',linewidth=0.8)
for i in np.arange(0, N,10):
    if vx_result[i,day*144]<0:
        ax.scatter(y_result[i,day*144],x_result[i,day*144],c='r',marker='o',edgecolors='r',s=9,zorder=2)
    else:
        ax.scatter(y_result[i,day*144],x_result[i,day*144],c='black',marker='o',edgecolors='r',s=9,zorder=2)

ax.scatter(ce_result[-3,day*144],ce_result[-4,day*144],c='black',s=20)
ax.plot([y_result[-1,day*144],y_result[0,day*144]],[x_result[-1,day*144],x_result[0,day*144]],'black',linewidth=4)#首尾连接段

ax.text(0,-0.3E8,  '第 '+str(day)+' 天', fontsize=12,ha="center", va="center",  color="black")

plt.show
plt.savefig('..\图像结果\算'+name+'-'+'某一时刻.svg',bbox_inches='tight')

#%%绘制某一时刻的应力分布图像----------
fig, ax = plt.subplots(figsize=(3,8))
plt.xlim(-x_m*0.2, x_m*0.2)  #由于我使用的坐标系的问题，要旋转90度来画图
plt.ylim(-x_m*0.2, x_m*1.2)

img = plt.imread('earth.png')#画地球

day=100
color_result=np.maximum((f_result[:,day*144])/(5E6),0)*5
#color_result=color_result.astype(np.int)
ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='equal')#画地球
ax.imshow(imgm,extent=(-4e6,4e6,x_m-4e6,x_m+4e6),aspect='equal')#画月亮
#ax.plot(y_result[:,day*144], x_result[:,day*144],'k',linewidth=0.8)
ax.scatter(y_result[:,day*144],x_result[:,day*144],c=color_result,s=2,zorder=2,cmap=plt.cm.rainbow)

ax.scatter(ce_result[-3,day*144],ce_result[-4,day*144],c='black',s=20)
ax.plot([y_result[-1,day*144],y_result[0,day*144]],[x_result[-1,day*144],x_result[0,day*144]],'black',linewidth=4)#首尾连接段
ax.text(0,-0.3E8,  '第 '+str(day)+' 天', fontsize=12,ha="center", va="center",  color="black")
plt.show
plt.savefig('..\图像结果\算'+name+'-'+'某一时刻应力分布.svg',bbox_inches='tight')
'''

#%%四合一图-------------------
fig, ax = plt.subplots(2,2,figsize=(12,5))

x_ce_result=np.mean(ce_result[-4][-10000:-1])

uu2=[]
for i in np.arange(0,t_result[-1]+10+864000*2,864000*2) :
   uu2.append('day'+str(int(i/86400)))

ax1=plt.subplot(2,2,1)
ax1.plot(t_result,v_mean,color='black', linewidth=2)
ax1.set_xticks(np.arange(0,t_result[-1]+10+864000*2,864000*2),labels=uu2)
ax1.text(0.5, 0.1,  'Extreme speed: v ='+str(int(np.mean(v_mean[80*144:120*144])))+'m/s',transform = ax1.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax1.set_ylabel('speed (m/s)')


ax2=plt.subplot(2,2,3)
ax2.plot(t_result,kai_result,color='black', linewidth=2)
ax2.set_xticks(np.arange(0,t_result[-1]+10+864000*2,864000*2),labels=uu2)
ax2.text(0.5, 0.1,  'Open distance: d ='+str(int(np.mean(kai_result[80*144:120*144])/1000))+'km',transform = ax2.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax2.set_ylabel('open distance(m)')


ax3=plt.subplot(2,2,2)
ax3.plot(t_result,f_mean,color='black', linewidth=2)
ax3.set_xticks(np.arange(0,t_result[-1]+10+864000*2,864000*2),labels=uu2)
ax3.text(0.5, 0.1,  'Mean inner force: f ='+str(int(np.mean(f_mean[80*144:120*144])/1000))+'kN',transform = ax3.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax3.set_ylabel('mean inner force(N)')

ax4=plt.subplot(2,2,4)
ax4.plot(t_result,f_bigest,color='black', linewidth=2)
ax4.set_xticks(np.arange(0,t_result[-1]+10+864000*2,864000*2),labels=uu2)
ax4.text(0.05, 0.2,  'Bigest inner force: f ='+str(int(np.mean(f_bigest[80*144:120*144])/1000))+'kN',transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")

ax4.set_ylabel('bigest inner force(N)')
ax4.text(0.5, 0.3,  'Coefficient of payload mass: \u03BC ='+str(c_mass_p),transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax4.text(0.5, 0.2,  'Mass of one payload: m$_p$ ='+str(int(mass_p))+'kg',transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax4.text(0.5, 0.1, 'Length of tether: L= '+str((x_ms-x_end)/1000)+'km',transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax4.text(0.5, 0.5,  '$\mathdefault{x_{end}}$='+str(int(x_end/1000))+'km',transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")
ax4.text(0.5, 0.4,  'x$_{ce}$ ='+str(int(x_ce_result/1000))+'km',transform = ax4.transAxes, fontsize=10,ha="left", va="center",  color="black")


plt.tight_layout()
ax1.text(float(keys['alx']), float(keys['aly']), 'a', transform=ax1.transAxes,fontsize=30, ha="center", va="center",  color="darkgrey")
ax2.text(float(keys['blx']), float(keys['bly']), 'b', transform=ax2.transAxes,fontsize=30, ha="center", va="center",  color="darkgrey")
ax3.text(float(keys['clx']), float(keys['cly']), 'c', transform=ax3.transAxes,fontsize=30, ha="center", va="center",  color="darkgrey")
ax4.text(float(keys['dlx']), float(keys['dly']), 'd', transform=ax4.transAxes,fontsize=30, ha="center", va="center",  color="darkgrey")
plt.show
plt.savefig('../../论文图片\四合一'+str(x_end/1E7)+'-'+str(c_mass_p)+'-.svg',bbox_inches='tight')


#%%以下就要画动画啦，动画就是先画个图，然后在图的基础上更新-----
fig, ax = plt.subplots(figsize=(10,10))
plt.xlim(-x_m*0.2, x_m*0.2)  #由于我使用的坐标系的问题，要旋转90度来画图
plt.ylim(-x_m*0.2, x_m*1.2)
img = plt.imread('earth.png')#画地球
imgm = plt.imread('moon3.png')#画地球

ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='equal')#画地球
ax.imshow(imgm,extent=(-4e6,4e6,x_m-4e6,x_m+4e6),aspect='equal')#画月亮
basketcolor = np.where(vx_result <0, 'black', vx_result)
basketcolor[basketcolor != 'black'] =  'r'


def animate(i):
    line, = ax.plot(y_result[:,0], x_result[:,0],'k',linewidth=0.8) #画出初始图形
    line.set_ydata(x_result[:,i*2])  # update the data
    line.set_xdata(y_result[:,i*2])
    link, =ax.plot([y_result[-1,i*2],y_result[0,i*2]],[x_result[-1,i*2],x_result[0,i*2]],'r',linewidth=3)#首尾连接段
    basket =ax.scatter(y_result[0:N:20,i*2],x_result[0:N:20,i*2],c=basketcolor[0:N:20,i*2],marker='o',edgecolors='r',s=9,zorder=2)
    riqi=ax.text(0,-0.3E8,  '第 '+str('{:.3f}'.format(i/144))+' 天', fontsize=12,ha="center", va="center",  color="black")
    #climber=ax.scatter(yc_result[i*2],xc_result[i*2],c='r')
    ce=ax.scatter(ce_result[-3,i*2],ce_result[-4,i*2],c='r',s=3)
    return line,ce,link,basket,riqi
# Init only required for blitting to give a clean slate.
def init():
    line, = ax.plot(y_result[:,0], x_result[:,0],'k',linewidth=0.8) #画出初始图形

    #basket, =ax.scatter(y_result[:,0],x_result[:,0],c=basketcolor[:,0],marker='o',edgecolors='r',s=9,zorder=2)

    #line.set_ydata(x_result[:,0])
    #line.set_xdata(y_result[:,0])
    #line2, = ax.plot(y_result_kong[:,0], x_result_kong[:,0],'k',linewidth=0.8,color='red') #画出初始图形
    #line2.set_ydata(x_result_kong[:,0])
    #line2.set_xdata(y_result_kong[:,0])
    ce=ax.scatter(ce_result[-3,0],ce_result[-4,0],c='r',s=3)
    link, =ax.plot([y_result[-1,0],y_result[0,0]],[x_result[-1,0],x_result[0,0]],'r',linewidth=3)#首尾连接段
    riqi=ax.text(0,-0.3E8,  '第 '+str(0)+' 天', fontsize=12,ha="center", va="center",  color="black")
    #climber=ax.scatter(yc_result[0],xc_result[0],c='r')
    return line,ce,link,riqi
# call the animator.  blit=True means only re-draw the parts that have changed.
# blit=True dose not work on Mac, set blit=False
# interval= update frequency

ani = animation.FuncAnimation(fig=fig,
                              func=animate,
                              interval=1,
                              frames=int(len(t_result)),
                              init_func=init,
                              blit=True)
#ani.save('..\图像结果\算'+name+'-'+'动画.mp4',writer='ffmpeg')
#ani.save("tmp.gif", writer='pillow')

plt.show()

