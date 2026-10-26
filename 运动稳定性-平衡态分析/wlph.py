import numpy as np
import math
import scipy.integrate as si #引入积分用的函数
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import solve_ivp
from numpy import array,mat,sin,cos,exp
from scipy.optimize import root
N=20

num =np.arange(0,N+1,1) #质点序号数组
#for i  in num: print(num[i])
print ("num=",num)
#-----基本几何参数--------------------
mu_e = 3.986E14 #地球引力常数，单位m^3*s^(-2)
mu_m = 4.9028E12 #月球引力常数，单位m^3*s^(-2)
mu_m =0
omega = 2*math.pi/(27.32*24*3600) #月球公转角速度，单位S^(-1)
omega = 0
epsilon=0.0  #月球公转偏心率
theta_0=0.00
#--这几个几何参数仅用于计算绳子的参数--
x_e=0 #地心，单位m
y_e=0
x_m = 3.84403E8 #月心距离地心的距离，单位m
y_m=0
x_ms = 3.7826E8+4403000 #月面距离地心的距离，单位m
# paras = [mu_e,mu_m,omega,x_e,x_m,x_ms]
# print(paras)
x_end= 0.5E8 #天梯端点所处的位置，单位m,放到前面参数表去了
l_s=(x_ms-x_end)/N #每段天梯的长度，单位m


#%%%材料参数----------
rho = 1E3 #材料密度，单位kg*m^(-3)
S = 6000000000 #材料抗拉强度，单位N*m^(-2)，即Pa
c=0.01#阻尼系数
E=9E11 #材料的弹性模量
k=np.ones(N+1) #胡克系数初始化
xi=np.arange(x_ms,x_end-l_s,-l_s)  #各节点的x轴位置
#print("xi 这个数组是这样的 ", xi)

#%%%以下是研究变截面的情况-----------------
def Area(x): #定义x处的横截面积Area
    # 变截面（恒定应力）情况下的截面积表达式
    return 0.01*math.exp((-mu_e/(x-x_e)+mu_m/(x-x_m)-(x**2)*(omega**2)/2)*rho/S)
#for x in xi:
#    print(x,"处的截面积为",Area(x))

#%%%每个点的质量----------
def mass2(x): #定义某一段的质量
    return sp.integrate.quad(Area,x,x+l_s)[0]*rho #使用scipy的积分
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

x=(xi-xi[0])*math.cos(theta_0)+x_ms*(1-epsilon)
y=(xi[0]-xi)*math.sin(theta_0)
x0=x
y0=y
xy0=np.append(x0[1:],y0[1:])
#这是把第0个点贴回来，就是月面那个点
#%%五力平衡的定义-------------
def F5(xy):  #定义五力平衡非线性代数方程组
    #从第1个点到第N个点，不考虑0点月面
    #这是把第0个点贴回来，就是月面那个点
    FX=np.zeros(N+1) #初始化x方向合力
    FY=np.zeros(N+1) #初始化y方向合力
    x=np.append(xi[0],xy[:N]) #初始化
    y=np.append(0,xy[N:]) #初始化
    de2=(x-x_e)**2+(y-y_e)**2 #初始化到地球距离平方
    de=np.sqrt(de2) #初始化到地球距离
    dm2=(x-x_m)**2+(y-y_m)**2 #初始化到月球距离平方
    dm=np.sqrt(dm2) #初始化到月球距离
    ll=np.append(0,np.sqrt(np.diff(x)**2+np.diff(y)**2)) #初始化各段新长度
    #lr=np.append(0,np.maximum(0,1-l_s/ll[1:])) #初始化各段伸长率
    lr=np.append(0,1-l_s/ll[1:]) #初始化各段伸长率


    for i in np.arange(1,N): #从1到N-1
        FX[i]=mu_e*mass[i]/de2[i]*(x_e-x[i])/de[i]\
              +mu_m*mass[i]/dm2[i]*(x_m-x[i])/dm[i]\
              +mass[i]*x[i]*omega**2\
              +k[i]*lr[i]*(x[i-1]-x[i])\
              +k[i+1]*lr[i+1]*(x[i+1]-x[i])
        FY[i]=mu_e*mass[i]/de2[i]*(y_e-y[i])/de[i]\
              +mu_m*mass[i]/dm2[i]*(y_m-y[i])/dm[i]\
              +mass[i]*y[i]*omega**2\
              +k[i]*lr[i]*(y[i-1]-y[i])\
              +k[i+1]*lr[i+1]*(y[i+1]-y[i])
    #最后一个点的写法不同，如下
    FX[N]=mu_e*mass[N]/de2[N]*(x_e-x[N])/de[N]\
          +mu_m*mass[N]/dm2[N]*(x_m-x[N])/dm[N]\
          +mass[N]*x[N]*omega**2\
          +k[N]*lr[N]*(x[N-1]-x[N])
    FY[N]=mu_e*mass[N]/de2[N]*(y_e-y[N])/de[N]\
          +mu_m*mass[N]/dm2[N]*(y_m-y[N])/dm[N]\
          +mass[N]*y[N]*omega**2\
          +k[N]*lr[N]*(y[N-1]-y[N])
    #第一个点直接就是0算求

    FXY=np.append(FX[1:],FY[1:]) #贴成一个向量，从第1个点开始，不管0点
    return FXY  #输出的是，一个标量

result=root(F5,xy0,method='lm')
print("初值=",xy0)
print("结果x=",result.x[:N])
print("结果y=",result.x[N:])


#%%绘制截图-圆形----------
fig,ax = plt.subplots(figsize=(10,10))
ax.set_xlim(-x_m*1.2,x_m*1.2)
ax.set_ylim(-x_m*1.2,x_m*1.2)

#img = plt.imread('earth.png')#画地球
#ax.imshow(img,extent=(-1e7,1e7,-1e7,1e7),aspect='auto')#画地球

ax.plot(result.x[:N],result.x[N:],'k',linewidth=0.8)
ax.plot(x,y,'g',linewidth=2)
#%%雅克比矩阵------------
yakebi=np.zeros((4*N+4,4*N+4))


x=np.append(x_ms,result.x[:N])
y=np.append(0,result.x[N:])
for i in np.arange(1,N,1):
    yakebi[i*4+2,i*4]=2*mu_e*mass[i]/((np.absolute(x_e-x[i]))**3)\
                   +2*mu_m*mass[i]/((np.absolute(x_m-x[i]))**3)\
                   +mass[i]*omega**2+k[i]+k[i+1]+c*omega
    yakebi[i*4+2,(i-1)*4]=k[i]
    yakebi[i*4+2,(i+1)*4]=k[i+1]
    yakebi[i*4+2,i*4+2]=-c
    yakebi[i*4+3,i*4+1]=-mu_e*mass[i]/((np.absolute(x_e-x[i]))**3)\
                       -mu_m*mass[i]/((np.absolute(x_m-x[i]))**3)\
                       +k[i]+k[i]*l_s/np.absolute(x[i-1]-x[i])\
                       +k[i+1]+k[i+1]*l_s/np.absolute(x[i+1]-x[i])\
                       +mass[i]*omega**2+c*omega
    yakebi[i*4+3,(i-1)*4+1]=k[i]+k[i]*l_s/np.absolute(x[i-1]-x[i])
    yakebi[i*4+3,(i+1)*4+1]=k[i+1]+k[i+1]*l_s/np.absolute(x[i+1]-x[i])
    yakebi[i*4+3,i*4+3]=-c
    yakebi[i*4,i*4+2]=1
    yakebi[i*4+1,i*4+3]=1
#最后一个N点不同
yakebi[N*4+2,N*4]=2*mu_e*mass[N]/((np.absolute(x_e-x[N]))**3)\
               +2*mu_m*mass[N]/((np.absolute(x_m-x[N]))**3)\
               +mass[N]*omega**2+k[N]+c*omega
yakebi[N*4+2,(N-1)*4]=k[N]
yakebi[N*4+2,N*4+2]=-c
yakebi[N*4+3,N*4+1]=-mu_e*mass[N]/((np.absolute(x_e-x[N]))**3)\
                   -mu_m*mass[N]/((np.absolute(x_m-x[N]))**3)\
                   +k[N]+k[N]*l_s/np.absolute(x[N-1]-x[N])\
                   +mass[N]*omega**2+c*omega
yakebi[N*4+3,(N-1)*4+1]=k[N]+k[N]*l_s/np.absolute(x[N-1]-x[N])
yakebi[N*4+3,N*4+3]=-c
yakebi[N*4,N*4+2]=1
yakebi[N*4+1,N*4+3]=1

yakebi=yakebi[4:,4:]
print(yakebi)


tezhengzhi=np.linalg.eig(yakebi)[0].reshape(N,4)
print('雅克比矩阵的特征值为',tezhengzhi)
