# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 22:09:11 2022

@author: Www
"""

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

C_e = -4800000
C_m = 380000000
mu_e = 3.986E14
mu_m = 4.9028E12
omega = 2.6616996246359261859816393057264E-6
rhoxA=0.001 #意思是rho乘以A


fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111, projection='3d')

xtick=np.array([0,50,100,150,200,250,300,350])*1000000
xticklabels=np.array([0,50,100,150,200,250,300,350]).astype(int).astype(str)
ax.set_xticks(xtick)
ax.set_xticklabels(xticklabels)
ax.set_xlabel('x  ($10^6$ m)')
ax.set_ylabel("radius of section (m)")
ax.set_yticks(np.linspace(-2, 2, 10))
ax.set_yticklabels([])
ax.set_zticks(np.linspace(-2, 2, 10))
ax.set_zticklabels([])

n=100
t = np.linspace(10000000, 377000000, n)

#第8种材料，fiberglass
rho=2.4E3
s=3.43E9
#%%---设定3万公里处的圆面半径是0.1m
def radius(x,rho,s):
    return np.sqrt(np.exp((-mu_e/(x - C_e) - x**2*omega**2/2 + mu_m/(x - C_m))*rho/s)/np.pi)

t8=radius(t,rho=2.4E3,s=3.43E9)/radius(1E8,rho=2.4E3,s=3.43E9)*0.1
t9=radius(t,rho=1.4E3,s=3.6E9)/radius(1E8,rho=1.4E3,s=3.6E9)*0.1
t11=radius(t,rho=1E3,s=60E9)/radius(1E8,rho=1E3,s=60E9)*0.1
t12=radius(t,rho=1E3,s=6E9)/radius(1E8,rho=1E3,s=6E9)*0.1




u = np.linspace(0,2*np.pi,100)  # 把圆分按角度为50等分
h = np.linspace(0,1,100)        # 把高度1均分为20份
z8 = np.outer(np.sin(u),t8)  # x值重复20次
y8 = np.outer(np.cos(u),t8)  # y值重复20次
x8 = np.outer(np.ones(len(u)),t)   # x，y 对应的高度
z9 = np.outer(np.sin(u),t9)  # x值重复20次
y9 = np.outer(np.cos(u),t9)  # y值重复20次
x9 = np.outer(np.ones(len(u)),t)   # x，y 对应的高度
z11 = np.outer(np.sin(u),t11)  # x值重复20次
y11 = np.outer(np.cos(u),t11)  # y值重复20次
x11 = np.outer(np.ones(len(u)),t)   # x，y 对应的高度
z12 = np.outer(np.sin(u),t12)  # x值重复20次
y12 = np.outer(np.cos(u),t12)  # y值重复20次
x12 = np.outer(np.ones(len(u)),t)   # x，y 对应的高度

ax.set_xlim(0,380000000)
#ax.set_xlim(-2.5,2.5)
ax.set_ylim(-2,2)
ax.set_zlim(-2,2)
#ax.set_aspect('equal')

ax.plot_surface(x8, y8-1, z8,alpha=0.99, cmap=cm.copper)
ax.plot_surface(x9, y9-0.2, z9,alpha=0.99, cmap=cm.copper)
ax.plot_surface(x11, y11+0.5, z11,alpha=0.99, cmap=cm.copper)
ax.plot_surface(x12, y12+1, z12,alpha=0.99, cmap=cm.copper)



z1 = np.linspace(-1.5, 1.5, 100)
y1 = np.linspace(-1.5, 1.5, 100)
z11, y11 = np.meshgrid(z1, y1)

ax.plot_surface(1E8,y11,z11,color='g', alpha=0.2) 
#ax.text(1E8, 0, 1.5, '$x=10^8 m$', (1,1,0),ha='left',va='bottom')
'''
ax.text(2E8,-1,0.2, '8) Fiberglass', (1,0,0),ha='left',va='bottom')
ax.text(2E8,-0.5,0.2, '9) Kevlar', (1,0,0),ha='left',va='bottom')
ax.text(2E8,0.5,0.2, '11) Carbon Nanotube', (1,0,0),ha='left',va='bottom')
ax.text(2E8,1,0.2, '12) Carbon Nanotube', (1,0,0),ha='left',va='bottom')
'''

fig.savefig("Cone.png", dpi=100,transparent = False)

plt.show()