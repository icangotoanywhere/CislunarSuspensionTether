# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 17:36:21 2024

@author: Www
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 定义常数和参数
G = 6.67430e-11  # 引力常数 (m^3 kg^-1 s^-2)
rho = 2.7e3  # 小行星密度 (kg/m^3)
omega = 2 * np.pi / (12 * 3600)  # 小行星的角速度 (rad/s)

# 定义动力学方程
def orbital_siphon_dynamics(t, y, params):
    # y 是一个包含位置和速度的向量，形状为 (6,)，其中前3个元素是位置，后3个元素是速度
    # params 包含模型参数，例如质量、弹性系数等
    m = params['mass']
    k = params['spring_constant']
    c = params['damping_coefficient']

    # 计算离心力
    centrifugal_force = omega**2 * y[:3]

    # 计算弹簧力
    F_spring = -k * np.linalg.norm(y[3:] - y[:3])

    # 计算阻尼力
    F_damping = -c * y[3:]

    # 计算加速度
    a = np.zeros(6)
    a[:3] = (centrifugal_force + F_spring + F_damping) / m
    a[3:] = np.zeros(3)  # 假设没有外部力作用在速度上

    return a

# 定义时间积分函数
def integrate_dynamics(t_span, y0, params):
    return solve_ivp(orbital_siphon_dynamics, t_span, y0, args=(params,))

# 初始条件
t_span = (0, 86400)  # 模拟时间跨度，这里是一天
y0 = np.zeros(6)  # 初始位置和速度，例如：位置 (0, 0, 10) 米，速度 (0, 0, 0) m/s
params = {'mass': 1e3, 'spring_constant': 1e5, 'damping_coefficient': 1e3}

# 进行数值积分
solution = integrate_dynamics(t_span, y0, params)

# 绘制结果
fig, axs = plt.subplots(2, 1, figsize=(10, 8))  # 创建一个包含两个子图的图形对象

# 绘制位置随时间的变化
axs[0].plot(solution.t, solution.y[:3, 0], label='Position')  # 只绘制前三个维度，即位置
axs[0].legend()
axs[0].set_xlabel('Time (s)')
axs[0].set_ylabel('Position (m)')

# 绘制速度随时间的变化
axs[1].plot(solution.t, solution.y[3:, 0], label='Velocity')  # 绘制后三个维度，即速度
axs[1].legend()
axs[1].set_xlabel('Time (s)')
axs[1].set_ylabel('Velocity (m/s)')

plt.tight_layout()  # 调整子图间距
plt.show()
