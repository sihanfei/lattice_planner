import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols
import math
"""
l = k0 + k1 * s + k2 * s^2 + k3 * s^3 + k4 * s^4 + k5 * s^5
dl/ds(c) =  k1     + 2*k2*s^1 + 3*k3*s^2 + 4*k4*s^3 + 5*k5*s^4 : 曲率c. 限定曲率,实质限定了前后轮转角 : 转弯半径为 轴长/tan(前轮转角theta) : 基于简化自行车模型 : 曲率c = tan(theta)/轴长
ddl/dds(k) =         2*k2     + 6*k3*s^1 + 12*k4*s^2 + 20*k5*s^3 : 曲率变化率k. 限定曲率变化率,实质限定了转角速度 d(theta)/ ds
----
s = j0 + j1 * t + j2 * t^2 + j3 * t^3 + j4 * t^4 + j5 * t^5
ds/dt =  j1     + 2*j2*t^1 + 3*j3*t^2 + 4*j4*t^3 + 5*j5*t^4 : 速度 v
dds/ddt =         2*j2     + 6*j3*t^1 + 12*j4*t^2 + 20*j5*t^3 : 加速度 a
"""

# 起终点在圆弧上,起始状态 (s0=1,l0=0,c0=1/10,k0=0),终点状态(s1=17,l1=0,c1=1/10,k1=0)


def solveNOrderFunction(*args):
    """
    根据输入参数自动生成n阶方程组,并求解方程系数
    输入系数组成形式: y=f(x)
        [[x0,y0,dy0,ddy0...],[xi,yi,dyi,ddyi...]...]
    输出
        [k0,k1...kn], order
    """
    order = 0
    error = 0
    # 与输入参数个数等阶数
    for value_i in args:
        order += len(value_i) - 1

    matrix_x = np.zeros((order, order))
    matrix_y = np.zeros((order, 1))
    index = 0

    for value_i in args:
        x = value_i[0]  # x的值,ls方程中的s,st方程中的t
        for row_i, y_diff_i in enumerate(value_i[1:]):
            # [y0,dy0,...,y1,dy1,...].T
            matrix_y[index, 0] = y_diff_i  # Y
            # [[x0**0   x0**1     x0**2 ... x0**(order-1)       x0**order]
            #  [0     1*x0**0   2*x0**1 ...               order*x0**(order-1)]
            #  [0     0         2*1*x0**0 ..    order*(order-1)*x0**(order-2)]
            #  [x1**0   x1**1 ...]
            #  []]
            for col_i in range(row_i, order):
                if row_i == 0:
                    matrix_x[index, col_i] = x**col_i
                else:
                    matrix_x[index, col_i] = col_i * matrix_x[index -
                                                              1, col_i - 1]
            index += 1
    try:
        coef = np.linalg.solve(matrix_x, matrix_y)  # 系数,列向量
        y = np.dot(matrix_x, coef)
        print('y=', y)
        error = np.sqrt(sum((y - matrix_y)**2))
    # print('value x coef', np.dot(value_s, coef))
    except RuntimeError:
        coef = []
        order = 0
        error = 0
    finally:
        return coef, order, error


def getNOrderOutput(xdata, coef):
    """
    根据输入的xdata数据,求以coef为系数的n阶函数的输出值
    """
    ydata = []
    order = len(coef)
    for x in xdata:
        y = sum([coef[i] * x**i for i in range(order)])
        ydata.extend(y)
    return ydata


if __name__ == "__main__":
    s0 = (-3.2751, 5.055, 3.4671, -4.5613)
    s1 = (-1.5783, 6.3587)
    s2 = (3.1471, -4.5666, 1.429, -9.5802)

    # s0 = (501.0, 0.0, 1 / 5.0, 0.0)
    # s1 = (510.0, 1.0)
    # s2 = (523.0, 1.0)
    # s3 = (570.0, 1.0, 1 / 5.0, 0.0)
    coef, order, error = solveNOrderFunction(s0, s1, s2)
    print('error=', error)
    s = np.linspace(s0[0], s2[0], math.ceil(s2[0] + 1 - s0[0]) * 10)

    l = getNOrderOutput(s, coef)
    # print(l)
    plt.plot(s, l)
    plt.plot([s0[0]], [s0[1]], 'o')
    plt.plot([s1[0]], [s1[1]], 'o')
    plt.plot([s2[0]], [s2[1]], 'o')

    plt.show()
