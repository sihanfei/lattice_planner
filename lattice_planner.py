import numpy as np
import matplotlib.pyplot as plt
from sympy import diff, symbols
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


def QuinticSL(slist, k):
    l = k[0] + k[1] * slist**1 + k[2] * slist**2 + k[3] * slist**3 + k[
        4] * slist**4 + k[5] * slist**5
    return l


def SLFunction(s):
    """
    input: s
    return: s_matrix
    """
    return [[1.0, s, s**2, s**3, s**4, s**5],
            [0.0, 1.0, 2.0 * s, 3.0 * s**2, 4.0 * s**3, 5.0 * s**4],
            [0.0, 0.0, 2.0, 6.0 * s, 12.0 * s**2, 20.0 * s**3]]


def QuinticSLSolve(**kwargs):
    s0 = 1.0
    if 's0' in kwargs:
        s0 = kwargs['s0']

    l0 = 0.0
    if 'l0' in kwargs:
        l0 = kwargs['l0']

    c0 = 0.0
    if 'c0' in kwargs:
        c0 = kwargs['c0']

    k0 = 0.0
    if 'k0' in kwargs:
        k0 = kwargs['k0']

    s1 = 30.0
    if 's1' in kwargs:
        s1 = kwargs['s1']

    l1 = 0.0
    if 'l1' in kwargs:
        l1 = kwargs['l1']

    c1 = 0.0
    if 'c1' in kwargs:
        c1 = kwargs['c1']

    k1 = 0.0
    if 'k1' in kwargs:
        k1 = kwargs['k1']

    kS0 = SLFunction(s0)

    kS1 = SLFunction(s1)

    kS = np.mat(np.vstack((kS0, kS1)))

    rS = np.mat([l0, c0, k0, l1, c1, k1]).T

    return np.linalg.solve(kS, rS)


def SLTraject(order, *args):
    """
    input:
        1.方程的阶数
        2.单个sl控制点参数,组成方式 s/l/dl/ddl
    output:
        [s^i] x k(nx1) = l'(nx1)
    """
    print(len(args))
    for v in args:
        print(v)
    pass


def slNorder(order):
    x = symbols('x', real=True)
    return [x**i for i in range(order)]


def solveNOrderFunction(*args):
    """
    根据输入参数自动生成n阶方程组,并求解方程系数
    输入系数组成形式: y=f(x)
        [[x0,y0,dy0,ddy0...],[xi,yi,dyi,ddyi...]...]
    输出
        [k0,k1...kn], order
    """
    order = 0
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

    print(matrix_x)
    print(matrix_y)

    coef = np.linalg.solve(matrix_x, matrix_y)  # 系数,列向量
    # print('value x coef', np.dot(value_s, coef))
    return coef, order


def solveNOrderLS(*args):
    """
    给出一系列的ls控制点参数, l = f(s)
    input:[(s0,l0,dl0,ddl0),(s1,l1,dl1,ddl1),(s2,l2,dl2),(s3,l3)]
        每个控制点的参数可以不是4个,但是需要按照 s/l/dl/ddl的顺序来给
    output:
        求解出来的sl方程的系数以及阶数
    """
    order = sum([(len(point_parameter) - 1)
                 for point_parameter in args])  # 因为阶数只和l/dl/ddl相关,所以要减掉s的数量
    value_l = np.zeros((order, 1))  # l/dl/ddl 列向量
    index_l = 0  # 列向量索引
    value_s = np.zeros((order, order))  # nxn的s矩阵

    for point in args:
        s = point[0]  # s的值
        """
        [[s^0   s^1     s^2     s^3 ... s^order],
         [0     1s^0    2s^1    3s^2... order s^(order-1)],
         ....]]
        """
        x = symbols('x', real=True)
        func = [x**i for i in range(order)
                ]  # [s^0   s^1     s^2     s^3 ... s^order]
        for value in point[1:]:
            value_l[index_l, 0] = value
            value_s[index_l, :] = [
                func_i.evalf(subs={x: s}) for func_i in func
            ]  # 求[s^0   s^1     s^2     s^3 ... s^order]实际值
            index_l += 1
            func = [diff(func_i, x, 1) for func_i in func]  # 求1阶导数表达式
        pass
    pass
    # print('value_l', value_l)
    coef = np.linalg.solve(value_s, value_l)  # 系数,列向量
    # print('value x coef', np.dot(value_s, coef))
    return coef, order


def NOrderSL(s, coef):
    l = []
    order = len(coef)
    for si in s:
        li = sum([coef[i] * si**i for i in range(order)])
        l.extend(li)
    return l


def NOrderFunction(xdata, coef):
    ydata = []
    order = len(coef)
    for x in xdata:
        y = sum([coef[i] * x**i for i in range(order)])
        ydata.extend(y)
    return ydata


if __name__ == "__main__":
    s0 = (501.0, 0.0, 1 / 5.0, 0.0)
    s1 = (510.0, 1.0)
    s2 = (523.0, 1.0)
    s3 = (570.0, 1.0, 1 / 5.0, 0.0)
    coef, order = solveNOrderFunction(s0, s1, s2, s3)
    s = np.linspace(501.0, 570.0, 5000)

    l = NOrderFunction(s, coef)
    # print(l)
    plt.plot(s, l)
    plt.plot([s0[0]], [s0[1]], 'o')
    plt.plot([s1[0]], [s1[1]], 'o')
    plt.plot([s2[0]], [s2[1]], 'o')
    plt.plot([s3[0]], [s3[1]], 'o')

    coef, order = solveNOrderLS(s0, s1, s2, s3)
    l = NOrderSL(s, coef)
    plt.plot(s, l, '-')

    plt.show()
