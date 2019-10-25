# coding : utf-8

import numpy as np
import math
from sympy import diff, symbols

import matplotlib.pyplot as plt

import lattice_planner as lp


class AdjustmentControlPointStruct(list):
    def __init__(self, x=0.0, y=0.0, theta=0.0):
        super().__init__([x, y, theta])
        self.x = self[0]
        self.y = self[1]
        self.theta = self[2]

    def set_x(self, x):
        self[0] = x
        self.x = x

    def set_y(self, y):
        self[1] = y
        self.y = y

    def set_theta(self, theta):
        self[2] = theta
        self.theta = theta


def makeAdjustmentLine(start,
                       end,
                       min_radius=0,
                       max_angle_velocity=40,
                       step=0.1):
    """
    根据给定点的坐标(x,y)/朝向(东北天),以及车辆最小转弯半径/最大转动角速度,计算连接给定点的缓行道路线
    input:
      start(x,y,theta=angle)
      end(x,y,theta=angle)
      min_radius(m)
      max_angle_velocity(angle/s)
      step(离散距离)
    output:
      xy_data
    theory:
      y0=f(x0)
      dy/dx|x0 = tan(theta)
      所以可以是4阶
        y=k0+k1*x^1+k2*x^2+k3*x^3
    """
    coef, order, error = lp.solveNOrderFunction(
        (start.x, start.y, np.tan(start.theta * np.pi / 180)),
        (end.x, end.y, np.tan(end.theta * np.pi / 180)))
    print(coef, order, error)
    x_data = np.arange(start.x, end.x, step * (1 if start.x < end.x else -1))

    y_data = lp.getNOrderOutput(x_data, coef)
    return x_data, y_data


if __name__ == "__main__":
    CP0 = AdjustmentControlPointStruct(0, 0, 0)
    CP1 = AdjustmentControlPointStruct(-10, 10, 90)
    x, y = makeAdjustmentLine(CP0, CP1)
    # print(x, y)
    plt.plot(x, y)
    plt.show()
