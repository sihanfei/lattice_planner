# lattice_planner
lattice plan algrithm with python

# theory 常规车辆(非常规车辆的运动学方程不同，曲率计算方法不同)
1.d-s : d = f5(s) 5阶曲线 d = a0 + a1 s^1 + a2 s^2 + a3 s^3 + a4 s^4 + a5 s^5
对于起始点，取s=0，d=d0，车辆前轮转角theta，则曲率为 tan(theta)/l(l为前后轮轴距)，曲率导数为 (车轮转角速度)omiga*ctan(theta)/l
有 d0 = a0
   d0'(曲率) = a1 
   d0''(曲率导数) = 2 a2

低于结束点，取s=S，d=d1，且令
有 d1 = a0 + a1 S^1 + a2 S^2 + a3 S^3 + a4 S^4 + a5 S^5
   d1' = a1 + 2 a2 S^1 + 3 a3 S^2 + 4 a4 S^3 + 5 a5 S^4
   d1'' = 2 a2 + 6 a3 S^1 + 12 a4 S^2 + 20 a5 S^3

