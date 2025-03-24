from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
import math
import numexpr as ne
dt = 0.1


def rk4(x: float, y: float, differential_equation: Callable[[float, float], float]):
    k1 = differential_equation(x,y)
    k2 = differential_equation(x + dt/2, y + dt*(k1/2))
    k3 = differential_equation(x + dt/2, y+ dt*(k2/2))
    k4 = differential_equation(x+dt, y+dt*k3)
    y_new = y + (dt/6) * (k1+2*k2+2*k3+k4)
    x_new = x + dt
    return np.array([x_new, y_new])

def euler(x: float, y: float, differential_equation: Callable[[float, float], float]):
    return np.array([x+dt, y+dt*differential_equation(x,y)])

def diff(x: float, y:float) -> float:
    return y+x


iterations = 200

starting_point = np.array([0,0.1])

runge_values = np.empty((iterations+1, 2))
euler_values = np.empty((iterations+1, 2))

runge_values[0] = starting_point
euler_values[0] = starting_point

for i in range(1,iterations):
    last = i -1
    last_runge_point = runge_values[last]
    r_x = last_runge_point[0]
    r_y = last_runge_point[1]

    last_euler_point = euler_values[last]
    e_x = last_euler_point[0]
    e_y = last_euler_point[1]
    
    rk4_point = rk4(r_x, r_y, differential_equation=diff)
    euler_point = euler(e_x, e_y, differential_equation=diff)
    runge_values[i] = rk4_point
    euler_values[i] = euler_point

fig, ax = plt.subplots()

dts = runge_values[:,0]

ys = ne.evaluate("1.1*exp(dts)-dts-1")


ax.plot(dts, runge_values[:, 1], 'b', label='RK4')
ax.plot(dts, euler_values[:,1], 'r', label='Euler')
ax.plot(dts, ys, 'g', label='Sand funktion')


ax.set_xlim(0,4)
ax.set_ylim(0, 30)
ax.legend()
plt.show()