from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
import math
import numexpr as ne
import scipy.integrate as integrate
import scipy.special as special

pi = 3.14159
d = 0.067 # ball diameter
m = 0.058 # ball mass 
g = 9.82 # gravity
rho = 1.2041 # air density
n = 3000 # spin
alpha = (10/180)*pi
w = (d/2)*(n*2*pi/60)

dt = 0.005


def D(v):
    C_d = 0.508 + 1/(
        (22.503+ 4.196/((w/v)**2.5))**0.4
    )
    return C_d * ((pi*(d**2))/(8*m*g)) * rho * (v**2)

def M(v):
    C_l = 1/(2.022+(0.981/(w/v)))
    return C_l * ((pi*(d**2))/(8*m*g)) * rho * (v**2)

def diff(t: float, v: float) -> float:
    return ((math.sin(t)+D(v))/(math.cos(t)+M(v))) * v

def rk4(x: float, y: float, differential_equation: Callable[[float, float], float]):
    k1 = differential_equation(x,y)
    k2 = differential_equation(x + dt/2, y + dt*(k1/2))
    k3 = differential_equation(x + dt/2, y+ dt*(k2/2))
    k4 = differential_equation(x+dt, y+dt*k3)
    y_new = y + (dt/6) * (k1+2*k2+2*k3+k4)
    x_new = x + dt
    return np.array([x_new, y_new])

iterations = 2000

starting_point = np.array([0,1])
velocities = np.empty(iterations+1)

velocities[0] = 30 # start velocity

runge_values = np.empty((iterations+1, 2))


def calc_x( t,v):
    return -(1/g)*integrate.quad(
        lambda t : ((v**2*math.cos(t))/(math.cos(t)+M(v)))
         ,alpha, t)

def calc_y(t,v):
    return starting_point[0][1] - 1/g * integrate.quad(
        lambda t : ((v**2*math.sin(t))/(math.cos(t)+M(v))), alpha, t
    )

for i in range(1,iterations):
    last = i -1
    last_runge_point = runge_values[last]
    last_velocity = velocities[last]
    r_x = last_runge_point[0]
    r_y = last_runge_point[1]

    
    rk4_point = rk4(iterations*dt, last_velocity, differential_equation=diff) # gives a new time and velocity
    time = rk4_point[0]
    velocity = rk4_point[1]
    velocities[i] = velocity
    runge_values[i] = [calc_x(time, velocity), calc_y(time,velocity)]

fig, ax = plt.subplots()

ax.plot(runge_values[:,0], runge_values[:, 1], 'b', label='RK4')

plt.show()