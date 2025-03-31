from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
import math
import numexpr as ne
import scipy.integrate as integrate

pi = 3.14159
d = 0.067 # ball diameter
m = 0.058 # ball mass 
g = 9.81 # gravity
rho = 1.2041 # air density
n = 1800 # spin

i = 1

alpha = (10/180)*pi

tau = alpha

w = (d/2)*(n*6.28/60)

v = 30

dtau = 0.005

x = 0

y0 = 0.8
y = y0

def D():
    C_d = 0.508 + 1/(
        (22.503+ 4.196/(ratio()**2.5))**0.4
    )
    return C_d * ((pi*(d**2))/(8*m*g)) * rho * (v**2)

def M():
    C_l = 1/(2.022+(0.981/ratio()))
    return C_l * ((pi*(d**2))/(8*m*g)) * rho * (v**2)

def ratio():
    return w/v

def y_integrand(integrator):
    return (v**2*math.sin(integrator))/(math.cos(integrator)+M())

def x_integrand(integrator):
    return (v**2*math.cos(integrator))/(
        math.cos(integrator)+M()
    )

def calc_y():
    return y0 - (1/g)* integrate.quad(y_integrand, alpha, tau)[0]

def calc_x():
    return -1/g * integrate.quad(x_integrand, alpha, tau)[0]

def diff(t, v):
    return (
        math.sin(t)+D()
    )/(math.cos(t)+M())*v

points = []
velocities = []

while y > 0:
    y = calc_y()
    x = calc_x()
    points.append([x,y])
    v = v - dtau * diff(tau, v)
    velocities.append(v)
    tau = tau - dtau
    i = i + 1
    print(f"Height is: {y}")


points = np.array(points)
fig, ax = plt.subplots()

#ax.plot(points[:,0], points[:, 1], 'b', label='y')
ax.plot(points[:, 0], np.array(velocities)/10, 'g', label='v')
plt.legend()
plt.show()