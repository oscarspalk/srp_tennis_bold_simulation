from typing import Callable
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import math
import numexpr as ne
import mpl_toolkits.mplot3d.art3d as art3d


f_b = np.array([0, 0.0017]) # newton
f_t = np.array([0, -0.55])  # Newton tennis ball mass
spin = 0#53.33 # Hz
dt = 0.01
m = 0.056 # kg
# x is time, y is velocity
def rk4(x: float, y: np.array, differential_equation: Callable[[float, np.array], float]):
    k1 = differential_equation(x,y)
    k2 = differential_equation(x + dt/2, y + dt*(k1/2))
    k3 = differential_equation(x + dt/2, y+ dt*(k2/2))
    k4 = differential_equation(x+dt, y+dt*k3)
    y_new = y + (dt/6) * (k1+2*k2+2*k3+k4)
    x_new = x + dt
    return (x_new, y_new)

def vec_len(v):
    return math.sqrt(v[0]**2+v[1]**2)

def f_m(v):
    hatted_v = np.array([-v[1], v[0]])
    hat_vec_norm = hatted_v/vec_len(hatted_v)
    return hat_vec_norm * -6.18*math.pow(10, -6)*(spin-4.8)* vec_len(v)**2 # newton

def f_d(v):
    return -9.36*math.pow(10,-4)* v**2 # newton

def diff(x: float, y: np.array) -> float:
    return (f_t+f_b+f_d(y)+f_m(y))/m

iterations = 2000
starting_velocity = np.array([20, 4]) # meters per second

runge_t_values = np.empty(iterations+1)
runge_values = np.empty((iterations+1,2))

x_y = np.empty((iterations+1, 2))

runge_t_values[0] = 0
runge_values[0] = starting_velocity

x_y[0] = [0, 0.9]

break_index = 0

for i in range(1,iterations):
    # udregn nye hastigheder
    last = i -1
    last_runge_point = runge_values[last]
    last_time = runge_t_values[last]
    rk4_point = rk4( last_time,last_runge_point, differential_equation=diff)
    runge_t_values[i] = rk4_point[0]
    runge_values[i] = rk4_point[1]

    # udregn nyt sted
    last_place = x_y[last]
    place_point = last_place + rk4_point[1]*dt
    x_y[i] = place_point

    if(x_y[i][1] < 0):
        break_index = i
        break

dts = runge_t_values[:break_index]

fig = plt.figure(figsize=(11,7))

ax = fig.add_subplot(projection='3d', computed_zorder=False)

fig.tight_layout()

tennis_court = Rectangle((0,0), 23.77, 8.23, color='green',zorder=0)
tennis_net = Rectangle((0,0), 8.23, 1.06, color='brown', zorder=1)
ax.add_patch(tennis_court)
ax.add_patch(tennis_net)
art3d.pathpatch_2d_to_3d(tennis_court, z=0, zdir="z")
art3d.pathpatch_2d_to_3d(tennis_net, z=23.77/2, zdir="x")

ax.plot(x_y[:break_index, 0], np.full_like(dts, 2),x_y[:break_index, 1],  'b', label='Bold',zorder=2)
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_zlabel("z (m)")
ax.set_box_aspect(aspect=(23.77, 8.23, 10))
ax.set_xlim(0,23.77)
ax.set_ylim(0, 8.23)
ax.set_zlim(0,3)


plt.savefig("sim_3d")

plt.show()
fig = plt.figure(figsize=(11,7))

ax = fig.add_subplot()

ax.plot(x_y[:break_index, 0], x_y[:break_index, 1], 'b')
ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")
ax.set_xlim(0, 23.77)
ax.set_ylim(0,3)

tennis_net = Rectangle((23.77/2,0), 0.1, 1.06, color='brown', zorder=1, label='Net')
ax.add_patch(tennis_net)
ax.legend()
def annot_max(x,y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    text= "x={:.3f}, y={:.3f}".format(xmax, ymax)
    if not ax:
        ax=plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops=dict(arrowstyle="->",connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data',textcoords="axes fraction",
              arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    ax.annotate(text, xy=(xmax, ymax), xytext=(0.94,0.96), **kw)

annot_max(x_y[:break_index, 0],x_y[:break_index, 1])

plt.savefig("sim_2d")

plt.show()