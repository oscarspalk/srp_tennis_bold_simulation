import numpy
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot, cm

nx = 1001
ny = 1001
nt = 10
nit = 50 
c = 1

length_y = 10
length_x = 50

dx = length_x / (nx - 1)
dy = length_y / (ny - 1)
x = numpy.linspace(0, length_x, nx)
y = numpy.linspace(0, length_y, ny)
X, Y = numpy.meshgrid(x, y)

fig = pyplot.figure(figsize = (11,7), dpi=100)

scaling = 1

u = numpy.zeros((ny, nx), )
un = numpy.zeros((ny, nx))

v = numpy.zeros((ny, nx))
vn = numpy.zeros((ny, nx))

p = numpy.ones((ny, nx))
pn = numpy.ones((ny, nx))

b = numpy.zeros((ny, nx))

ax = fig.add_subplot(projection='3d')

xs =X[::scaling]
ys = Y[::scaling]
ps = p[::scaling]


def circle(x, y, r, arr, val):
    radiusSquared = r**2
    for i in range(x-r, x+r):
        for j in range(y-r, y+r):
            distanceToCenter = (i-x)**2+(j-y)**2
            if radiusSquared > distanceToCenter:
                # we are inside the circle
                arr[i,j] = val
         
            
circle(500, 500, 150, p, 2)

surf = ax.plot_surface(xs,ys,ps ,cmap=cm.viridis,
            linewidth=0, antialiased=False)


ax.set_zlabel("pressure")
ax.set_xlabel("x")
ax.set_ylabel("y")


pyplot.show()