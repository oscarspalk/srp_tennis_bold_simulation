import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

## Creates a circle at a given point in a 2D array
def circle(x, y, r, arr, val):
    radiusSquared = r**2
    for i in range(x-r, x+r):
        for j in range(y-r, y+r):
            distanceToCenter = (i-x)**2+(j-y)**2
            if radiusSquared > distanceToCenter:
                # we are inside the circle
                arr[i,j] = val

wind_speed = 5.

def build_up_b(rho, dt, dx, dy, u, v):
    b = numpy.zeros_like(u)
    b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) +
                                      (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))-
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    


   
    
    # Static BC Pressure @ cylinder, p = 0

    return b

def pressure_poisson_periodic(p, dx, dy):
    pn = numpy.empty_like(p)
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) /
                         (2 * (dx**2 + dy**2)) -
                         dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        
        # Wall boundary conditions, pressure
        p[-1, :] =p[-2, :]  # dp/dy = 0 at y = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
    
        # Static BC Pressure @ cylinder, p = 0
        circle(500, 500, 150, p, 0)

    return p

##variable declarations
nx = 1001
ny = 1001
nt = 100
nit = 50 
c = 1

length_y = 500
length_x = 500

dx = length_x / (nx - 1)
dy = length_y / (ny - 1)
x = numpy.linspace(0, length_x, nx)
y = numpy.linspace(0, length_y, ny)
X, Y = numpy.meshgrid(x, y)


##physical variables
rho = 1
nu = .1
F = 0

sigma = 0.05
dt =(dx*sigma)/wind_speed

#initial conditions
u = numpy.full((ny, nx), wind_speed, dtype=numpy.float64)
un = numpy.zeros((ny, nx),dtype=numpy.float64)

v = numpy.zeros((ny, nx),dtype=numpy.float64)
vn = numpy.zeros((ny, nx),dtype=numpy.float64)

p = numpy.ones((ny, nx),dtype=numpy.float64)
pn = numpy.ones((ny, nx),dtype=numpy.float64)

b = numpy.zeros((ny, nx),dtype=numpy.float64)

udiff = 1
stepcount = 0

print(f"dt: {dt}")

while stepcount < nt:
    un = u.copy()
    vn = v.copy()

    b = build_up_b(rho, dt, dx, dy, u, v)
    p = pressure_poisson_periodic(p, dx, dy)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx * 
                    (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * 
                    (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                     dt / (2 * rho * dx) * 
                    (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                     nu * (dt / dx**2 * 
                    (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                     dt / dy**2 * 
                    (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + 
                     F * dt)

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                     un[1:-1, 1:-1] * dt / dx * 
                    (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     vn[1:-1, 1:-1] * dt / dy * 
                    (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                     dt / (2 * rho * dy) * 
                    (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                     nu * (dt / dx**2 *
                    (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                     dt / dy**2 * 
                    (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # Wall BC: u,v = 0 @ y = 0,2
    u[0, :] = 0
    u[-1, :] = 0
    v[0, :] = 0
    v[-1, :]=0
    
    # Static BC Pressure @ cylinder, u,v = 0
    circle(500, 500, 150, u, 0)
    circle(500, 500, 150, v, 0)

    udiff = (numpy.sum(u) - numpy.sum(un)) / numpy.sum(u)
    stepcount += 1
    print(f"at step {stepcount}")

scaling = 20
fig = pyplot.figure(figsize = (7,7), dpi=100)
pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)  
pyplot.colorbar()
pyplot.streamplot(X,Y,u,v)

pyplot.show()