from ufl import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from dune.ufl import Constant, DirichletBC
from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.fem.plotting import plotPointData

name = 'streamtracer'
file = 'input/'+name+'.{:04d}.png'
N = 2

print("Read images")
images = []
for i in range(1, N+1):
    image = mpimg.imread(file.format(i))
    image = image[:, :, 0]
    images += [image]

Y, X = images[0].shape

print("Construct grid with {}x{} elements".format(X, Y))
grid = structuredGrid([0,0], [X,Y], [X-1,Y-1])
cspace = lagrange(grid, order=1)

c = []
for i in range(N):
    ch = cspace.interpolate(0, name="c{:03d}".format(i))
    ch.as_numpy[:] = np.flip(images[i], 0).flatten()
    c += [ch]

vspace = lagrange(grid, order=1, dimRange=2)
v = TrialFunction(vspace)
vv = TestFunction(vspace)

# Horn-Schunck method: (c_t + \nabla c \cdot \mathbf{v}) \nabla c - \alpha^2 \Delta \mathbf{v} = 0
print("Solve Horn-Schunck method")
A = 0
for i in range(N-1):
    A += ((c[i+1] - c[i]) + dot(grad(c[i]), v)) * dot(grad(c[i]), vv) * dx

alpha = Constant(10, name="alpha")
A += alpha**2 * inner(grad(v), grad(vv)) * dx

x = SpatialCoordinate(vspace)
bottom = DirichletBC(vspace, zero(2), x[1] < 1e-6)
top = DirichletBC(vspace, zero(2), x[1] > Y-1e-6)

scheme = galerkin([A == 0, bottom, top], solver=("suitesparse", "umfpack"))
vh = vspace.function(name="vh")
scheme.solve(vh)

pointdata = {"v": [vh[0], vh[1], 0]}
print("Write output")
for i in range(N):
    pointdata["c"] = c[i]
    grid.writeVTK("output/"+name+"-"+str(i), pointdata=pointdata)

fig = plt.figure(figsize=(8,4.8))
plotPointData(vh[0], figure=fig, gridLines=None, colorbar=None)
plt.savefig("output/"+name+".png", dpi=300)
