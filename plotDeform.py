import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
#from copy import deepcopy

def MaxEnt( params, x, Xa):
    n = len(params);
    na = len(Xa);
    nx = len(x);
    p = np.zeros(na);
    z = 0
    beta = np.array(params[nx:]).reshape((nx, nx));
    for i, xa in enumerate(Xa):
        h = beta.dot(x-xa);
        #h = np.fabs(h-np.rint(h));
        x2 = -np.sum(np.square(h));
        x1 = np.dot(params[0:nx], x-xa);
        z += np.exp(x2+x1)
        p[i]= np.exp(x2+x1);
    return p/z;

fnames = ["MaxEntAtoms.txt", "MaxEntNodes.txt", "MaxEntParams.txt"];
data = dict.fromkeys(range(len(fnames)),None);
for i, v in enumerate(fnames):
    d = np.loadtxt(v);
    data[i] = d;

X = data[0]
Y = data[1]
L = data[2]

N = len(X)
D = 2;
n = len(Y)

x0 = np.array([d[0:D] for d in X]).reshape((N,D))
a0 = [d[D:] for d in X] 
y0 = [d[0:D] for d in Y]
b0 = [d[D:] for d in Y]

hull = ConvexHull(x0);

xu0 = [p[0::D] for p in a0]
xv0 = [p[1::D] for p in a0]

yu0 = [p[0::D] for p in b0]
yv0 = [p[1::D] for p in b0]


fnames = ["MaxEntAtoms1.txt", "MaxEntNodes1.txt", "MaxEntParams1.txt"];
data1 = dict.fromkeys(range(len(fnames)),None);
for i, v in enumerate(fnames):
    d = np.loadtxt(v);
    data1[i] = d;

x = data1[0]
y = data1[1]
l = data1[2]

x1 = np.array([d[0:D] for d in x]).reshape((N,D))
a1 = [d[D:] for d in x] 
y1 = [d[0:D] for d in y]
b1 = [d[D:] for d in y]

hull1 = ConvexHull(x1);

xu1 = [p[0::D] for p in a1]
xv1 = [p[1::D] for p in a1]

yu1 = [p[0::D] for p in b1]
yv1 = [p[1::D] for p in b1]


fig, (ax0, ax1) = plt.subplots(2, sharex=True, sharey=True)

for simplex in hull.simplices:
    ax0.plot(x0[simplex, 0], x0[simplex, 1], 'k-o')

ax0.quiver([d[0] for d in x0], [d[1] for d in x0], [u[0] for u in xu0], [u[1] for u in xu0], units = 'xy', scale = 1, color = 'k', width = 5e-3, );
ax0.quiver([d[0] for d in x0], [d[1] for d in x0], [v[0] for v in xv0], [v[1] for v in xv0], units = 'xy', scale = 1, color = 'k', width = 5e-3, );

ax0.scatter([d[0] for d in x0], [d[1] for d in x0], color='k');

ax0.quiver([d[0] for d in y0], [d[1] for d in y0], [u[0] for u in yu0], [u[1] for u in yu0], units = 'xy', scale = 1, color = 'b', width = 5e-3, )
ax0.quiver([d[0] for d in y0], [d[1] for d in y0], [v[0] for v in yv0], [v[1] for v in yv0], units = 'xy', scale = 1, color = 'b', width = 5e-3, )

ax0.scatter([d[0] for d in y0], [d[1] for d in y0], color='b');

ax0.grid(True)

for simplex in hull1.simplices:
    ax1.plot(x1[simplex, 0], x1[simplex, 1], 'k-o')

ax1.quiver([d[0] for d in x1], [d[1] for d in x1], [u[0] for u in xu1], [u[1] for u in xu1], units = 'xy', scale = 1, color = 'k', width = 5e-3, minshaft = 1, minlength = 1)
ax1.quiver([d[0] for d in x1], [d[1] for d in x1], [v[0] for v in xv1], [v[1] for v in xv1], units = 'xy', scale = 1, color = 'k', width = 5e-3, minshaft = 1, minlength = 1)

ax1.scatter([d[0] for d in x1], [d[1] for d in x1], color='k');

ax1.quiver([d[0] for d in y1], [d[1] for d in y1], [u[0] for u in yu1], [u[1] for u in yu1], units = 'xy', scale = 1, color = 'g', width = 5e-3, minshaft = 1, minlength = 1)
ax1.quiver([d[0] for d in y1], [d[1] for d in y1], [v[0] for v in yv1], [v[1] for v in yv1], units = 'xy', scale = 1, color = 'g', width = 5e-3, minshaft = 1, minlength = 1)

ax1.scatter([d[0] for d in y1], [d[1] for d in y1], color='g');

ax1.grid(True)


plt.show();

