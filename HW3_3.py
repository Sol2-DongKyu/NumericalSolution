import numpy as np
import math
import matplotlib.pyplot as plt

import numpy as np
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
a = 10
b = 8/3
r = 28

def dX(X,Y,Z):
    return -a*X+a*Y

def dY(X,Y,Z):
    return -X*Z+r*X-Y

def dZ(X,Y,Z):
    return X*Y-b*Z

# time step(h) & total number(N)
h = 0.01 # 0.01, 0.1, 1, 10
t_last = 70
N = int(t_last/h)


# data base
t_s = [0]
X_s = [-0.001]
Y_s = [-0.001]
Z_s = [0]
dX_s = [dX(X_s[0],Y_s[0],Z_s[0])]

# Iteration
for i in range(0,N): # i = 0,1,2,...,N-1
    X_i = X_s[i]
    Y_i = Y_s[i]
    Z_i = Z_s[i]

    k1_X = X_i + h / 2 * dX(X_i,Y_i,Z_i)
    k1_Y = Y_i + h / 2 * dY(X_i, Y_i, Z_i)
    k1_Z = Z_i + h / 2 * dZ(X_i, Y_i, Z_i)

    k2_X = X_i + h / 2 * dX(k1_X, k1_Y, k1_Z)
    k2_Y = Y_i + h / 2 * dY(k1_X, k1_Y, k1_Z)
    k2_Z = Z_i + h / 2 * dZ(k1_X, k1_Y, k1_Z)

    k3_X = X_i + h * dX(k2_X, k2_Y, k2_Z)
    k3_Y = Y_i + h * dY(k2_X, k2_Y, k2_Z)
    k3_Z = Z_i + h * dZ(k2_X, k2_Y, k2_Z)

    X_next = X_i + h * (1/6*dX(X_i,Y_i,Z_i) + 1/3*dX(k1_X,k1_Y,k1_Z) + 1/3*dX(k2_X,k2_Y,k2_Z) + 1/6*dX(k3_X,k3_Y,k3_Z))
    Y_next = Y_i + h * (1/6*dY(X_i,Y_i,Z_i) + 1/3*dY(k1_X,k1_Y,k1_Z) + 1/3*dY(k2_X,k2_Y,k2_Z) + 1/6*dY(k3_X,k3_Y,k3_Z))
    Z_next = Z_i + h * (1/6*dZ(X_i,Y_i,Z_i) + 1/3*dZ(k1_X,k1_Y,k1_Z) + 1/3*dZ(k2_X,k2_Y,k2_Z) + 1/6*dZ(k3_X,k3_Y,k3_Z))

    # append solutions
    X_s.append(X_next)
    Y_s.append(Y_next)
    Z_s.append(Z_next)
    t_s.append(t_s[i]+h)
    dX_s.append(dX(X_next ,Y_next ,Z_next))


# Plot
fig, axes = plt.subplots(nrows=3, ncols=1)

axes[0].plot(t_s, X_s, label='X')
axes[1].plot(t_s, Y_s, label='Y')
axes[2].plot(t_s, Z_s, label='Z')
plt.xlabel('time', fontsize=15)
for i in range(0, len(axes)):
    axes[i].legend(loc='upper right')

fig.tight_layout()

# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

# Plot
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(X_s, Y_s, Z_s, lw=0.5)
ax.scatter(X_s[0], Y_s[0], Z_s[0], lw=0.5, marker="X", c='r')
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Lorenz Attractor")

# Plot
fig, ax = plt.subplots()

ax.plot(X_s, dX_s, lw=0.5)
ax.scatter(X_s[0], dX_s[0], lw=0.5, marker="X", c='r')
ax.set_xlabel("X")
ax.set_ylabel("dX/dt")

plt.show()
