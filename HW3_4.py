import numpy as np
import matplotlib.pyplot as plt
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import copy

# derivative
K = 900
def dR(R,F):
    return 2*R*(1-R/K)-0.01*R*F #2*R-0.01*R*F

def dF(R,F):
    return -F+0.01*R*F

def AB2(i,h,R_s,F_s):
    R_next = R_s[i] + 3 / 2 * h * dR(R_s[i], F_s[i]) - 1 / 2 * h * dR(R_s[i - 1], F_s[i - 1])
    F_next = F_s[i] + 3 / 2 * h * dF(R_s[i], F_s[i]) - 1 / 2 * h * dF(R_s[i - 1], F_s[i - 1])
    return R_next, F_next

def LF(i,h,R_s,F_s):
    R_next = R_s[i] + 3 / 2 * h * dR(R_s[i], F_s[i]) - 1 / 2 * h * dR(R_s[i - 1], F_s[i - 1])
    F_next = F_s[i] + 3 / 2 * h * dF(R_s[i], F_s[i]) - 1 / 2 * h * dF(R_s[i - 1], F_s[i - 1])
    return R_next, F_next

# time step(h) & total number(N)
h = 0.005
t_last = 100
N = int(t_last/h)

# data base
t_s = [0]
R_0 = [300]
F_0 = [150]


# initial scheme beginning
## 1) Explicit Euler
R_s_AB2_initEE = copy.deepcopy(R_0)
F_s_AB2_initEE = copy.deepcopy(F_0)
R_s_LF_initEE = copy.deepcopy(R_0)
F_s_LF_initEE = copy.deepcopy(F_0)

R_1_EE = R_0[0] + h * dR(R_0[0],F_0[0])
F_1_EE = F_0[0] + h * dF(R_0[0],F_0[0])

R_s_AB2_initEE.append(R_1_EE)
F_s_AB2_initEE.append(F_1_EE)
R_s_LF_initEE.append(R_1_EE)
F_s_LF_initEE.append(F_1_EE)


## 2) RK2
R_s_AB2_initRK2 = copy.deepcopy(R_0)
F_s_AB2_initRK2 = copy.deepcopy(F_0)
R_s_LF_initRK2 = copy.deepcopy(R_0)
F_s_LF_initRK2 = copy.deepcopy(F_0)

k1_R = R_0[0] + h / 2 * dR(R_0[0],F_0[0])
k1_F = F_0[0] + h / 2 * dF(R_0[0],F_0[0])
R_1_RK2 = R_0[0] + h * dR(k1_R,k1_F)
F_1_RK2 = F_0[0] + h * dF(k1_R,k1_F)

R_s_AB2_initRK2.append(R_1_RK2)
F_s_AB2_initRK2.append(F_1_RK2)
R_s_LF_initRK2.append(R_1_RK2)
F_s_LF_initRK2.append(F_1_RK2)

t_s.append(t_s[0]+h)

# Iteration
for i in range(1,N): # i = 1,2,...,N

    # AB2
    R_next_AB2_initEE, F_next_AB2_initEE = AB2(i,h,R_s_AB2_initEE,F_s_AB2_initEE)
    R_next_AB2_initRK2, F_next_AB2_initRK2 = AB2(i, h, R_s_AB2_initRK2, F_s_AB2_initRK2)

    # LF
    R_next_LF_initEE, F_next_LF_initEE = AB2(i,h,R_s_LF_initEE,F_s_LF_initEE)
    R_next_LF_initRK2, F_next_LF_initRK2 = AB2(i,h,R_s_LF_initRK2,F_s_LF_initRK2)

    # append solutions
    R_s_AB2_initEE.append(R_next_AB2_initEE)
    F_s_AB2_initEE.append(F_next_AB2_initEE)
    R_s_AB2_initRK2.append(R_next_AB2_initRK2)
    F_s_AB2_initRK2.append(F_next_AB2_initRK2)

    R_s_LF_initEE.append(R_next_LF_initEE)
    F_s_LF_initEE.append(F_next_LF_initEE)
    R_s_LF_initRK2.append(R_next_LF_initRK2)
    F_s_LF_initRK2.append(F_next_LF_initRK2)
    t_s.append(t_s[i]+h)
    print()


# Scheme 1: AB2 - EE
# Plot: (F,R)
fig, ax = plt.subplots()
plt.title('AB2 - EE')
ax.plot(F_s_AB2_initEE, R_s_AB2_initEE, lw=0.5)
plt.xlabel('Foxes', fontsize=15)
plt.ylabel('Rabbits', fontsize=15)

# Plot: time ~ F,R
fig, axes = plt.subplots(nrows=2, ncols=1)
axes[0].set_title('AB2 - EE')
axes[0].plot(t_s, R_s_AB2_initEE, label='Rabbit')
axes[1].plot(t_s, F_s_AB2_initEE, label='Foxes')
plt.xlabel('time', fontsize=15)
for i in range(0, len(axes)):
    axes[i].legend(loc='upper right')
fig.tight_layout()


# Scheme 2: AB2 - RK2
# Plot: (F,R)
fig, ax = plt.subplots()
ax.plot(F_s_AB2_initRK2, R_s_AB2_initRK2, lw=0.5)
plt.title('AB2 - RK2')
plt.xlabel('Foxes', fontsize=15)
plt.ylabel('Rabbits', fontsize=15)

# Plot: time ~ F,R
fig, axes = plt.subplots(nrows=2, ncols=1)
axes[0].set_title('AB2 - RK2')
axes[0].plot(t_s, R_s_AB2_initRK2, label='Rabbit')
axes[1].plot(t_s, F_s_AB2_initRK2, label='Foxes')
plt.xlabel('time', fontsize=15)
for i in range(0, len(axes)):
    axes[i].legend(loc='upper right')
fig.tight_layout()


# Scheme 3: LF - EE
# Plot: (F,R)
fig, ax = plt.subplots()
ax.plot(F_s_LF_initEE, R_s_LF_initEE, lw=0.5)
plt.title('LF - EE')
plt.xlabel('Foxes', fontsize=15)
plt.ylabel('Rabbits', fontsize=15)

# Plot: time ~ F,R
fig, axes = plt.subplots(nrows=2, ncols=1)
axes[0].set_title('LF - EE')
axes[0].plot(t_s, R_s_LF_initEE, label='Rabbit')
axes[1].plot(t_s, F_s_LF_initEE, label='Foxes')
plt.xlabel('time', fontsize=15)
for i in range(0, len(axes)):
    axes[i].legend(loc='upper right')
fig.tight_layout()


# Scheme 4: LF - RK2
# Plot: (F,R)
fig, ax = plt.subplots()
plt.title('LF - RK2')
ax.plot(F_s_LF_initRK2, R_s_LF_initRK2, lw=0.5)
plt.xlabel('Foxes', fontsize=15)
plt.ylabel('Rabbits', fontsize=15)

# Plot: time ~ F,R
fig, axes = plt.subplots(nrows=2, ncols=1)
axes[0].set_title('LF - RK2')
axes[0].plot(t_s, R_s_LF_initRK2, label='Rabbit')
axes[1].plot(t_s, F_s_LF_initRK2, label='Foxes')
plt.xlabel('time', fontsize=15)
for i in range(0, len(axes)):
    axes[i].legend(loc='upper right')
fig.tight_layout()


plt.show()