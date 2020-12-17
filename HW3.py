from math import *
import numpy as np
import matplotlib.pyplot as plt

# definition of function
def y_diff(t,y):
    return -0.5*tanh(t)*y+exp(-0.5*sqrt(t))

# 1) Explicit method
def EE(n, h,yn, func):
    return yn+h*func(n*h,yn)

# 2) Implicit method
def IE(n, h, yn):
    return 1/(1+0.5*h*tanh((n+1)*h))*yn

# 3) TR method
def TR(n,h,yn):
    return (1-1/4*tanh(n*h))/(1+1/4*tanh((n+1)*h))*yn

# 4) RK2
def RK2(n,h,yn,func):
    prd = yn+h/2*func(yn,n*h)
    return yn+h*func(prd,0.5*h*(2*n+1))

y_init = 1
all_y_EEs = []
all_y_IEs = []
all_y_TRs = []
all_y_RK2 = []

h_s = [0.5, 3.9, 5] # 0.5, 3.9, 5
all_t = []
for i in range(0,len(h_s)):
    y_EE = [y_init]
    y_IE = [y_init]
    y_TR = [y_init]
    y_RK2 = [y_init]

    t = [0]
    h = h_s[i]
    N = floor(30/h)
    for n in range(0,N): # n = 0,1,2,...,N-1
        tn_post = (n+1)*h
        t.append(tn_post)

        # 1) EE
        yn = y_EE[n]
        yn_post = EE(n,h,yn,y_diff)
        y_EE.append(yn_post)
        # 2) IE
        yn = y_IE[n]
        yn_post = IE(n,h,yn)
        y_IE.append(yn_post)
        # 3) TR
        yn = y_TR[n]
        yn_post = TR(n,h,yn)
        y_TR.append(yn_post)
        # 4) RK2
        yn = y_RK2[n]
        # yn_post = RK2(n,h,yn,y_diff)
        yn_post = TR(n, h, yn)
        y_RK2.append(yn_post)

    all_t.append(t)
    all_y_EEs.append(y_EE)
    all_y_IEs.append(y_IE)
    all_y_TRs.append(y_TR)
    all_y_RK2.append(y_RK2)

fig1, ax1 = plt.subplots()
plt.title('Explicit Method')
ax1.plot(all_t[0],all_y_EEs[0], label='h=0.5')
ax1.plot(all_t[1],all_y_EEs[1], label='h=3.9')
ax1.plot(all_t[2],all_y_EEs[2], label='h=5')
ax1.legend()

fig2, ax2 = plt.subplots()
plt.title('Implicit Method')
ax2.plot(all_t[0],all_y_IEs[0], label='h=0.5')
ax2.plot(all_t[1],all_y_IEs[1], label='h=3.9')
ax2.plot(all_t[2],all_y_IEs[2], label='h=5')
ax2.legend()

fig3, ax3 = plt.subplots()
plt.title('Trapezoidal Method')
ax3.plot(all_t[0],all_y_TRs[0], label='h=0.5')
ax3.plot(all_t[1],all_y_TRs[1], label='h=3.9')
ax3.plot(all_t[2],all_y_TRs[2], label='h=5')
ax3.legend()

fig4, ax4 = plt.subplots()
plt.title('RK2 Method')
ax4.plot(all_t[0],all_y_RK2[0], label='h=0.5')
ax4.plot(all_t[1],all_y_RK2[1], label='h=3.9')
ax4.plot(all_t[2],all_y_RK2[2], label='h=5')
ax4.legend()

# plt.show()

# exact_t = np.linspace(0,30,100)
# exact_y = []
# for temp_t in exact_t:
#     exact_y.append(np.sqrt(2/(exp(temp_t)+exp(-(temp_t)))))
# fig5, ax5 = plt.subplots()
# plt.title('Exact solution')
# ax5.plot(exact_t,exact_y)
# plt.show()