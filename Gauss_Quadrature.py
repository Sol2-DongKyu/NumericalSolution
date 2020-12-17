import pandas as pd
import os
from math import *
import numpy as np
dataPath = 'D:\\NaverCloud\\대학원\\20년 2학기 수업\\수치해석\\HW2'

def f(x):
    return 2*pow(x,2)*exp(2*x) + 1/(pow(x-1,2)+0.008) + 1/sqrt(4-pow(x-1,2)) + log(4*x+1)

## N = 8
df_n8 = pd.read_csv(os.path.join(dataPath, 'GQtable_n8.csv'),
                     header=None)
## N = 16
df_n16 = pd.read_csv(os.path.join(dataPath, 'GQtable_n16.csv'),
                     header=None)
## N = 32
df_n32 = pd.read_csv(os.path.join(dataPath, 'GQtable_n32.csv'),
                     header=None)
N = 32
df = df_n32

w = df[1]
absc = df[2]

I = 0
for i in range(0,N):
    print('absc:{}, weight:{}'.format(absc[i],w[i]))
    I += 0.5*f(0.5+0.5*absc[i])*w[i]

error_ = abs((I-21.2946)/21.2946)*100
print('Integral: {}, Error:{}'.format(I,log2(error_)))