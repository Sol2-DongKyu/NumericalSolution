from math import *
import numpy as np
import matplotlib.pyplot as plt

global Xs, Fs


def f(x):
    return 2 * pow(x, 2) * exp(2*x) +1/( pow((x-1),2) +0.008 ) + \
        1/(sqrt(4-pow((x-1),2))) + log(4*x+1)
    # return 5*pow(x,3)+2*pow(x,2)

def adaptive_simpson( f, a, b, tol ):
    """
    Evaluates the integral of f(x) on [a,b].

    USAGE:
        s = adaptive_simpson( f, a, b, tol )

    INPUT:
        f         - function to integrate
        a, b      - left and right endpoints of interval of integration
        tol       - desired upper bound on allowable error

    OUTPUT:
        float     - the value of the integral

    NOTES:
        Integrates the function f(x) on [a,b] with an error bound of
        given by tol using an adaptive Simpson's rule to approximate
        the value of the integral.  Notice that this is not very
        efficient -- it is recursive and recomputes many function
        values that have already been computed.  The code in this
        function is meant to be clear and explain the ideas behind
        adaptive schemes rather than to be an efficient implementation
        of one.

    AUTHOR:
        Jonathan Senning <jonathan.senning@gordon.edu>
        Gordon College, Department of Mathematics and Computer Science
        Written in Octave February 19, 1999
        Converted to Python March 5, 2008
    """

    # Theory says the factor to multiply the tolerance by should be 15, but
    # since that assumes that the fourth derivative of f is fairly constant,
    # we want to be a bit more conservative...

    tol_factor = 15.0

    h = 0.5 * ( b - a )

    x0 = a
    x1 = a + 0.5 * h
    x2 = a + h
    x3 = a + 1.5 * h
    x4 = b

    f0 = f( x0 )
    f1 = f( x1 )
    f2 = f( x2 )
    f3 = f( x3 )
    f4 = f( x4 )

    s0 = (2.0 * h) * ( f0 + 4.0 * f2 + f4 ) / 6.0
    s1 = (2.0 * h) * ( f0 + 4.0 * f1 + 2.0 * f2 + 4.0 * f3 + f4 ) / 12.0

    if abs( s0 - s1 ) >= tol_factor * tol:
        s = adaptive_simpson( f, x0, x2, 0.5 * tol ) + \
            adaptive_simpson( f, x2, x4, 0.5 * tol )
    else:
        s = s1 + ( s1 - s0 ) / 15.0
        for i in [x0,x1,x2,x3,x4]:
            Xs.append(i)
        for i in [f0,f1,f2,f3,f4]:
            Fs.append(i)
    return s


# print('{}\n{}'.format(f(2),f(4)))
count = 1;
for tol in np.power(0.1,[i for i in range(1,5)]): #[0.1~0.0001]
    Xs = []
    Fs = []
    ax = []
    print('tolerance: {}'.format(tol))
    ans1 = adaptive_simpson( f, 0, 1, tol )
    print('Integral: {}\n'.format(ans1))
    print('Xs : {}\n'.format(len(Xs)))
    print('Fs : {}\n'.format(len(Fs)))
    fig = plt.figure(count)
    plt.title('f(x) with tolerance:%.4f'%tol,{'fontsize':15})
    plt.xlabel('x',{'fontsize':15})
    ax = plt.subplot(111)
    ax.plot(Xs,Fs, linestyle='--', marker='o', color='b', label='Selected X')
    ax.legend(loc='upper left')
    count += 1
    plt.show()
# ans = adaptive_simpson( f, a, b, tol );
# np.linspace(start = 0, stop = 100, num = 5)