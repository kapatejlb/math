from sympy import *
from matplotlib import mlab
import pylab
import mpmath
import random

import matplotlib.pyplot as plt

def finding_lincomb(quality, c, x, p, q, f):
    print("Linear combination type:")
    linear_combination = 0
    for i in range(quality):
        linear_combination += c[i] * (x ** (i)) * (1 - x**2)
    print(linear_combination)

    print("\nFirst derivative:")
    y_d1 = linear_combination.diff(x)
    print(y_d1)

    print("\nSecond derivative:")
    y_d2 = y_d1.diff(x)
    print(y_d2)

    #y"+py'+qy=f(x)
    print("\nAfter substitution we get:")
    y = y_d2 + p * y_d1 + q * linear_combination - f
    print(y)
    return y, linear_combination

def plot_graphic(linear_combination, x, xmin, xmax, dx):
    f = lambdify(x, linear_combination, 'numpy')
    # xmin = -10000
    # xmax = 10000
    # dx = 0.1
    xlist = mlab.frange(xmin, xmax, dx)
    ylist = [f(x) for x in xlist]
    # pylab.plot(xlist, ylist)
    #
    # pylab.axhline(0, color='black')
    # pylab.axvline(0, color='black')
    # pylab.show()
    return xlist, ylist

def colloc_dots(quality):
    collocation_dots = []
    for i in range(quality):
        collocation_dots.append(random.randrange(-999, 999)/1000)
    return collocation_dots

def Galerkin_method(quality, x, p, q, f, c):

    y, linear_combination = finding_lincomb(quality, c, x, p, q, f)

    print("\nOur f's are:")
    fi = []
    for i in range(quality):
        fi.append((x ** (i)) * (1 - x**2))
    print(fi)

    print("\nAfter that gain system:")
    system = []
    for i in range(quality):
        system.append(integrate(fi[i]*y, (x, -1, 1)))
    print(system)

    print("\nAnswer:")
    coeffs = linsolve(system, c)
    coeffs = next(iter(coeffs))

    for i in range(quality):
        linear_combination = linear_combination.copy().subs({c[i] : coeffs[i]})

    print(linear_combination)

    xlist, ylist = plot_graphic(linear_combination, x, -10000, 10000, 0.1)
    return xlist, ylist
def discrete_MLS(quality, x, p, q, f, c):
    N = int(input("\nEnter N- num of dots:"))
    colloc_dos = colloc_dots(N)

    y, linear_combination = finding_lincomb(quality, c, x, p, q, f)

    print("\nNow we get S:")
    s = 0
    for i in range(N):
        s += lambdify(x, y.copy()**2, 'numpy')(colloc_dos[i])
    print(s)

    print("\nSystem after differentiation:")
    system = []
    for i in range(quality):
        system.append(s.copy().diff(c[i]))
    print(system)

    print("\nAnswer:")
    coeffs = linsolve(system, c)
    coeffs = next(iter(coeffs))

    for i in range(quality):
        linear_combination = linear_combination.copy().subs({c[i] : coeffs[i]})

    print(linear_combination)

    xlist, ylist = plot_graphic(linear_combination, x, -10000, 10000, 0.1)
    return xlist, ylist
def collocation(quality, x, p, q, f, c):
    collocation_dots = colloc_dots(quality)

    y, linear_combination = finding_lincomb(quality, c, x, p, q, f)

    print("\nAnd then gain system:")
    system = []
    for i in range(quality):
        system.append(y.copy().subs(x, collocation_dots[i]))

    print(system)

    print("\nAnswer:")
    coeffs = linsolve(system, c)
    coeffs = next(iter(coeffs))

    for i in range(quality):
        linear_combination = linear_combination.copy().subs({c[i] : coeffs[i]})

    print(linear_combination)

    xlist, ylist = plot_graphic(linear_combination, x, -10000, 10000, 0.1)
    return xlist, ylist
def least_square_method(quality, x, p, q, f, c):

    y, linear_combination = finding_lincomb(quality, c, x, p, q, f)

    print("\nIntegral:")
    integral = integrate(y**2, (x, -1, 1))
    print(integral)

    print("\nSystem after differentiation:")
    system = []
    for i in range(quality):
        system.append(integral.copy().diff(c[i]))

    print(system)

    print("\nAnswer:")
    coeffs = linsolve(system, c)
    coeffs = next(iter(coeffs))

    for i in range(quality):
        linear_combination = linear_combination.copy().subs({c[i] : coeffs[i]})

    print(linear_combination)

    xlist, ylist = plot_graphic(linear_combination, x, -10000, 10000, 0.1)
    return xlist, ylist

quality = 3
x = Symbol("x")
p = 0.0
q = 1.0 + x**2
f = -1.0
c = []
for i in range(quality):
    c.append(Symbol("c"+str(i+1)))

xlist, ylist = collocation(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'b')
xlist, ylist = least_square_method(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'r')
xlist, ylist = discrete_MLS(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'c')
xlist, ylist = Galerkin_method(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'y')

plt.legend(("Collocation", "Least Square Method", "Discrete MLS", "Galerkin Method"))
pylab.axhline(0, color='black')
pylab.axvline(0, color='black')
plt.show()


print("---------------------------------------------------------------")
quality = 5
x = Symbol("x")
a = sin(mpmath.radians(13))
b = cos(mpmath.radians(13))
p = 0/a
q = (1 + b * x**2)/a
f = -1/a
c = []
for i in range(quality):
    c.append(Symbol("c"+str(i+1)))

xlist, ylist = collocation(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'b')
xlist, ylist = least_square_method(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'r')
xlist, ylist = discrete_MLS(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'c')
xlist, ylist = Galerkin_method(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'y')

plt.legend(("Collocation", "Least Square Method", "Discrete MLS", "Galerkin Method"))
pylab.axhline(0, color='black')
pylab.axvline(0, color='black')
plt.show()


print("---------------------------------------------------------------")
quality = 5
x = Symbol("x")
a = sin(mpmath.radians(13*x))
b = cos(mpmath.radians(13*x))
p = 0/a
q = (1 + b * x**2)/a
f = -1/a
c = []
for i in range(quality):
    c.append(Symbol("c"+str(i+1)))

xlist, ylist = collocation(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'b')
# xlist, ylist = least_square_method(quality, x, p, q, f, c)
# plt.plot(xlist, ylist, 'r')
xlist, ylist = discrete_MLS(quality, x, p, q, f, c)
plt.plot(xlist, ylist, 'c')
# xlist, ylist = Galerkin_method(quality, x, p, q, f, c)
# plt.plot(xlist, ylist, 'y')

plt.legend(("Collocation", "Discrete MLS"))
pylab.axhline(0, color='black')
pylab.axvline(0, color='black')
plt.show()
