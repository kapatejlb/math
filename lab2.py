from sympy import *
from matplotlib import mlab
import matplotlib.pyplot as plt
import math

def find_difference(ylist1, ylist2, meta):
    difference = []
    ylist2 = ylist2[::2]
    for i in range(len(ylist1)):
        difference.append(math.fabs(ylist1[i] - ylist2[i]))

    # print(str(max(difference))+"-----------------------------------------------------------------=------------------------------------------------------------------")
    # difference = max(difference)
    return max(difference) > meta

def find_dots(x1, x2, n, x, q, f):
    h = (x2 - x1) / (n)
    print("h = " + str(h))

    xlist = mlab.frange(x1, x2, h)

    print("X's:")
    print(len(xlist))
    print(xlist)

    # y" = -(1 + x^2)y - 1
    # ay" + (1 + bx^2)y = -1

    yks = []
    for i in range(n+1):
        yks.append(Symbol("y_d" + str(i + 1)))

    q = lambdify(x, q, 'numpy')
    f = lambdify(x, f, 'numpy')
    system = []

    # for i in range(n - 2):
    #     system.append((yks[i + 2] - 2 * yks[i + 1] + yks[i])/(h**2) - q(xlist[i + 1]) * yks[i + 1] - f(xlist[i + 1]))
    #     #system.append(yks[i] - (2 + h ** 2 * q(xlist[i + 1])) * yks[i + 1] + yks[i + 2] - h ** 2 * f(xlist[i + 1]))
                                            #^^^^^деление на нолb
    k = 1
    while k!=n:
        system.append(yks[k-1] - (2 - h**2 * q(xlist[k]))*yks[k] + yks[k+1] - f(xlist[k])*h**2)
        k += 1

    print("\nGain system:")
    for el in system:
        print(el)

    system[0] = system[0].subs(yks[0], 0)
    system[-1] = system[-1].subs(yks[-1], 0)
    print("\nThen we substitute the boundary points and get solvable system :")
    for el in system:
        print(el)


    coeffs = linsolve(system, yks[1 : n])
    coeffs = next(iter(coeffs))
    print("\nAnswer:")
    print(coeffs)

    ylist = [0]
    for el in coeffs:
        ylist.append(el)
    ylist.append(0)

    print("\nYk's :")
    print(ylist)

    return xlist, ylist

def func(meta, x1, x2, n, x, q, f):
    ylist = [1]
    ylist2 = [0]
    while find_difference(ylist, ylist2, meta):
        xlist, ylist = find_dots(x1, x2, n, x, q, f)
        n *= 2
        xlist2, ylist2 = find_dots(x1, x2, n, x, q, f)
        # print(ylist)
        # print(ylist2)

        plt.plot(xlist, ylist)
        plt.plot(xlist2, ylist2)
        plt.show()

meta = 10**-3
x1 = -1
x2 = 1
n = 3
x = Symbol("x")
a = 1
b = 1
p = 0/a
q = (1 + b * x**2)/a
f = -1/a

yks = []
for i in range(n + 1):
    yks.append(Symbol("y_d" + str(i + 1)))

qwe = yks[-2] + p*yks[-3] + q*yks[-4] - f
print(qwe)

func(meta, x1, x2, n, x, q, f)

print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
meta = 10**-3
x1 = -1
x2 = 1
n = 3
x = Symbol("x")
a = sin(13)
b = cos(13)
p = 0/a
q = (1 + b * x**2)/a
f = -1/a
func(meta, x1, x2, n, x, q, f)

# print("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n")
# meta = 10**-1
# x1 = -1
# x2 = 1
# n = 3
# x = Symbol("x")
# a = sin(13*x)
# b = cos(13*x)
# p = 0/a
# q = (1 + b * x**2)/a
# f = -1/a
# func(meta, x1, x2, n, x, q, f)