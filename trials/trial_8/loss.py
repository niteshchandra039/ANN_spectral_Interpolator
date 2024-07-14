from cmath import exp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit




epoch_loss = pd.read_csv('epoch_vs_loss.csv')
ep = epoch_loss['epoch'][1:]
ls = epoch_loss['loss'][1:]
x = np.log10(ep)
y = np.log10(ls)

x_ = x[10:1000]
y_ = y[10:1000]

#f = interpolate.UnivariateSpline(x_, y_)
xnew = np.log10(np.arange(10**1, 10**6, 1e03))

#ynew = f(xnew)   # use interpolationfunction returned by `interp1d`
plt.plot(x, y)

#plt.show()


def func(x, a, b, c, d): # simple quadratic example
    return a*np.exp(-b*(x-c)) + d

popt, pcov = curve_fit(func, x_, y_, bounds=(0, [1., 100., 0.5, 0.10]))

ynew = func(xnew, *popt)
plt.plot(xnew, ynew, 'r-',
         label='fit: a=%5.3f, b=%5.3f, c=%5.3f, d=%5.3f' % tuple(popt))

print(popt)
#popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
#popt

#plt.plot(xdata, func(xdata, *popt), 'g--',
#         label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))
#plt.xlabel('x')
#plt.ylabel('y')
plt.legend()
#plt.xlim(min(x_), max(x_))
#plt.xscale('log')
#plt.yscale('log')
plt.show()