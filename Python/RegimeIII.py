import matplotlib.pyplot as plt
import numpy as np

D = 0.9
xm = 2.5;
ym = 2.0

D2 = D*D
D4 = D2*D2
D6 = D4*D2

def f1(t):
    return 1/t**3.

def f2(t):
    return t**3

def f3(t):
    return 1./t**1.5

def f4(t):
    return 1./D2/t**2.

t1 = np.arange(1.e-3, 1., 1.e-3)
t2 = np.arange(1., xm, 1.e-3)
t3 = np.arange(1., 1./D4, 1.e-3)
t4 = np.arange(1./D4, xm, 1.e-3)

fig = plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

plt.vlines(x=1, ymin=0, ymax=1, linewidth=2, color='blue')
plt.vlines(x=1/D4, ymin=0, ymax=D6, linewidth=2, color='blue')
plt.hlines(y=D6, xmin=1./D4, xmax=xm, linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')
plt.plot(t4, f4(t4), linewidth=2, color='blue')

plt.text(1.7, 1.5, 'VR',  fontsize=25)
plt.text(0.95, 1.5, 'VI',  fontsize=25)
plt.text(0.5, 1., 'I',  fontsize=25)
plt.text(1.2, 0.35, 'RI',  fontsize=25)
plt.text(2.1, 0.35, 'DR',  fontsize=25)
plt.text(1.7, 0.1, 'SC',  fontsize=25)

plt.text(0.55, 1.5, '$P=\hat{t}^{-3}$', fontsize=20, color='red')
plt.text(1.17, 1.5, '$P=\hat{t}^{\,3}$', fontsize=20, color='red')
plt.text(1.2, 0.78, '$P=\hat{t}^{-3/2}$', fontsize=20, color='red')
plt.text(1.2, 0.15, '$\hat{t}=D^{-4}$', fontsize=20, color='red')
plt.text(1.9, 0.56, '$P=D^6$', fontsize=20, color='red')
plt.text(2.0, 0.09, '$P=\hat{t}^{-2} D^{-2}$', fontsize=20, color='red')
plt.text(0.78, 0.15, '$\hat{t}=1$', fontsize=20, color='red')

plt.plot([xm-0.02,xm], [1.,1.], linewidth=2, color='black')
plt.plot([1./D4,1./D4], [ym-0.02,ym], linewidth=2, color='black')

plt.text(xm+0.02, 0.955, '$P=1$', fontsize=20, color='black')

arr_width = .02
plt.arrow(0., 1.25, 2.0, 0.,  width = 0.0025, head_width = 3 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 0.74, 2.0, 0.,  width = 0.0025, head_width = 3 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 0.3, 2.3, 0.,  width = 0.0025, head_width = 3 * arr_width, head_length = 3 * arr_width, color='black') 

plt.text(2.07, 1.22, '1', fontsize=20, color='black')
plt.text(2.07, 0.71, '2', fontsize=20, color='black')
plt.text(2.37, 0.27, '3', fontsize=20, color='black')

plt.xlim(0., xm)
plt.ylim(0., ym)

plt.xlabel('$\hat{t}$', fontsize=18)
plt.ylabel('$P$', fontsize=18)

plt.savefig("RegimeIII.pdf")

#plt.show()
