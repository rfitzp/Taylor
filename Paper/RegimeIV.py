import matplotlib.pyplot as plt
import numpy as np

xm = 2.5
ym = 5.

D = 1.2

D2 = D*D
DM2 = 1./D2
D3 = D2*D
D4 = D2*D2
D6 = D4*D2

def f1(t):
    return 1./t**3.

def f2(t):
    return t**3.

def f3(t):
    return D4*t

def f4(t):
    return D2/t

def f5(t):
    return D2*t**2.

def f6(t):
    return D2/t**4.

def f7(t):
    return 1./D2/t**2.

t1 = np.arange(1.e-3, 1./D, 1.e-3)
t2 = np.arange(D2, xm, 1.e-3)
t3 = np.arange(1./D, D2, 1.e-3)
t4 = np.arange(1./D, 1., 1.e-3)
t5 = np.arange(1., D2, 1.e-3)
t6 = np.arange(1., D2, 1.e-3)
t7 = np.arange(D2, xm, 1.e-3)

fig = plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

plt.hlines(y=D6, xmin=D2, xmax=xm, linewidth=2, color='blue')
plt.vlines(x=D2, ymin=0, ymax=1./D6, linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')
plt.plot(t4, f4(t4), linewidth=2, color='blue')
plt.plot(t5, f5(t5), linewidth=2, color='blue')
plt.plot(t6, f6(t6), linewidth=2, color='blue')
plt.plot(t7, f7(t7), linewidth=2, color='blue')

plt.text(2.0, 3.5, 'VR',  fontsize=25)
plt.text(1.08, 3.5, 'VI',  fontsize=25)
plt.text(0.3, 3.5, 'I',  fontsize=25)
plt.text(1.065, 2., 'DI',  fontsize=20)
plt.text(2.0, 1.5, 'DR',  fontsize=25)
plt.text(1.5, 0.03, 'SC',  fontsize=20)

plt.text(0.67, 3.5, '$P=\hat{t}^{-3}$', fontsize=20, color='red')
plt.text(1.28, 3.5, '$P=\hat{t}^{\,3}$', fontsize=20, color='red')
plt.text(1.28, 2.19, '$P=D^2 \hat{t}^{\,2}$', fontsize=20, color='red')
plt.text(1.26, 0.07, '$\hat{t}=D^2$', fontsize=13, color='red')
plt.text(1.7, 3.02, '$P=D^6$', fontsize=20, color='red')
plt.text(2.16, 0.18, '$P=D^{-2}\,\hat{t}^{-2}$', fontsize=13, color='red')
plt.text(0.80, 2.35, '$P=D^4 \hat{t}$', fontsize=20, color='red')
plt.text(0.49, 1.36, '$P=D^2 \hat{t}^{-1}$', fontsize=20, color='red')
plt.text(1.17, 0.8, '$P=D^2 \hat{t}^{-4}$', fontsize=20, color='red')

plt.plot([xm-0.02,xm], [D2,D2], linewidth=2, color='black')
plt.plot([xm-0.02,xm], [D3,D3], linewidth=2, color='black')
plt.plot([xm-0.02,xm], [1./D6,1./D6], linewidth=2, color='black')
plt.plot([1/D,1/D], [0., 0.04], linewidth=2, color='black')
plt.plot([1.,1.], [0., 0.04], linewidth=2, color='black')
#plt.plot([D,D], [3.93, 4.01], linewidth=2, color='black')

plt.text(xm+0.01, D2-0.08, '$P=D^2$', fontsize=20, color='black')
plt.text(xm+0.01, D3-0.08, '$P=D^3$', fontsize=20, color='black')
plt.text(xm+0.01, 1./D6-0.08, '$P=D^{-6}$', fontsize=17, color='black')
#plt.text(0.62, 4.02, '$R=D^{-2}$', fontsize=20, color='black')
#plt.text(0.93, 4.02, '$R=1$', fontsize=20, color='black')
plt.text(1/D-0.13, 0-0.23, '$\hat{t}=D^{-1}$', fontsize=15, color='black')

arr_width = .02
plt.arrow(0., 4.5,  2.0, 0.,  width = 0.0025, head_width = 6 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 2.7,  2.0, 0.,  width = 0.0025, head_width = 6 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 1.7,  1.6, 0.,  width = 0.0025, head_width = 6 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 1.2,  2.0, 0.,  width = 0.0025, head_width = 6 * arr_width, head_length = 3 * arr_width, color='black')
plt.arrow(0., 0.28, 2.0, 0.,  width = 0.0025, head_width = 6 * arr_width, head_length = 3 * arr_width, color='black')

plt.text(2.07, 4.42, '1', fontsize=20, color='black')
plt.text(2.07, 2.62, '4', fontsize=20, color='black')
plt.text(1.67, 1.62, '5', fontsize=20, color='black')
plt.text(2.07, 1.12, '6', fontsize=20, color='black')
plt.text(2.07, 0.2,  '7', fontsize=20, color='black')

plt.xlim(0., xm)
plt.ylim(0., ym)

plt.xlabel('$\hat{t}$', fontsize=18)
plt.ylabel('$P$', fontsize=18)

plt.savefig("RegimeIV.pdf")

#plt.show()
