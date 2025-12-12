import matplotlib.pyplot as plt
import numpy as np

D = 0.9

D2 = D*D
D4 = D2*D2
D6 = D4*D2

def f1(q):
    return q*q/D2

def f2(q):
    return q**1.5

def f3(q):
    return q**3.

def f4(q):
    return 1./q**3.

t1 = np.arange(0., D4, 1.e-3)
t2 = np.arange(D4, 1., 1.e-3)
t3 = np.arange(1., 2.**(1./3.), 1.e-3)
t4 = np.arange(2.**(-1./3.), 1., 1.e-3)

fig = plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

plt.hlines(y=D6, xmin=0, xmax=D4, linewidth=2, color='blue')
plt.vlines(x=D4, ymin=0, ymax=D6, linewidth=2, color='blue')
plt.vlines(x=1., ymin=0, ymax=1., linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')
plt.plot(t4, f4(t4), linewidth=2, color='blue')

plt.text(0.4, 1.4, 'VR',  fontsize=25)
plt.text(0.975, 1.4, 'VI',  fontsize=25)
plt.text(1.2, 0.5, 'I',  fontsize=25)
plt.text(0.8, 0.5, 'RI',  fontsize=25)
plt.text(0.2, 0.3, 'DR',  fontsize=25)
plt.text(0.5, 0.1, 'SC',  fontsize=25)

plt.text(1.17, 1.5, '$P=R^3$', fontsize=20, color='red')
plt.text(0.66, 1.5, '$P=R^{-3}$', fontsize=20, color='red')
plt.text(0.64, 0.75, '$P=R^{3/2}$', fontsize=20, color='red')
plt.text(0.665, 0.15, '$R=D^4$', fontsize=20, color='red')
plt.text(0.25, 0.56, '$P=D^6$', fontsize=20, color='red')
plt.text(0.01, 0.07, '$P=R^2 D^{-2}$', fontsize=20, color='red')
plt.text(1.01, 0.15, '$R=1$', fontsize=20, color='red')

plt.plot([1.38,1.41], [1.,1.], linewidth=2, color='black')

plt.text(1.41, 0.95, '$P=1$', fontsize=20, color='black')

plt.xlim(0., 1.4)
plt.ylim(0., 2.)

plt.xlabel('$R$', fontsize=18)
plt.ylabel('$P$', fontsize=18)

plt.savefig("RegimeI.pdf")

#plt.show()
