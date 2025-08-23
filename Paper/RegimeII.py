import matplotlib.pyplot as plt
import numpy as np

D = 1.2

D2 = D*D
DM2 = 1./D2
D3 = D2*D
D4 = D2*D2
D6 = D4*D2

def f1(q):
    return q*q/D2

def f2(q):
    return D2/q/q

def f3(q):
    return q**3.

def f4(q):
    return 1./q**3.

def f5(q):
    return D4/q

def f6(q):
    return D2*q

def f7(q):
    return D2*q**4.

t1 = np.arange(0., DM2, 1.e-3)
t2 = np.arange(DM2, 1, 1.e-3)
t3 = np.arange(D, 4.**(1./3.), 1.e-3)
t4 = np.arange(4.**(-1./3.), DM2, 1.e-3)
t5 = np.arange(DM2, D, 1.e-3)
t6 = np.arange(1., D, 1.e-3)
t7 = np.arange(DM2, 1., 1.e-3)

fig = plt.figure(figsize=(10,7))
plt.rc('xtick', labelsize=14)
plt.rc('ytick', labelsize=14)

plt.hlines(y=D6, xmin=0, xmax=DM2, linewidth=2, color='blue')
plt.vlines(x=DM2, ymin=0, ymax=1./D6, linewidth=2, color='blue')
plt.plot(t1, f1(t1), linewidth=2, color='blue')
plt.plot(t2, f2(t2), linewidth=2, color='blue')
plt.plot(t3, f3(t3), linewidth=2, color='blue')
plt.plot(t4, f4(t4), linewidth=2, color='blue')
plt.plot(t5, f5(t5), linewidth=2, color='blue')
plt.plot(t6, f6(t6), linewidth=2, color='blue')
plt.plot(t7, f7(t7), linewidth=2, color='blue')

plt.text(0.3, 3.5, 'VR',  fontsize=25)
plt.text(1.1, 3.5, 'VI',  fontsize=25)
plt.text(1.4, 1.0, 'I',  fontsize=25)
plt.text(0.97, 1.7, 'DI',  fontsize=25)
plt.text(0.4, 1.0, 'DR',  fontsize=25)
plt.text(0.55, 0.03, 'SC',  fontsize=20)

plt.text(1.34, 2.25, '$P=R^3$', fontsize=20, color='red')
plt.text(0.665, 3.5, '$P=R^{-3}$', fontsize=20, color='red')
plt.text(0.63, 1.5, '$P=D^2 R^{-2}$', fontsize=20, color='red')
plt.text(0.7, 0.1, '$R=D^{-2}$', fontsize=20, color='red')
plt.text(0.25, 3.04, '$P=D^6$', fontsize=20, color='red')
plt.text(0.35, 0.28, '$P=R^2 D^{-2}$', fontsize=20, color='red')
plt.text(0.91, 2.3, '$P=D^4 R^{-1}$', fontsize=20, color='red')
plt.text(1.08, 1.4, '$P=D^2 R$', fontsize=20, color='red')
plt.text(0.9, 0.8, '$P=D^2 R^4$', fontsize=20, color='red')

plt.plot([1.58,1.61], [D2,D2], linewidth=2, color='black')
plt.plot([1.58,1.61], [D3,D3], linewidth=2, color='black')
plt.plot([1.58,1.61], [1./D6,1./D6], linewidth=2, color='black')
plt.plot([0.694,0.694], [3.93, 4.01], linewidth=2, color='black')
plt.plot([1.,1.], [3.93, 4.01], linewidth=2, color='black')
plt.plot([D,D], [3.93, 4.01], linewidth=2, color='black')

plt.text(1.61, 1.35, '$P=D^2$', fontsize=20, color='black')
plt.text(1.61, 1.668, '$P=D^3$', fontsize=20, color='black')
plt.text(1.61, 0.25, '$P=D^{-6}$', fontsize=17, color='black')
plt.text(0.62, 4.02, '$R=D^{-2}$', fontsize=20, color='black')
plt.text(0.93, 4.02, '$R=1$', fontsize=20, color='black')
plt.text(1.13, 4.02, '$R=D$', fontsize=20, color='black')

plt.xlim(0., 1.6)
plt.ylim(0., 4.)

plt.xlabel('$R$', fontsize=18)
plt.ylabel('$P$', fontsize=18)

plt.savefig("RegimeII.pdf")

#plt.show()
