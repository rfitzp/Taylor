import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

QE  = 0.1
Qe  = 1.0
Pp  = 5.0
Sig = 10.
eVR = math.gamma(1./6.) /6.**(2./3.) /math.pi /math.gamma(5./6.) /Pp**(1./6.) /Sig
tau = 20.

def Ftaylor (t):

    if (t < 0.):
        return 0.
    else:
        fac1 = eVR /(1j * (QE + Qe) + eVR)
        fac2 = 1. - cmath.exp (- (1j * (QE + Qe) + eVR) * t)
        
        return fac1 * fac2

tt = np.arange (0., 100.01, 0.01)

Psir = []
Psii = []

for t in tt:

    Psi = Ftaylor(t) - Ftaylor(t-tau)

    Psir.append (Psi.real)
    Psii.append (Psi.imag)
        
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (1, 1, 1)

plt.xlim (tt[0]-0.02*tau, tt[-1])

mx = 1.05*max(Psir)

x = [0., 0., tau, tau]
y = [0., mx, mx,  0.]

plt.plot    (tt, Psir, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\hat{\Psi}_0)$")
plt.plot    (tt, Psii, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\hat{\Psi}_0)$")
plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (x, y,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

#plt.show ()
plt.savefig("Figure11.pdf")
