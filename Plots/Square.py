import math
import numpy as np
import matplotlib.pyplot as plt

x   = input("tau: ")
tau = float(x)

infile = open ("Taylor.out", "r")

t  = []
fr = []
fi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    t.append(c1)
    fr.append(c2)
    fi.append(c3)

dt = t[1] - t[0]
n  = int (tau/dt)

Fr = []
Fi = []
for i in range (len(t)):
    if i < n:
        Fr.append (fr[i])
        Fi.append (fi[i])
    else:
        Fr.append (fr[i]-fr[i-n])
        Fi.append (fi[i]-fi[i-n])
    
print ("dt = %10.3e tau = %10.3e n = %4d" % (dt, tau, n))    
                      
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (1, 1, 1)

plt.xlim (t[0]-0.02*tau, t[-1])

mx = 1.05*max(Fr)

x = [0., 0., tau, tau]
y = [0., mx, mx,  0.]

plt.plot    (t, Fr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\hat{\Psi}_0)$")
plt.plot    (t, Fi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\hat{\Psi}_0)$")
plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (x, y,  color = 'black', linewidth = 3, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

plt.show ()
#plt.savefig("Figure3.pdf")
