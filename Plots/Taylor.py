import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Taylor.out", "r")

t  = []
f  = []
fe = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    t.append(c1)
    f.append(c2)
    fe.append(c3)
                      
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (2, 1, 1)

plt.xlim (t[0], t[-1])

plt.plot    (t, f,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\Psi_0$)")
plt.plot    (t, fe, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\Psi_0$)")
plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 10.)

plt.plot    (t, f,  color = 'blue',  linewidth = 2, linestyle = 'solid', label=r"Re($\Psi_0$)")
plt.plot    (t, fe, color = 'red',   linewidth = 2, linestyle = 'solid', label=r"Im($\Psi_0$)")
plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

plt.show ()
