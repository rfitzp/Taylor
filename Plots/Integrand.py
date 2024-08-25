import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Integrand.out", "r")

s  = []
p  = []
Wr = []
Wi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[1])
    c2      = float(numbers[5])
    c3      = float(numbers[6])
    p.append(c1)
    Wr.append(c2)
    Wi.append(c3)
                       
fig = plt.figure (figsize = (8.0, 6.0))
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (2, 1, 1)

plt.xlim (p[0], p[-1])

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\overline{\Psi}_0)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\overline{\Psi}_0)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (-3., 3.)

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\overline{\Psi}_0)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\overline{\Psi}_0)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.tight_layout ();

plt.show ()
