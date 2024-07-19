import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Integrand.out", "r")

s  = []
p  = []
Wr = []
Wi = []
Vr = []
Vi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[1])
    c2      = float(numbers[2])
    c3      = float(numbers[3])
    c4      = float(numbers[4])
    c5      = float(numbers[5])
    p.append(c1)
    Wr.append(c2)
    Wi.append(c3)
    Vr.append(c4)
    Vi.append(c5)
                      
fig = plt.figure (figsize = (9.0, 8.0))
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (2, 2, 1)

plt.xlim (p[0], p[-1])

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,   color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\hat{\Delta}_s)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\hat{\Delta}_s)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (p[0], p[-1])

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,   color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Vr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(F_s)$")
plt.plot    (p, Vi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(F_s)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (-3., 3.)

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,   color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\hat{\Delta}_s)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\hat{\Delta}_s)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (-3., 3.)

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (0.,   color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Vr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(F_s)$")
plt.plot    (p, Vi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(F_s)$")

plt.xlabel (r'$\omega$', fontsize = "20")
plt.legend (fontsize = "15")

plt.tight_layout();

plt.show ()
