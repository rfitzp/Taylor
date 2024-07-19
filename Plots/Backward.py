import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Backward.out", "r")

p  = []
Wr = []
Wi = []
Vr = []
Vi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = float(numbers[4])
    p.append(c1)
    Wr.append(c2)
    Wi.append(c3)
    Vr.append(c4)
    Vi.append(c5)
                      
fig = plt.figure (figsize = (9.0, 8.0))
plt.rc ('xtick', labelsize = 17) 
plt.rc ('ytick', labelsize = 17)

plt.subplot (2, 2, 1)

plt.xlim (0., p[0])

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(W)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(W)$")

plt.xlabel (r'$p$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 2)

plt.xlim (0., p[0])

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Vr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(V)$")
plt.plot    (p, Vi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(V)$")

plt.xlabel (r'$p$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 3)

plt.xlim (0., 1.)
plt.ylim (-2., 2.)

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Wr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(W)$")
plt.plot    (p, Wi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(W)$")

plt.xlabel (r'$p$', fontsize = "20")
plt.legend (fontsize = "15")

plt.subplot (2, 2, 4)

plt.xlim (0., 1.)
plt.ylim (-2., 2.)

plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (p, Vr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(V)$")
plt.plot    (p, Vi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(V)$")

plt.xlabel (r'$p$', fontsize = "20")
plt.legend (fontsize = "15")

plt.tight_layout();

plt.show ()
