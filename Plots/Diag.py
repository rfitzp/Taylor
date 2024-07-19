import math
import numpy as np
import matplotlib.pyplot as plt

infile = open ("Diagnostic.out", "r")

t    = []
hmax = []
hmin = []
err  = []
rept = []
stp  = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[3])
    c5      = int(numbers[4])
    c6      = float(numbers[5])
    t.append(c1)
    hmin.append(c2)
    hmax.append(c3)
    stp.append(c6)

fig = plt.figure (figsize=(12.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot(1, 2, 1)

plt.xlim (t[0], t[-1])

plt.plot (t, hmin, color='blue', linewidth = 2, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize="20")
plt.ylabel(r'$\log_{10}(h)$', fontsize="20")

plt.subplot(1, 2, 2)

plt.xlim (t[0], t[-1])
#plt.ylim (0., 5.)

plt.plot (t, stp, color='red', linewidth = 2, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize="20")
plt.ylabel(r'$\log_{10}$(steps)', fontsize="20")

plt.show ()
