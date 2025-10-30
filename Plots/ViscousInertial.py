import math
import numpy as np
import matplotlib.pyplot as plt

#with open("Analytic.py") as f:
#    code = f.read ()
#exec (code)

infile = open ("Taylor.out", "r")

x = input("tmax: ")

tmax = float(x)

t  = []
fr = []
fi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    if c1 < tmax:
        t .append(c1)
        fr.append(c2)
        fi.append(c3)

infile = open ("Analytic.out", "r")

ti  = []
fir = []
fii = []
fvr = []
fvi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[7])
    c5      = float(numbers[8])
    if c1 < tmax:
        ti .append(c1)
        fir.append(c4)
        fii.append(c5)
        fvr.append(c2)
        fvi.append(c3)
    
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (1, 1, 1)

plt.xlim (0., tmax)

plt.plot    (t,  fr,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\Psi_0$)")
plt.plot    (t,  fi,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\Psi_0$)")
plt.plot    (ti, fvr, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fvi, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fir, color = 'blue',  linewidth = 2, linestyle = 'dotted')
plt.plot    (ti, fii, color = 'red',   linewidth = 2, linestyle = 'dotted')
plt.axhline (0.,      color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

plt.show ()
#plt.savefig("Taylor.pdf")
