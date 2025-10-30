import math
import numpy as np
import matplotlib.pyplot as plt

#with open("Analytic.py") as f:
#    code = f.read ()
#exec (code)

infile = open ("Taylor.out", "r")

t  = []
fr = []
fi = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    t .append(c1)
    fr.append(c2)
    fi.append(c3)

infile = open ("Analytic.out", "r")

ti  = []
fir = []
fii = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[3])
    c3      = float(numbers[4])
    ti .append(c1)
    fir.append(c2)
    fii.append(c3)
    
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (2, 1, 1)

plt.xlim (0., t[-1])

plt.plot    (t,  fr,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\Psi_0$)")
plt.plot    (t,  fi,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\Psi_0$)")
plt.plot    (ti, fir, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fii, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.axhline (0.,      color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.subplot (2, 1, 2)

plt.xlim (0., 20.)

plt.plot    (t,  fr,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\Psi_0$)")
plt.plot    (t,  fi,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\Psi_0$)")
plt.plot    (ti, fir, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fii, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.axhline (0.,      color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

plt.show ()
#plt.savefig("Taylor.pdf")
