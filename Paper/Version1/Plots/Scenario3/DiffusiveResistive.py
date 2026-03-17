import json
import math
import numpy as np
import matplotlib.pyplot as plt

with open ("Namelist.json", "r") as f:
    data = json.load (f)

QE    = data["QE"]
Qe    = data["Qe"]
Qi    = data["Qi"]
D     = data["D"]
Pphi  = data["Pphi"]
Pperp = data["Pperp"]
Sigma = data["Sigma"]
iotae = Qe / (Qe - Qi)

fac = 1.5

t1 = fac/Pphi**0.5/D

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
    c2      = float(numbers[11])
    c3      = float(numbers[12])
    ti .append(c1)
    fir.append(c2)
    fii.append(c3)

fontsize = 15
    
fig = plt.figure (figsize=(8.0, 8.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (2, 1, 1)

plt.xlim (0., t[-1])

plt.plot    (t,  fr,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\hat{\Psi}_0$)")
plt.plot    (t,  fi,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\hat{\Psi}_0$)")
plt.plot    (ti, fir, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fii, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.axhline (0.,      color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (t1,      color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = fontsize)
plt.legend(fontsize = fontsize)

plt.subplot (2, 1, 2)

plt.xlim (0., 10.)

plt.plot    (t,  fr,  color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\hat{\Psi}_0$)")
plt.plot    (t,  fi,  color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\hat{\Psi}_0$)")
plt.plot    (ti, fir, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fii, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.axhline (0.,      color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (t1,      color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = fontsize)
plt.legend(fontsize = fontsize)

plt.tight_layout ();

#plt.show ()
plt.savefig("Figure9.pdf")
