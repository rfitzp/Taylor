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

infile = open ("Taylor.out", "r")

tmax = 5.
fac  = 1.5

t1   = fac
t2   = fac/D**4
t3   = fac/Pphi**0.5/D

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

ti   = []
fir  = []
fii  = []
frir = []
frii = []
fscr = []
fsci = []

for line in infile: 

    numbers = line.split() 
    c1      = float(numbers[0])
    c2      = float(numbers[1])
    c3      = float(numbers[2])
    c4      = float(numbers[7])
    c5      = float(numbers[8])
    c6      = float(numbers[9])
    c7      = float(numbers[10])
    if c1 < tmax:
        ti  .append(c1)
        fir .append(c2)
        fii .append(c3)
        frir.append(c4)
        frii.append(c5)
        fscr.append(c6)
        fsci.append(c7)

fontsize = 15        

fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize = fontsize) 
plt.rc ('ytick', labelsize = fontsize)

plt.subplot (1, 1, 1)

plt.xlim (0., tmax)

plt.plot    (t,  fr,   color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"Re($\hat{\Psi}_0$)")
plt.plot    (t,  fi,   color = 'red',   linewidth = 2, linestyle = 'solid', label = r"Im($\hat{\Psi}_0$)")
plt.plot    (ti, fir,  color = 'blue',  linewidth = 2, linestyle = 'dotted')
plt.plot    (ti, fii,  color = 'red',   linewidth = 2, linestyle = 'dotted')
plt.plot    (ti, frir, color = 'blue',  linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, frii, color = 'red',   linewidth = 2, linestyle = 'dashed')
plt.plot    (ti, fscr, color = 'blue',  linewidth = 2, linestyle = 'dashdot')
plt.plot    (ti, fsci, color = 'red',   linewidth = 2, linestyle = 'dashdot')

plt.axhline (0.,       color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (t1,       color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (t2,       color = 'black', linewidth = 2, linestyle = 'dotted')
plt.axvline (t3,       color = 'black', linewidth = 2, linestyle = 'dotted')

plt.xlabel(r'$\hat{t}$', fontsize = fontsize)
plt.legend(fontsize = fontsize)

plt.tight_layout ();

#plt.show ()
plt.savefig("Figure8.pdf")
