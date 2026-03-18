import math
import numpy as np
import matplotlib.pyplot as plt

N   = 11
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
n  = int (tau/dt/N)

ff = []
for i in range (N):
    ff.append(math.sin(i*math.pi/(N-1)) + 0.05)

tt = []
f1 = []
for i in range (N):
    if i == 0:
        tt.append (i*tau/(N-1))
        f1.append (math.sin(i*math.pi/(N-1)) + 0.05)
        tt.append (i*tau/(N-1))
        f1.append (math.sin((i+1)*math.pi/(N-1)) + 0.05)
    elif i == N-1:
        pass
        #tt.append (i*tau/(N-1))
        #ff.append (math.sin((i+1)*math.pi/(N-1)) + 0.05)
    else:
        tt.append (i*tau/(N-1))
        f1.append (math.sin(i*math.pi/(N-1)) + 0.05)
        tt.append (i*tau/(N-1))
        f1.append (math.sin((i+1)*math.pi/(N-1)) + 0.05)

Fr = []
Fi = []
for i in range (len(t)):
    if i < n:       
        Fr.append (ff[0]*fr[i])
        Fi.append (ff[0]*fi[i])
    elif i >= n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n])) 
    elif i > 2*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]))
    elif i > 3*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]))
    elif i > 4*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n]))
    elif i > 5*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]))
    elif i > 6*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]) + ff[5]*(fr[i-5*n]-fr[i-6*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]) + ff[5]*(fi[i-5*n]-fi[i-6*n]))
    elif i > 7*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]) + ff[5]*(fr[i-5*n]-fr[i-6*n]) + ff[6]*(fr[i-6*n]-fr[i-7*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]) + ff[5]*(fi[i-5*n]-fi[i-6*n]) + ff[6]*(fi[i-6*n]-fi[i-7*n]))
    elif i > 8*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]) + ff[5]*(fr[i-5*n]-fr[i-6*n]) + ff[6]*(fr[i-6*n]-fr[i-7*n])
                   + ff[7]*(fr[i-7*n]-fr[i-8*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]) + ff[5]*(fi[i-5*n]-fi[i-6*n]) + ff[6]*(fi[i-6*n]-fi[i-7*n])
                   + ff[7]*(fi[i-7*n]-fi[i-8*n]))
    elif i > 9*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]) + ff[5]*(fr[i-5*n]-fr[i-6*n]) + ff[6]*(fr[i-6*n]-fr[i-7*n])
                   + ff[7]*(fr[i-7*n]-fr[i-8*n]) + ff[8]*(fr[i-8*n]-fr[i-9*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]) + ff[5]*(fi[i-5*n]-fi[i-6*n]) + ff[6]*(fi[i-6*n]-fi[i-7*n])
                   + ff[7]*(fi[i-7*n]-fi[i-8*n]) + ff[8]*(fi[i-8*n]-fi[i-9*n]))
    elif i > 10*n:
        Fr.append (ff[0]*(fr[i]-fr[i-n]) + ff[1]*(fr[i-n]-fr[i-2*n]) + ff[2]*(fr[i-2*n]-fr[i-3*n]) + ff[3]*(fr[i-3*n]-fr[i-4*n])
                   + ff[4]*(fr[i-4*n]-fr[i-5*n]) + ff[5]*(fr[i-5*n]-fr[i-6*n]) + ff[6]*(fr[i-6*n]-fr[i-7*n])
                   + ff[7]*(fr[i-7*n]-fr[i-8*n]) + ff[8]*(fr[i-8*n]-fr[i-9*n]) + ff[9]*(fr[i-9*n]-fr[i-10*n]))
        Fi.append (ff[0]*(fi[i]-fi[i-n]) + ff[1]*(fi[i-n]-fi[i-2*n]) + ff[2]*(fi[i-2*n]-fi[i-3*n]) + ff[3]*(fi[i-3*n]-fi[i-4*n])
                   + ff[4]*(fi[i-4*n]-fi[i-5*n]) + ff[5]*(fi[i-5*n]-fi[i-6*n]) + ff[6]*(fi[i-6*n]-fi[i-7*n])
                   + ff[7]*(fi[i-7*n]-fi[i-8*n]) + ff[8]*(fi[i-8*n]-fi[i-9*n]) + ff[9]*(fi[i-9*n]-fi[i-10*n]))
        
print ("dt = %10.3e tau = %10.3e n = %4d" % (dt, tau, n))
print (ff)
                      
fig = plt.figure (figsize=(8.0, 6.0))
plt.rc ('xtick', labelsize=17) 
plt.rc ('ytick', labelsize=17)

plt.subplot (1, 1, 1)

plt.xlim (t[0]-0.02*tau, t[-1])

tt.insert(0,  tt[0])
tt.append(tt[-1])
f1.insert(0,  0.)
f1.append(0.)

mx = 1.05*max(Fr)
y  = mx*np.asarray(f1)

plt.plot    (t, Fr, color = 'blue',  linewidth = 2, linestyle = 'solid', label = r"$Re(\hat{\Psi}_0)$")
plt.plot    (t, Fi, color = 'red',   linewidth = 2, linestyle = 'solid', label = r"$Im(\hat{\Psi}_0)$")
plt.axhline (0.,    color = 'black', linewidth = 2, linestyle = 'dotted')
plt.plot    (tt, y, color = 'black', linewidth = 3, linestyle = 'solid')

plt.xlabel(r'$\hat{t}$', fontsize = "20")
plt.legend(fontsize = "15")

plt.tight_layout ();

plt.show ()
#plt.savefig("Figure4.pdf")
