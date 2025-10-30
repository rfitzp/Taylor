import json
import numpy as np
import math
import cmath
from scipy.integrate import quad
from scipy.special import j0

pi   = math.pi
pisq = math.sqrt (pi)

# ####################################
# Read calculation data from JSON file
# ####################################

with open ("../Inputs/Namelist.json", "r") as f:
    data = json.load (f)

QE    = data["QE"]
Qe    = data["Qe"]
Qi    = data["Qi"]
D     = data["D"]
Pphi  = data["Pphi"]
Pperp = data["Pperp"]
Sigma = data["Sigma"]
iotae = Qe / (Qe - Qi)

tmax = data["tmax"]
Nt   = data["Nt"]

print ("\nQE   = %10.3e Qe = %10.3e Qi = %10.3e D = %10.3e Pphi = %10.3e Pperp = %10.3e Sigma = %10.3e"
       % (QE, Qe, Qi, D, Pphi, Pperp, Sigma))
print ("tmax = %10.3e Nt =  %4d" % (tmax, Nt))

# ########################
# Calculate epsilon values
# ########################

#Sigma *= 0.75

epsI  = 2. /pi /Sigma
epsVI = 1. /pisq /Pphi**0.25 /Sigma
epsDI = 2. * D /pi /(iotae * Pphi)**0.5 /Sigma
epsVR = math.gamma (1./6.) /6.**(2./3.) /pi /math.gamma (5./6.) /Pphi**(1./6.) /Sigma
epsDR = math.gamma (0.25) * D**0.5 /2./pi /math.gamma(0.75) / (iotae*Pperp)**0.25 /Sigma

print ("\nepsI = %10.3e epsVI = %10.3e epsDI = %10.3e epsVR = %10.3e epsDR = %10.3e"
       % (epsI, epsVI, epsDI, epsVR, epsDR))

# #########################
# Define analytic integrads
# #########################
def FunIr (tp, t):

    if QE + Qe == 0:
        fac1 = epsI
        fac2 = t - tp
    else:
        fac1 = epsI /( (QE + Qe) * 1j )
        fac2 = 1. - cmath.exp ( (QE + Qe) * 1j * (tp - t) )
    fac3 = cmath.exp (- (QE + 0.5*Qi) * 1j * tp )
    fac4 = j0 (0.5*Qi * tp)

    return (fac1 * fac2 * fac3 * fac4).real

def FunIi (tp, t):
    
    if QE + Qe == 0:
        fac1 = epsI
        fac2 = t - tp
    else:
        fac1 = epsI /( (QE + Qe) * 1j )
        fac2 = 1. - cmath.exp ( (QE + Qe) * 1j * (tp - t) )
    fac3 = cmath.exp (- (QE + 0.5*Qi) * 1j * tp )
    fac4 = j0 (0.5*Qi * tp)

    return (fac1 * fac2 * fac3 * fac4).imag

def FunVIr (tp, t):

    fac1 = epsVI /math.gamma (1.25)
    fac2 = tp**0.25 * cmath.exp (- (QE + Qe) * 1j * tp)

    return (fac1 * fac2).real

def FunVIi (tp, t):

    fac1 = epsVI /math.gamma (1.25)
    fac2 = tp**0.25 * cmath.exp (- (QE + Qe) * 1j * tp)

    return (fac1 * fac2).imag

def FunDRr (t):

    fac1 = epsDR / ( (QE + Qe) * 1j + epsDR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsDR ) * t)

    return (fac1 * fac2).real

def FunDRi (t):

    fac1 = epsDR / ( (QE + Qe) * 1j + epsDR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsDR ) * t)

    return (fac1 * fac2).imag

def FunVRr (t):

    fac1 = epsVR / ( (QE + Qe) * 1j + epsVR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsVR ) * t)

    return (fac1 * fac2).real

def FunVRi (t):

    fac1 = epsVR / ( (QE + Qe) * 1j + epsVR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsVR ) * t)

    return (fac1 * fac2).imag

# ############################
# Calculate analytic solutions
# ############################

tt = np.linspace (0., tmax, Nt)

with open ("Analytic.out", "w") as f:

    for t in tt:
    
        Ir, errr = quad (FunIr,  0., t, args = (t))
        Ii, erri = quad (FunIi,  0., t, args = (t))
        Jr, errr = quad (FunVIr, 0., t, args = (t))
        Ji, erri = quad (FunVIi, 0., t, args = (t))

        f.write ("%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n"
                 % (t, Ir, Ii, FunDRr (t), FunDRi (t), FunVRr (t), FunVRi (t), Jr, Ji))

f.close ()    


        
