import json
import numpy as np
import math
import cmath
from scipy.integrate import quad
from scipy.special import j0
from scipy.special import jv

pi   = math.pi
pisq = math.sqrt (pi)

# ####################################
# Read calculation data from JSON file
# ####################################

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

tmax = data["tmax"]
Nt   = data["Nt"]

print ("\nQE   = %10.3e Qe = %10.3e Qi = %10.3e D = %10.3e Pphi = %10.3e Pperp = %10.3e Sigma = %10.3e"
       % (QE, Qe, Qi, D, Pphi, Pperp, Sigma))
print ("tmax = %10.3e Nt =  %4d" % (tmax, Nt))

# ########################
# Calculate epsilon values
# ########################

epsI  = 2. /pi /Sigma
epsVI = 1. /pisq /Pphi**0.25 /Sigma
epsDI = 2. * D /pi /(iotae * Pphi)**0.5 /Sigma
epsRI = math.gamma(0.25) / (2. * pi * math.gamma(0.75) * Sigma)
epsSC = D /pi /Sigma
epsVR = math.gamma (1./6.) /6.**(2./3.) /pi /math.gamma (5./6.) /Pphi**(1./6.) /Sigma
epsDR = math.gamma (0.25) * D**0.5 /2./pi /math.gamma(0.75) / (iotae * Pperp)**0.25 /Sigma

print ("\nepsI = %10.3e epsVI = %10.3e epsDI = %10.3e epsRI = %10.3e epsSC = %10.3e epsVR = %10.3e epsDR = %10.3e\n"
       % (epsI, epsVI, epsDI, epsRI, epsSC, epsVR, epsDR))

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
    fac4 = j0 (0.5*abs(Qi) * tp)

    return (fac1 * fac2 * fac3 * fac4).real

def FunIi (tp, t):
    
    if QE + Qe == 0:
        fac1 = epsI
        fac2 = t - tp
    else:
        fac1 = epsI /( (QE + Qe) * 1j )
        fac2 = 1. - cmath.exp ( (QE + Qe) * 1j * (tp - t) )
    fac3 = cmath.exp (- (QE + 0.5*Qi) * 1j * tp )
    fac4 = j0 (0.5*abs(Qi) * tp)

    return (fac1 * fac2 * fac3 * fac4).imag

def FunVIr (tp, t):

    fac1 = epsVI /math.gamma (1.25)
    fac2 = tp**0.25 * cmath.exp (- (QE + Qe) * 1j * tp)

    return (fac1 * fac2).real

def FunVIi (tp, t):

    fac1 = epsVI /math.gamma (1.25)
    fac2 = tp**0.25 * cmath.exp (- (QE + Qe) * 1j * tp)

    return (fac1 * fac2).imag

def FunDIr (tp, t):

    fac1 = epsDI /math.gamma(1.5) 
    fac2 = tp**0.5 * cmath.exp ( - (QE + Qe) * 1j * tp)

    return (fac1 * fac2).real

def FunDIi (tp, t):

    fac1 = epsDI /math.gamma(1.5) 
    fac2 = tp**0.5 * cmath.exp ( - (QE + Qe) * 1j * tp)

    return (fac1 * fac2).imag

def FunRI1r (tpp, tp):

    fac1 = (tp - tpp + 1.e-15)**(-0.25) * (tpp + 1.e-15)**(-0.25)
    fac2 = cmath.exp ( 1j * (Qe - Qi/2.) * tpp )
    fac3 = jv (-0.25, 0.5*abs(Qi) * tpp)
    
    return (fac1 * fac2 * fac3).real

def FunRI1i (tpp, tp):

    fac1 = (tp - tpp + 1.e-15)**(-0.25) * (tpp + 1.e-15)**(-0.25)
    fac2 = cmath.exp ( 1j * (Qe - Qi/2.) * tpp )
    fac3 = jv (-0.25, 0.5*abs(Qi) * tpp)
    
    return (fac1 * fac2 * fac3).imag

def FunRIr (tp, t):

    fac1       = epsRI * abs(Qi)**0.25 /(2.*pi)**0.5
    fac2       = cmath.exp (- 1j * (QE + Qe) * tp)
    fac3, errr = quad (FunRI1r, 0., tp, args = (tp))
    fac4, erri = quad (FunRI1i, 0., tp, args = (tp))

    return (fac1 * fac2 * (fac3 + 1j * fac4)).real

def FunRIi (tp, t):

    fac1       = epsRI * abs(Qi)**0.25 /(2.*pi)**0.5
    fac2       = cmath.exp (- 1j * (QE + Qe) * tp)
    fac3, errr = quad (FunRI1r, 0., tp, args = (tp))
    fac4, erri = quad (FunRI1i, 0., tp, args = (tp))

    return (fac1 * fac2 * (fac3 + 1j * fac4)).imag

def FunSCr (tp, t):

    fac1 = epsSC /math.gamma(0.5) /(QE + Qe) /1j
    fac2 = 1. - cmath.exp ( (QE + Qe) * 1j * (tp - t))
    fac3 = cmath.exp ( - QE * 1j * tp) /(tp + 1.e-15)**0.5

    return (fac1 * fac2 * fac3).real

def FunSCi (tp, t):

    fac1 = epsSC /math.gamma(0.5) /(QE + Qe) /1j
    fac2 = 1. - cmath.exp ( (QE + Qe) * 1j * (tp - t))
    fac3 = cmath.exp ( - QE * 1j * tp) /(tp + 1.e-15)**0.5

    return (fac1 * fac2 * fac3).imag

def PsiVRr (t):

    fac1 = epsVR / ( (QE + Qe) * 1j + epsVR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsVR ) * t)

    return (fac1 * fac2).real

def PsiVRi (t):

    fac1 = epsVR / ( (QE + Qe) * 1j + epsVR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsVR ) * t)

    return (fac1 * fac2).imag

def PsiDRr (t):

    fac1 = epsDR / ( (QE + Qe) * 1j + epsDR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsDR ) * t)

    return (fac1 * fac2).real

def PsiDRi (t):

    fac1 = epsDR / ( (QE + Qe) * 1j + epsDR )
    fac2 = 1. - cmath.exp ( - ( (QE + Qe) * 1j + epsDR ) * t)

    return (fac1 * fac2).imag

# ############################
# Calculate analytic solutions
# ############################

tt = np.linspace (0., tmax, Nt)

with open ("Analytic.out", "w") as f:

    for t in tt:
    
        PsiIr,  errr = quad (FunIr,  0., t, args = (t))
        PsiIi,  erri = quad (FunIi,  0., t, args = (t))
        PsiVIr, errr = quad (FunVIr, 0., t, args = (t))
        PsiVIi, erri = quad (FunVIi, 0., t, args = (t))
        PsiSCr, errr = quad (FunSCr, 0., t, args = (t))
        PsiSCi, erri = quad (FunSCi, 0., t, args = (t))
        PsiDIr, errr = quad (FunDIr, 0., t, args = (t))
        PsiDIi, erri = quad (FunDIi, 0., t, args = (t))
        if t == 0.:
            PsiRIr = 0.
            PsiRIi = 0.
        else:
            PsiRIr, errr = quad (FunRIr, 0., t, args = (t))
            PsiRIi, errr = quad (FunRIi, 0., t, args = (t))

        print ("%9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e"
               % (t, PsiIr, PsiIi, PsiVIr, PsiVIi, PsiDIr, PsiDIi, PsiRIr, PsiRIi, PsiSCr, PsiSCi, PsiDRr (t), PsiDRi (t), PsiVRr (t), PsiVRi (t)))

        f.write ("%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n"
                 % (t, PsiIr, PsiIi, PsiVIr, PsiVIi, PsiDIr, PsiDIi, PsiRIr, PsiRIi, PsiSCr, PsiSCi, PsiDRr (t), PsiDRi (t), PsiVRr (t), PsiVRi (t)))

f.close ()    


        
