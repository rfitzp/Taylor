// Taylor.h

// ###############################################
// Class to solve viscid, two-fluid Taylor problem
// ###############################################

#pragma once

#define _CRT_SECURE_NO_DEPRECATE
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

#ifdef _WIN32
 #include <direct.h>
 #define mkdir _mkdir
#else
 #include <sys/stat.h>
 #include <sys/types.h>
#endif

#include <nlohmann/json.hpp>
#include <gsl/gsl_spline.h>

using namespace std;
using json = nlohmann::json;

// ############
// Class header
// ############
class Taylor
{
private:

  // ..................
  // Physics parameters
  // ..................
  double QE;    // Normalized ExB frequency (read from namelist)
  double Qe;    // Normalized electron diamagnetic frequency (read from namelist)
  double Qi;    // Normalized ion diamagnetic frequency (read from namelist)
  double D;     // Normalized ion sound radius (read from namelist)
  double Pphi;  // Normalized momentum diffusivity (read from namelist)
  double Pperp; // Normalized energy diffusivity (read from namelist)
  double iotae; // Ratio of electron diamagnetic frequency to total diamagnetic frequency 
  double Sigma; // S^1/3 /(-E_ss), where S is Lundquist number, and E_ss is tearing stability index (read from namelist)

  // .....................
  // Simulation parameters
  // .....................
  double tmax;  // Maximum value of normalized time (read from namelist)
  int    Nt;    // Number of equally spaced times between 0 and tmax at which to evaluate reconnected magnetic flux (read from namelist)
  
  // ....................................
  // Inverse Laplace transform parameters
  // ....................................
  double            sigma;    // Real part of g on Bromwich contour (read from namelist)
  double            omax;     // Bromwich contour runs from omega = -omax to omega = +omax (read from namelist)
  int               No;       // Number of grid points on Bromwich contour (read from namelist)
  double            omega;    // Imaginary part of g on Bromwich contour
  vector<double>    gg_i;     // Imaginary part of g on Bromwich contour
  vector<double>    Psib_r;   // Real part of Laplace transformed reconnected flux on Bromwich contour
  vector<double>    Psib_i;   // Imaginary part of Laplace transformed reconnected flux on Bromwich contour
  gsl_spline*       spline_r; // Spline interpolator for Psib_r
  gsl_spline*       spline_i; // Spline interpolator for Psib_i
  gsl_interp_accel* acc_r;    // Accelerator for Psib_r
  gsl_interp_accel* acc_i;    // Accelerator for Psib_i
    
  // .........................
  // Layer equation parameters
  // .........................
  double            pstart;    // Layer equations integrated from p = pstart to p = pend (read from namelist)
  double            pend;      // Layer equations integrated from p = pstart to p = pend (read from namelist)
  int               Np;        // Number of points on p grid (read from namelist)
  gsl_spline*       spline_Vr; // Spline interpolator for V_r
  gsl_spline*       spline_Vi; // Spline interpolator for V_i
  gsl_interp_accel* acc_Vr;    // Accelerator for V_r
  gsl_interp_accel* acc_Vi;    // Accelerator for V_i

  // ................
  // RK4/5 parameters
  // ................
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ...............................
  // Adaptive integration parameters
  // ...............................
  double acc;         // Integration accuracy
  double h0;          // Initial step size
  double hmin;        // Minimum step length
  double hmax;        // Maximum step length
  int    maxrept;     // Maximum step of recursion
  int    flag;        // Error calculation flag
  int    rhs_chooser; // Right-hand side choice switch
  
  // Misc
  int    count;
  double t;           // Normalized time
  complex<double> Im; // Square-root of -1: Im = complex<double> (0., 1.)
  
public:

  // Constructor
  Taylor ();

  // Solve problem
  void Solve ();

private:
  
  // Solve layer equations to obtain layer parameters
  tuple <complex<double>, complex<double>> GetLayerParameters ();
  
  // Evaluate right-hand sides of differential equations
  void Rhs (double x, vector<complex<double>>& y, vector<complex<double>>& dydx);
  // Adaptive step length RK4/RK5 integration routine
  void RK4RK5Adaptive (double& x, vector<complex<double>>& y, double& h,
		       double& t_err, double acc, double S, double T, int& rept,
		       int maxrept, double h_min, double h_max, int flag);
  // Fixed step length RK4/RK5 integration routine
  void RK4RK5Fixed (double& x, vector<complex<double>>& y, vector<complex<double>>& err, double h);

  // Read JSON namelist file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (const char* filename);
  // Check that directory exists and create it otherwise
  bool CreateDirectory (const char* path);
};
