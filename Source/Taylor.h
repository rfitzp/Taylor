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
  double QE;    // Normalized ExB frequency (read from JSON file)
  double Qe;    // Normalized electron diamagnetic frequency (read from JSON file)
  double Qi;    // Normalized ion diamagnetic frequency (read from JSON file)
  double D;     // Normalized ion sound radius (read from JSON file)
  double Pphi;  // Normalized momentum diffusivity (read from JSON file)
  double Pperp; // Normalized energy diffusivity (read from JSON file)
  double iotae; // Ratio of electron diamagnetic frequency to total diamagnetic frequency 
  double Sigma; // S^1/3 /(-E_ss), where S is Lundquist number, and E_ss is tearing stability index (read from JSON file)

  // .........................
  // Layer equation parameters
  // .........................
  double            pstart;    // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double            pend;      // Layer equations integrated from p = pstart to p = pend (read from JSON file)
  double            P3max;     // Value of Pmax[3] above which switch to low-D layer equations made (read from JSON file)
  int               Np;        // Number of points on p grid (read from JSON file)
  gsl_spline*       spline_Vr; // Spline interpolator for V_r
  gsl_spline*       spline_Vi; // Spline interpolator for V_i
  gsl_interp_accel* acc_Vr;    // Accelerator for V_r
  gsl_interp_accel* acc_Vi;    // Accelerator for V_i
  
  // ....................................
  // Inverse Laplace transform parameters
  // ....................................
  double            sigma;    // Real part of g on Bromwich contour (read from JSON file)
  double            omax;     // Bromwich contour runs from omega = -omax to omega = +omax (read from JSON file)
  int               No;       // Number of grid points on Bromwich contour (read from JSON file)
  double            omega;    // Imaginary part of g on Bromwich contour
  vector<double>    gg_i;     // Imaginary part of g on Bromwich contour
  vector<double>    Psib_r;   // Real part of Laplace transformed reconnected flux on Bromwich contour
  vector<double>    Psib_i;   // Imaginary part of Laplace transformed reconnected flux on Bromwich contour
  gsl_spline*       spline_r; // Spline interpolator for Psib_r
  gsl_spline*       spline_i; // Spline interpolator for Psib_i
  gsl_interp_accel* acc_r;    // Accelerator for Psib_r
  gsl_interp_accel* acc_i;    // Accelerator for Psib_i

  // .....................
  // Simulation parameters
  // .....................
  double tmax;  // Maximum value of normalized time (read from JSON file)
  int    Nt;    // Number of equally spaced times between 0 and tmax at which to evaluate reconnected magnetic flux (read from JSON file)
    
  // ................
  // RK4/5 parameters
  // ................
  double aa1, aa2, aa3, aa4, aa5, aa6, cc1, cc3, cc4, cc6, ca1, ca3, ca4, ca5, ca6;
  double bb21, bb31, bb32, bb41, bb42, bb43, bb51, bb52, bb53, bb54;
  double bb61, bb62, bb63, bb64, bb65;

  // ...............................
  // Adaptive integration parameters
  // ...............................
  double acc;         // Integration accuracy (read from JSON file)
  double h0;          // Initial step size (read from JSON file)
  double hmin;        // Minimum step length
  double hmax;        // Maximum step length (read from JSON file)
  int    maxrept;     // Maximum step of recursion
  int    flag;        // Error calculation flag
  int    rhs_chooser; // Right-hand side choice switch

  // ....
  // Misc
  // ....
  int             count, lowD;
  double          t;           // Normalized time
  complex<double> Im;          // Square-root of -1: Im = complex<double> (0., 1.)
  
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

  // Read JSON input file
  json ReadJSONFile (const string& filename);
  // Open new file for writing
  FILE* OpenFilew (const char* filename);
  // Check that directory exists and create it otherwise
  bool CreateDirectory (const char* path);
};
