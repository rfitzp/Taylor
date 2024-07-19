#include "Taylor.h"

// ###########
// Constructor
// ###########
Taylor::Taylor ()
{
  // ......................
  // Set physics parameters
  // ......................
  QE    = 0.01;
  Qe    = 1.;
  Qi    = -1.;
  D     = 3.02;
  Pphi  = 874.;
  Pperp = 287.;
  le    = 1.;
  Sigma = 10.;

  // .........................
  // Set simulation parameters
  // .........................
  tmax  = 300.;
  Nt    = 5000;

  // ...............................
  // Set Bromwich contour parameters
  // ...............................
  sigma = 0.01;
  omax  = 50.;
  No    = 5000;

  // .....................
  // Set p grid parameters
  // .....................
  pstart = 6.;
  pend   = 1.e-8;
  Np     = 5000;

  // ...................................
  // Set adaptive integration parameters
  // ...................................
  acc     = 1.e-12;
  h0      = 1.e-2;
  hmin    = 1.e-10;
  hmax    = 1.e-1;
  maxrept = 50;
  flag    = 2;

  // ....................
  // Set RK4/5 parameters
  // ....................
  aa1 = 0.;
  aa2 = 1. /5.;
  aa3 = 3. /10.;
  aa4 = 3. /5.;
  aa5 = 1.;
  aa6 = 7. /8.;
  
  cc1 = 37.  /378.;
  cc3 = 250. /621.;
  cc4 = 125. /594.;
  cc6 = 512. /1771.;
  
  ca1 = cc1 - 2825.  /27648.;
  ca3 = cc3 - 18575. /48384.;
  ca4 = cc4 - 13525. /55296.;
  ca5 =     - 277.   /14336.;
  ca6 = cc6 - 1.     /4.;
  
  bb21 = 1. /5.;
  
  bb31 = 3. /40.;
  bb32 = 9. /40.;
  
  bb41 =   3. /10.;
  bb42 = - 9. /10.;
  bb43 =   6. /5.;
  
  bb51 = - 11. /54.;
  bb52 =    5. /2.;
  bb53 = - 70. /27.;
  bb54 =   35. /27.;
  
  bb61 = 1631.  /55296.;
  bb62 = 175.   /512.;
  bb63 = 575.   /13824.;
  bb64 = 44275. /110592.;
  bb65 = 253.   /4096.;

  // ...........................
  // Set miscelaneous parameters
  // ...........................
  Im = complex<double> (0., 1.);
}

// #########################
// Function to solve problem
// #########################
void Taylor::Solve ()
{
  // ...............
  // Allocate memory
  // ...............
  gg_i.resize   (No + 1);
  Psib_r.resize (No + 1);
  Psib_i.resize (No + 1);
  
  // ..................................................................
  // Calculate Laplace transformed reconnected flux on Bromwich contour
  // ..................................................................
  FILE* file = OpenFilew ("Plots/Integrand.out");
  printf ("\n");
  for (int j = 0; j <= No; j++)
    {
      // Grid points are crowded around real axis in omega-space
      double fj = - 1. + double (j) * 2. /double (No);
      omega     = omax * fj*fj*fj;
      
      // Calculate Pbar on grid points
      complex<double> Delta, Fs;
      
      tie (Delta, Fs) = GetLayerParameters ();
      
      complex<double> denom = Sigma * Delta + 1.;
      complex<double> Pbar  = Fs /denom /complex<double> (sigma, omega);
      
      // Store Pbar data
      gg_i  [j] = omega;
      Psib_r[j] = real(Pbar);
      Psib_i[j] = imag(Pbar);
      if (j%100 == 0)
	printf ("j = %4d omega = %10.3e Delta = (%10.3e, %10.3e) Fs = (%10.3e, %10.3e) Pbar = (%10.3e, %10.3e)\n",
		j, omega, real(Delta), imag(Delta), real(Fs), imag(Fs), real(Pbar), imag(Pbar));
      
      fprintf (file, "%e %e %e %e %e %e %e %e\n",
	       sigma, omega, real(Delta), imag(Delta), real(Fs), imag(Fs), real(Pbar), imag(Pbar));
    }
  fclose(file);

  // .................................
  // Perform example layer calculation
  // .................................
  omega = 0.;
  complex<double> Delta, Fs;
  tie (Delta, Fs) = GetLayerParameters ();
  
  // .............................................................................
  // Set up interpolation for real and imaginary parts of Pbar on Bromwich contour
  // .............................................................................
  acc_r = gsl_interp_accel_alloc ();
  acc_i = gsl_interp_accel_alloc ();
  
  spline_r = gsl_spline_alloc (gsl_interp_cspline, No + 1);
  spline_i = gsl_spline_alloc (gsl_interp_cspline, No + 1);
  
  gsl_spline_init (spline_r, gg_i.data(), Psib_r.data(), No + 1);
  gsl_spline_init (spline_i, gg_i.data(), Psib_i.data(), No + 1);
  
  // ..............................................................
  // Inverse Laplace transform Laplace-transformed reconnected flux
  // ..............................................................
  file        = OpenFilew ("Plots/Taylor.out");
  FILE* file1 = OpenFilew ("Plots/Diagnostic.out");

  double                  om, h, t_err;
  int                     rept;
  vector<complex<double>> y(1);
  vector<complex<double>> err(1);
  rhs_chooser = 2;
  
  for (int i = 0; i <= Nt; i++)
    {
      t = double (i) * tmax /double (Nt);
      
      om    = - omax;
      h     = h0;
      count = 0;
      y[0]  = 0.;
      
      double hmin_ = 1.e6, hmax_ = 0., err_max = 0.; int reptmax_ = 0, stepcount = 0;
      do
	{
	  RK4RK5Adaptive (om, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag);
	  
	  if (h > hmax_)
	    hmax_ = h;
	  if (h < hmin_)
	    hmin_ = h;
	  if (reptmax_ < rept)
	    reptmax_ = rept;
	  if (err_max < t_err)
	    err_max = t_err;
	  stepcount++;
	}
      while (om + h < omax);
      RK4RK5Fixed (om, y, err, omax - om);
      
      printf ("t = %10.3e  Psi0 = (%10.3e, %10.3e)  h = (%9.2e, %9.2e)  err = %9.2e  rept = %2d  steps = %9.2e\n",
	      t, real(y[0]), imag(y[0]), hmin_, hmax_, err_max, reptmax_, double(stepcount));
      fprintf (file, "%e %e %e\n",
	       t, real(y[0]), imag(y[0]));
      fprintf(file1, "%e %e %e %e %d %e\n",
	      t, log10(hmin_), log10(hmax_), log10(err_max), reptmax_, log10(double(stepcount)));
      fflush(file); fflush(file1);
    }
  fclose(file); fclose(file1);
  
  // ........
  // Clean up
  // ........
  gsl_spline_free (spline_r);
  gsl_spline_free (spline_i);
  
  gsl_interp_accel_free (acc_r);
  gsl_interp_accel_free (acc_i);
}

// #############################################################
// Function to solve layer equations and obtain layer parameters
// #############################################################
tuple <complex<double>, complex<double>> Taylor::GetLayerParameters ()
{
  // ........
  // Define g
  // ........
  complex<double> g = complex<double> (sigma, omega);

  // .....................................................
  // Define data arrays for interpolation of V(p) function
  // .....................................................
  vector<double> pp (Np + 1);
  vector<double> VVr(Np + 1);
  vector<double> VVi(Np + 1);

  // ..............
  // Determine Pmax
  // ..............
  complex<double> gE  = g + Im * QE;
  complex<double> gEe = g + Im * (QE + Qe);
  complex<double> gEi = g + Im * (QE + Qi);
  complex<double> gPD = Pperp + (g + Im * (QE + Qi) * D*D);
  double          PS  = Pphi + Pperp;
  double          PP  = Pphi * Pperp;
  double          PD  = Pphi * D*D /le;
  
  vector<double> Pmax(6);
  Pmax[0] = pow (abs (gEe),          0.5);
  Pmax[1] = pow (abs (gEi * PS /PP), 0.5);
  Pmax[2] = pow (abs (gE * gEi /PP), 0.25);
  Pmax[3] = pow (abs (gPD /PD),      0.5);
  Pmax[4] = pow (abs (gEe /PD),      0.25);
  Pmax[5] = pow (abs (PD /PP),       0.25);
  
  double PMAX = 1.;
  for (int i = 0; i < 6; i++)
    if (Pmax[i] > PMAX)
      PMAX = Pmax[i];

  // ............
  // Setup p grid
  // ............
  for (int j = 0; j <= Np; j++)
    {
      // Grid points are crowded around p=0 in p-space
      double pj = double(j) / double(Np);
      pp[j]     = (pend + (pstart * PMAX - pend) * pj*pj*pj);
    }

  // ...................................................................
  // Integrate layer equations backward in p to calculate V(p) and Delta
  // ...................................................................
  complex<double>         alpha, beta, gamma, X;
  double                  p, h, t_err;
  int                     rept; count = 0;
  vector<complex<double>> y(2);
  vector<complex<double>> dydp(2);
  vector<complex<double>> err(2);
  rhs_chooser = 0;
  
  alpha = - gEe;
  beta  = PP /PD;
  gamma = beta * (1. + gEi * PS /PP - gPD /PD);
  X     = (gamma - sqrt (beta) * (1. - sqrt (beta) * alpha)) / (2. * sqrt (beta));
  
  p       = pstart * PMAX;
  h       = - h0;
  y[0]    = X - sqrt (beta) * p*p;
  y[1]    = 1. + X - sqrt (beta) * p*p;
  VVr[Np] = y[1].real();
  VVi[Np] = y[1].imag();
  
  FILE* file1 = OpenFilew ("Plots/Backward.out");
  for (int j = Np - 1; j >= 0; j--)
    {
      do
	{
	  RK4RK5Adaptive (p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag);
	}
      while (p > pp[j]);
      RK4RK5Fixed (p, y, err, pp[j] - p);

      fprintf (file1, "%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
	       p, y[0].real(), y[0].imag(), y[1].real(), y[1].imag(), log10(-h), log10(t_err));
      
      VVr[j] = y[1].real();
      VVi[j] = y[1].imag();
    }
  fclose (file1);
  
  Rhs (p, y, dydp);
  complex<double> Delta = M_PI /dydp[0];

  // .....................................
  // Set up interpolation of V(p) function
  // .....................................
  acc_Vr = gsl_interp_accel_alloc ();
  acc_Vi = gsl_interp_accel_alloc ();
  
  spline_Vr = gsl_spline_alloc (gsl_interp_cspline, Np + 1);
  spline_Vi = gsl_spline_alloc (gsl_interp_cspline, Np + 1);
  
  gsl_spline_init (spline_Vr, pp.data(), VVr.data(), Np + 1);
  gsl_spline_init (spline_Vi, pp.data(), VVi.data(), Np + 1);

  // ......................................................
  // Integrate layer equations forward in p to calculate Fs
  // ......................................................
  h           = h0;
  p           = pp[0];
  count       = 0;
  y[0]        = complex<double> (0., 0.);
  y[1]        = complex<double> (0., 0.);
  rhs_chooser = 1;

  FILE* file2 = OpenFilew ("Plots/Forward.out");
  for (int j = 1; j < Np-1; j++)
    {
      do
	{
	  RK4RK5Adaptive (p, y, h, t_err, acc, 0.95, 2., rept, maxrept, hmin, hmax, flag);
	}
      while (p < pp[j]);
      RK4RK5Fixed (p, y, err, pp[j] - p);
      
      fprintf (file2, "%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
	       p, y[0].real(), y[0].imag(), y[1].real(), y[1].imag(), log10(h), log10(t_err));
    }
  fclose (file2);
  
  complex<double> Fs = (Delta /gEe /M_PI) * y[1];

  // ........
  // Clean up
  // ........
  gsl_spline_free (spline_Vr);
  gsl_spline_free (spline_Vi);
  
  gsl_interp_accel_free (acc_Vr);
  gsl_interp_accel_free (acc_Vi);

  // .......................
  // Return layer parameters
  // .......................
  return make_tuple (Delta, Fs);
}

// ###################################
// Right hand sides of layer equations
// ###################################
void Taylor::Rhs (double x, vector<complex<double>>& y, vector<complex<double>>& dydx)
{
  // ...............................
  // Define layer equation variables
  // ...............................
  complex<double> W = y[0];
  complex<double> V = y[1];
  complex<double> g (sigma, omega);
  double          p = x;

  // .......................
  // Set calculation factors
  // .......................
  double p2 = p*p;
  double p3 = p*p2;
  double p4 = p2*p2;
  double D2 = D*D;
  
  if (rhs_chooser == 0)
    {
      // ............................................................
      // Right-hand sides for backward integration of layer equations
      // ............................................................
      complex<double> gE  = g + Im * QE;
      complex<double> gEe = g + Im * (QE + Qe);
      complex<double> gEi = g + Im * (QE + Qi);
      complex<double> gPD = Pperp + (g + Im * (QE + Qi) * D*D);
      double          PS  = Pphi + Pperp;
      double          PP  = Pphi * Pperp;
      double          PD  = Pphi * D*D /le;
      
      complex<double> A = p2 /(gEe + p2); 
      complex<double> B = gE * gEi + gEi * PS * p2 + PP * p4;
      complex<double> C = gEe + gPD * p2 + PD * p4;
      
      complex<double> AA = (gEe - p2) /(gEe + p2);
      complex<double> BB = (gEi * PS + 2.*PP * p2) /(gE * gEi + gEi * PS * p2 + PP * p4);
      complex<double> CC = (gPD + 2. * PD * p2) /(gEe + gPD * p2 + PD * p4);
      
      dydx[0] = - AA * W /p - W*W /p + B * p3 /A /C;
      dydx[1] = 2. * p * (BB - CC) * V + 3. * V /p - V*V /p + B * p3 /A /C;
    }
  else if (rhs_chooser == 1)
    {
      // ...........................................................
      // Right-hand sides for forward integration of layer equations
      // ...........................................................
      complex<double> gEe = g + Im * (QE + Qe);
      double          Vr  = gsl_spline_eval (spline_Vr, p, acc_Vr);
      double          Vi  = gsl_spline_eval (spline_Vi, p, acc_Vi);
      complex<double> V   = complex<double> (Vr, Vi);
      
      dydx[0] = V /p;
      dydx[1] = exp (y[0]);
    }
  else if (rhs_chooser == 2)
    {
      // ......................................................
      // Right-hand side for inverse Laplace transform equation
      // ......................................................
      double          psib_r = gsl_spline_eval (spline_r, x, acc_r);
      double          psib_i = gsl_spline_eval (spline_i, x, acc_i);
      complex<double> Pbar   = complex<double> (psib_r, psib_i);
      
      complex<double> expf = exp (sigma * t) * complex<double> (cos(x * t), sin(x * t));
      
      dydx[0] = expf * Pbar /2./M_PI;
    }
}

// ###################################
// Adaptive step length RK4/5 function
// ###################################
void Taylor::RK4RK5Adaptive (double& x, vector<complex<double>>& y, double& h,
			    double& t_err, double acc, double S, double T, int& rept,
			    int maxrept, double h_min, double h_max, int flag)
{
  // Allocate memory
  int neqns = y.size();
  
  vector<complex<double>> y0(neqns);
  vector<complex<double>> Err(neqns);
  double                  hin = h;
  
  // Save initial data
  double x0 = x;
  for (int i = 0; i < neqns; i++)
    y0[i] = y[i];
  
  // Take RK4/RK5 step 
  RK4RK5Fixed (x, y, Err, h);
  
  // Calculate truncation error
  t_err = 0.;
  double err, err1, err2;
  if (flag == 0)
    {
      // Use absolute truncation error 
      for (int i = 0; i < neqns; i++)
	{
	  err = abs(Err[i]);
	  t_err = (err > t_err) ? err : t_err;
	}
    }
  else if (flag == 1)
    {
      // Use relative truncation error
      for (int i = 0; i < neqns; i++)
	{
	  err = abs(Err[i] /y[i]);
	  t_err = (err > t_err) ? err : t_err;
	}
    }
  else
    {
      // Use mixed truncation error 
      for (int i = 0; i < neqns; i++)
	{
	  err1 = abs(Err[i] /y[i]);
	  err2 = abs(Err[i]);
	  err = (err1 < err2) ? err1 : err2;
	  t_err = (err > t_err) ? err : t_err;
	}
    }
  
  // Prevent small truncation error from rounding to zero
  if (t_err < 1.e-15)
    t_err = 1.e-15;
  
  // Calculate new step length
  double h_est;
  if (acc > t_err)
    h_est = S * h * pow(fabs(acc /t_err), 0.20);
  else
    h_est = S * h * pow(fabs(acc /t_err), 0.25);
  
  // Prevent step length from changing by more than factor T
  if (h_est /h > T)
    h *= T;
  else if (h_est / h < 1. /T)
    h /= T;
  else
    h = h_est;
  
  // Prevent step length from exceeding h_max
  h = (fabs(h) > h_max) ? h_max * h / fabs(h) : h;
  
  // Prevent step length from falling below h_min
  if (fabs(h) < h_min)
    {
      if (h >= 0.)
	h = h_min;
      else
	h = -h_min;
    }
  
  // Check if truncation error acceptable
  if ((t_err <= acc) || (count >= maxrept))
    {
      // If truncation error acceptable take step 
      rept  = count;
      count = 0;
    }
  else
    {
      // If truncation error unacceptable repeat step 
      count++;
      x = x0;
      for (int i = 0; i < neqns; i++)
	y[i] = y0[i];
      RK4RK5Adaptive (x, y, h, t_err, acc, S, T, rept,
		      maxrept, h_min, h_max, flag);
    }
}

// ################################
// Fixed step length RK4/5 function
// ################################
void Taylor::RK4RK5Fixed (double& x, vector<complex<double>>& y, vector<complex<double>>& err, double h)
{
  // Allocate memory
  int neqns = y.size();
  
  vector<complex<double>> dydx(neqns);
  vector<complex<double>> k1(neqns);
  vector<complex<double>> k2(neqns);
  vector<complex<double>> k3(neqns);
  vector<complex<double>> k4(neqns);
  vector<complex<double>> k5(neqns);
  vector<complex<double>> k6(neqns);
  vector<complex<double>> f(neqns);
  
  // First stage
  Rhs (x, y, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k1[i] = h * dydx[i];
      f[i]  = y[i] + bb21 * k1[i];
    }
  
  // Second stage
  Rhs (x + aa2 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k2[i] = h * dydx[i];
      f[i]  = y[i] + bb31 * k1[i] + bb32 * k2[i];
    }
  
  // Third stage
  Rhs (x + aa3 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k3[i] = h * dydx[i];
      f[i]  = y[i] + bb41 * k1[i] + bb42 * k2[i] + bb43 * k3[i];
    }
  
  // Fourth stage
  Rhs (x + aa4 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k4[i] = h * dydx[i];
      f[i]  = y[i] + bb51 * k1[i] + bb52 * k2[i] + bb53 * k3[i] + bb54 * k4[i];
    }
  
  // Fifth stage
  Rhs (x + aa5 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k5[i] = h * dydx[i];
      f[i]  = y[i] + bb61 * k1[i] + bb62 * k2[i] + bb63 * k3[i] + bb64 * k4[i] + bb65 * k5[i];
    }
  
  // Sixth stage
  Rhs (x + aa6 * h, f, dydx);
  for (int i = 0; i < neqns; i++)
    {
      k6[i] = h * dydx[i];
    }
  
  // Actual step 
  for (int i = 0; i < neqns; i++)
    {
      y[i]   = y[i] + cc1 * k1[i] + cc3 * k3[i] + cc4 * k4[i]               + cc6 * k6[i];
      err[i] =        ca1 * k1[i] + ca3 * k3[i] + ca4 * k4[i] + ca5 * k5[i] + ca6 * k6[i];
    }
  x += h;
}

// #####################################
// Function to open new file for writing
// #####################################
FILE* Taylor::OpenFilew (const char* filename)
{
  FILE* file = fopen (filename, "w");
  if (file == NULL)
    {
      printf ("OpenFilew: Error opening data-file: %s\n", filename);
      exit (1);
    }
  return file;
}
