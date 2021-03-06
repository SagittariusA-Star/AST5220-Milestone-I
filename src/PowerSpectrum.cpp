#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  // Setting up arrays of wavenumbers k with log-spacing
  
  Vector k_array(n_k);
  
  const double dk = (log10(Constants.k_max) - log10(Constants.k_min)) / (n_k - 1.0);
  for(int ik = 0; ik < n_k; ik++){
    k_array[ik] = log10(Constants.k_min) + ik * dk;
    k_array[ik] = pow(10, k_array[ik]);
  }
  
  Vector log_k_array = log(k_array);
  
  // Generating splines for Bessel functions
  generate_bessel_function_splines();

  // Computing line of sight integral
  line_of_sight_integration(k_array);

  // Solving for temperature power spectrum
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size()); // Vector of Bessel spline objects
  int N         = 1e4;  // Grid resolution
  Vector x_arr  = Utils::linspace(0, 5e3, N); // Vector of arguments
  Vector j_ell_x (N);                         // Vector of Bessel functions
  
  // Creating spline for Bessel functions for later use:

  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    for (int j = 0; j < N; j++){      
      j_ell_x[j] = Utils::j_ell(ell, x_arr[j]);
    }
    j_ell_splines[i].create(x_arr, j_ell_x, "Spline of j_ell(x)");
  }
  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));
  
  // Set up initial conditions for Theta_ell
  Vector Theta_ell_ic{0};
  double k;   // Dummy variable for temporary saving current wavenumber
  int N = 5e2;  // Number of x-values to use
  Vector x_array = Utils::linspace(Constants.x_start, Constants.x_end, N);
  
  for(int ik = 0; ik < k_array.size(); ik++){
    // Store the result for Source_ell(k) in results[ell][ik]
    k = k_array[ik];
    
    for (int i = 0; i < ells.size(); i++){  
      // Setting up integral as ODE
      ODEFunction dTheta_elldx = [&](double x, const double *y, double *dydx){      
        double S      = source_function(x, k);
        double eta    = cosmo -> eta_of_x(x);
        double eta0   = cosmo -> eta_of_x(0);
        double j_ell  = j_ell_splines[i](k * (eta0 - eta));
        dydx[0] = S * j_ell;
        return GSL_SUCCESS;
      };

      ODESolver ode;
      // Setting accuracy parameters
      double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
      ode.set_accuracy(hstart, abserr, relerr);
      
      // Solving ODE and extracting solution array from ODEsolver
      ode.solve(dTheta_elldx, x_array, Theta_ell_ic, gsl_odeiv2_step_rkf45);
      auto Theta_ell_today = ode.get_data_by_xindex(N - 1);
      
      // Saving result to array
      result[i][ik] = Theta_ell_today[0];
    }
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  //thetaT_ell_of_k_spline = std::vector<Spline>(nells);

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert -> get_Source_T(x, k);
  };

  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  // Spline the result and store it in thetaT_ell_of_k_spline
  thetaT_ell_of_k_spline.create(ells, k_array, thetaT_ell_of_k);
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================

Vector PowerSpectrum::solve_for_cell(
  Vector & log_k_array,
  Spline2D & f_ell_spline,
  Spline2D & g_ell_spline){
  const int nells      = ells.size();

  // Setting up array for later saving of results
  Vector result(nells);
  // Initial condition for integral solved as ODE
  Vector Cell_ic{0};

  // Perform integral for each ell
  for (int i = 0; i < nells; i++){
    // Setting up integral as ODE
    ODEFunction dCelldlogk = [&](double logk, const double *y, double *dydx){      
          double k = exp(logk);
          double P_primordial = primordial_power_spectrum(k);
          double Theta_ell_sq = f_ell_spline(ells[i], k) * g_ell_spline(ells[i], k);
          dydx[0] = 4 * M_PI * P_primordial * Theta_ell_sq;
          return GSL_SUCCESS;
        };

    ODESolver ode;

    // Setting up accuracy parameters
    double hstart = 1e-3, abserr = 1e-10, relerr = 1e-10;
    ode.set_accuracy(hstart, abserr, relerr);
    
    // Solving ODE and extracting solution array from ODEsolver
    ode.solve(dCelldlogk, log_k_array, Cell_ic, gsl_odeiv2_step_rkf45);
    auto Cell = ode.get_data_by_xindex(log_k_array.size() - 1);
    // Saving result
    result[i] = Cell[0];
  }
  return result;
}


double PowerSpectrum::primordial_power_spectrum(const double k) const{
  // Function returns primordial power spectrum as function of wavenumber k
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

double PowerSpectrum::get_Delta_M(const double x, const double k) const{
  // Function returns comoving density contrast of matter as function of time x
  // and scale k 
  double Phi = pert -> get_Phi(x, k); // Metric potential perturbation 
  double c = Constants.c;         
  return 2.0 * c * c * k * k * Phi / (3.0 * OmegaM0 * exp(-x) * H0 * H0);
}

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  // Function returns matter power spectrum as a function of time x and scale k [Mpc]. 
  // The returned value is in (Mpc)^3.

  double DeltaM = get_Delta_M(x, k_mpc);
  double P_primordial = primordial_power_spectrum(k_mpc);
  double pofk = std::fabs(DeltaM) * std::fabs(DeltaM) * P_primordial;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}
double PowerSpectrum::get_Theta_ell_of_k(const double k, const double ell) const{
  return thetaT_ell_of_k_spline(ell, k);
}
double PowerSpectrum::get_Theta_ell_of_k_sq(const double k, const double ell) const{
  double Theta_ell = thetaT_ell_of_k_spline(ell, k); 
  return Theta_ell * Theta_ell;
}
//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output(std::string filename) const{
  // Output in standard units of muK^2
  std::ofstream fp(filename.c_str());
  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);
  auto print_data = [&] (const double ell) {
    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline(ell) * normfactor  << " ";
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
 
  std::string filename2 = "matter_power_spec.txt";
  std::ofstream fp2(filename2.c_str());
  Vector k_arr = Utils::linspace(Constants.k_min, Constants.k_max, 1e3);
  double factor;
  double TCMB_factor = 1e6 * cosmo->get_TCMB();
  double TCMB_factor_sq = TCMB_factor * TCMB_factor;
  auto print_data2 = [&] (const double k) {
    factor = (2 * M_PI * M_PI) / (k * k * k);
    fp2 << k * Constants.Mpc / cosmo -> get_h()  << " ";
    fp2 << get_matter_power_spectrum(0, k) * factor * pow(cosmo -> get_h() / Constants.Mpc, 3) << " ";
    fp2 << get_Theta_ell_of_k(k, 10)  << " ";
    fp2 << get_Theta_ell_of_k(k, 25)  << " ";
    fp2 << get_Theta_ell_of_k(k, 50)  << " ";
    fp2 << get_Theta_ell_of_k(k, 100)  << " ";
    
    fp2 << get_Theta_ell_of_k_sq(k, 10) / (k * Constants.Mpc / cosmo -> get_h()) << " ";
    fp2 << get_Theta_ell_of_k_sq(k, 25) / (k * Constants.Mpc / cosmo -> get_h()) << " ";
    fp2 << get_Theta_ell_of_k_sq(k, 50) / (k * Constants.Mpc / cosmo -> get_h()) << " ";
    fp2 << get_Theta_ell_of_k_sq(k, 100) / (k * Constants.Mpc / cosmo -> get_h()) << " ";
    fp2 << "\n";
  };
  std::for_each(k_arr.begin(), k_arr.end(), print_data2);
}

