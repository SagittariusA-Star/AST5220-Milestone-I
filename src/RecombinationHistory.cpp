#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================    

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array(npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector ne_arr(npts_rec_arrays);

  x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  
  const double OmegaB = cosmo -> get_OmegaB(0);
  
  // Calculate recombination history
  double nH_current;
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;
    nH_current              = rho_crit0 * OmegaB / Constants.m_H * exp(-3 * x_array[i]);
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;      

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      //...
      //...
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;
    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...

      
      Vector Xe_ic{Xe_arr[i-1]};
      Vector x_current{x_array[i-1], x_array[i]};
      ODESolver ode;
      ode.solve(dXedx, x_current, Xe_ic, gsl_odeiv2_step_rkf45);
      
      auto all_data = ode.get_data();

      Xe_arr[i] = all_data[1][0];
      ne_arr[i] = nH_current * Xe_arr[i];
      
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  Vector log_Xe = log(Xe_arr);
  Vector log_ne = log(ne_arr);

  log_Xe_of_x_spline.create(x_array, log_Xe, "Spline of log Xe(x)");
  log_ne_of_x_spline.create(x_array, log_ne, "Spline of log ne(x)");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  //const double OmegaB      = cosmo->get_OmegaB();
  //...
  //...

  const double OmegaB = cosmo -> get_OmegaB(0);
  
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;
  double nH = OmegaB * rho_crit0 / (Constants.m_H * exp(3 * x)) ;
  double Xe_fraction;

  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  //...
  //...
  double Tb    = kBT(x);
  Xe_fraction  = 1 / nH * (m_e * Tb) / (2 * M_PI * hbar * hbar)
                        * sqrt((m_e * Tb) / (2 * M_PI * hbar * hbar))
                        * exp(-epsilon_0 / Tb);
  if (Xe_fraction > 1e10 ){
    Xe = 1;
  }
  
  else {
    Xe = 0.5 * (- Xe_fraction + sqrt(Xe_fraction * Xe_fraction + 4 * Xe_fraction));
  }

  ne = Xe * nH;

  return std::pair<double,double>(Xe, ne);
}

double RecombinationHistory::kBT(double x) const{
  // Computing baryon temperature energy as function of x
  return Constants.k_b * Constants.TCMB * exp(- x);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...
  double OmegaB       = cosmo -> get_OmegaB(0);
  double H            = cosmo -> H_of_x(x);        
  double Tb           = kBT(x);
  
  /*
  const double alpha  = sqrt(3 * m_e * m_e * c * c * sigma_T 
                               / (8 * M_PI * hbar * hbar));
  */
  double nH           = OmegaB * rho_crit0 / (Constants.m_H * exp(3 * x));
  double phi_2        = 0.448 * log(epsilon_0 / Tb);
  double alpha_2      = 8 / sqrt(3 * M_PI) * sigma_T * c
                          * sqrt(epsilon_0 / Tb) * phi_2;
  double beta_2       = alpha_2 * (m_e * Tb / (2 * M_PI * hbar * hbar))
                                * sqrt(m_e * Tb / (2 * M_PI * hbar * hbar))
                                * exp(-epsilon_0 / (3 * Tb));
  double beta         = alpha_2 * (m_e * Tb / (2 * M_PI * hbar * hbar))
                                * sqrt(m_e * Tb / (2 * M_PI * hbar * hbar))
                                * exp(-epsilon_0 / Tb);
  double n1s          = (1 - X_e) * nH; 
  double lambda_alpha = H * (3 * epsilon_0) * (3 * epsilon_0) * (3 * epsilon_0)
                            / (64 * M_PI * M_PI * n1s * hbar * hbar * hbar * c * c * c);
  double Cr           = (lambda_2s1s + lambda_alpha) / (lambda_2s1s + lambda_alpha + beta_2); 
  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  dXedx[0] = Cr / H * ( beta * (1 - X_e) - nH * alpha_2 * X_e * X_e);

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 10000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Electron density from spline
    double ne = ne_of_x(x);
    // Hubble parameter from BackgroundComplogy class
    double H = cosmo -> H_of_x(x);
    // Set the derivative for photon optical depth
    dtaudx[0] = - ne * Constants.sigma_T * Constants.c / H;
    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...
  
  // Solving ODE for tau(x)
  double tau_init = 1e5;
  Vector tau_ic{tau_init};

  ODESolver ode_tau;
  //double hstart = 1e-15, abserr = 1e-20, relerr = 1e-20;
  //ode_tau.set_accuracy(hstart, abserr, relerr);

  ode_tau.solve(dtaudx, x_array, tau_ic, gsl_odeiv2_step_rkf45);

  auto all_data_tau = ode_tau.get_data();
  Vector tau_values(npts);
  Vector dtaudx_values(npts);
  Vector ddtaudxdx_values(npts);
  double _dtaudx;
  double _tau;
  
  // Filling array with tau values
  for (int i = 0; i < all_data_tau.size(); i++){
      tau_values[i] = all_data_tau[i][0] - all_data_tau[npts - 1][0];
      dtaudx(x_array[i], &tau_values[i], &_dtaudx);
      dtaudx_values[i] = _dtaudx;
  }

  // Createing a spline of eta(x)
  tau_of_x_spline.create(x_array, tau_values, "Spline of tau(x)"); 
  dtaudx_of_x_spline.create(x_array, dtaudx_values, "Spline of dtaudx(x)"); 

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Vector g_tilde_values(npts);
  for (int i = 0; i < npts; i++){
      _tau    = tau_of_x_spline(x_array[i]);
      _dtaudx = dtaudx_of_x_spline(x_array[i]);
      g_tilde_values[i] = - _dtaudx * exp(- _tau);
  }
  g_tilde_of_x_spline.create(x_array, g_tilde_values, "Spline of g_tilde(x)");
  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  //return tau_of_x_spline.deriv_x(x);
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...
  //return tau_of_x_spline.deriv_xx(x);
  return dtaudx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";
    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";
    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

