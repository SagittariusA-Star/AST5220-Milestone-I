#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");
  //===================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //=================================================================== 
  Vector k_array(n_k);
  const double dk = (log10(Constants.k_max) - log10(Constants.k_min)) / (n_k - 1.0);
  const double c = Constants.c;
  double Hp;
  const double H0       = cosmo -> get_H0();
  const double OmegaR0  = cosmo -> get_OmegaR(0);
  double dtaudx;
  int len_tc;       // Length of tight coupling
  int x_full_index; // Index translating current index of x_full to that of x_all.
  Vector y_tight_coupling (Constants.n_ell_tot_tc);
  Vector x_all = Utils::linspace(x_start, x_end, n_x);  // All times, full + tight coupling.
  Vector x_tc;                  // Times for tight coupling 
  Vector x_full;                // Times for full solution
  Vector Psi (n_x * n_k);       // Time-time metric perturbation
  Vector Phi (n_x * n_k);       // Space-space metric perturbation
  Vector delta_cdm (n_x * n_k); // Dark matter density contrast
  Vector delta_b (n_x * n_k);   // Baryon density contrast
  Vector v_cdm (n_x * n_k);     // Dark matter velocity perturbation
  Vector v_b (n_x * n_k);       // Baryon velocity perturbation
  Vector Delta_cdm (n_x * n_k); // Co-moving dark matter density contrast 

  Vector2D Thetas = Vector2D(Constants.n_ell_theta, Vector(n_x * n_k));
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){
    k_array[ik] = log10(Constants.k_min) + ik * dk;
    k_array[ik] = pow(10, k_array[ik]);
    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k);

    for (int ix = 0; ix < n_x; ix++){
      if (x_all[ix] >= x_end_tight){
        len_tc = ix;
        x_tc = Utils::linspace(x_start, x_all[ix - 1], ix);
        x_full = Utils::linspace(x_all[ix - 1], x_end, n_x - ix + 1);
        break;
      }
    }
    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){      
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    ODESolver ode_tc;

    // Solving ODE and extracting solution array from ODEsolver
    ode_tc.solve(dydx_tight_coupling, x_tc, y_tight_coupling_ini);
    auto all_data_tc = ode_tc.get_data();

    for (int jx = 0; jx < len_tc; jx++){
      Hp                        = cosmo -> Hp_of_x(x_tc[jx]);
      dtaudx                    = rec   -> dtaudx_of_x(x_tc[jx]);
      
      Phi[jx + n_x * ik]        = all_data_tc[jx][Constants.ind_Phi_tc];
      delta_cdm[jx + n_x * ik]  = all_data_tc[jx][Constants.ind_deltacdm_tc];
      delta_b[jx + n_x * ik]    = all_data_tc[jx][Constants.ind_deltab_tc];
      v_cdm[jx + n_x * ik]      = all_data_tc[jx][Constants.ind_vcdm_tc];
      v_b[jx + n_x * ik]        = all_data_tc[jx][Constants.ind_vb_tc];
      Thetas[0][jx + n_x * ik]  = all_data_tc[jx][Constants.ind_start_theta_tc];
      Thetas[1][jx + n_x * ik]  = all_data_tc[jx][Constants.ind_start_theta_tc + 1];
      Thetas[2][jx + n_x * ik]  = - 20.0 * c * k / (45.0 * Hp * dtaudx) * Thetas[1][jx + n_x * ik];
        
      for (int ell = 3; ell < Constants.n_ell_theta; ell++){
          Thetas[ell][jx + n_x * ik] = - ell / (2.0 * ell + 1.0) * c * k / (Hp * dtaudx) * Thetas[ell - 1][jx + n_x * ik];
        }

      Psi[jx + n_x * ik] = - Phi[jx + n_x * ik]; 
                           - 12.0 * H0 * H0 / (c * c * k * k * exp(2 * x_tc[jx]))
                           * OmegaR0 * Thetas[2][jx + n_x * ik];
      Delta_cdm[jx + n_x * ik] = delta_cdm[jx + n_x * ik] - 3 * Hp / (c * k) * v_cdm[jx + n_x * ik];
      
    }

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    y_tight_coupling = ode_tc.get_data_by_xindex(len_tc - 1);
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };
    ODESolver ode_full;

    // Solving ODE and extracting solution array from ODEsolver
    ode_full.solve(dydx_full, x_full, y_full_ini);
    auto all_data_full = ode_full.get_data();
    
    x_full_index = 0;
    for (int jx = len_tc - 1; jx < n_x; jx++){
      Hp                        = cosmo -> Hp_of_x(x_full[x_full_index]);

      Phi[jx + n_x * ik]        = all_data_full[x_full_index][Constants.ind_Phi];
      delta_cdm[jx + n_x * ik]  = all_data_full[x_full_index][Constants.ind_deltacdm];
      delta_b[jx + n_x * ik]    = all_data_full[x_full_index][Constants.ind_deltab];
      v_cdm[jx + n_x * ik]      = all_data_full[x_full_index][Constants.ind_vcdm];
      v_b[jx + n_x * ik]        = all_data_full[x_full_index][Constants.ind_vb];
        
      for (int ell = 0; ell < Constants.n_ell_theta; ell++){
          Thetas[ell][jx + n_x * ik] = all_data_full[x_full_index][Constants.ind_start_theta + ell];
        }

      Psi[jx + n_x * ik] = - Phi[jx + n_x * ik]; 
                           - 12.0 * H0 * H0 / (c * c * k * k * exp(2 * x_full[x_full_index]))
                           * OmegaR0 * Thetas[2][jx + n_x * ik];
      Delta_cdm[jx + n_x * ik] = delta_cdm[jx + n_x * ik] - 3 * Hp / (c * k) * v_cdm[jx + n_x * ik];
      x_full_index++;
    }
  }
  Utils::StartTiming("integrateperturbation");
  
  delta_cdm_spline.create(x_all, k_array, delta_cdm, "delta_cdm spline");
  Delta_cdm_spline.create(x_all, k_array, Delta_cdm, "delta_cdm spline");
  delta_b_spline.create(x_all, k_array, delta_b, "delta_b spline");
  v_cdm_spline.create(x_all, k_array, v_cdm, "v_cdm spline");
  v_b_spline.create(x_all, k_array, v_b, "v_b spline");
  Phi_spline.create(x_all, k_array, Phi, "Phi spline");
  Psi_spline.create(x_all, k_array, Psi, "Psi spline");
  for (int ell = 0; ell < Constants.n_ell_theta; ell++){
    Theta_spline[ell].create(x_all, k_array, Thetas[ell], "Theta ell");
  }  
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  
  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];

  // Defining convinient short-hand notation for some quantities
  double Hp       = cosmo -> Hp_of_x(x);
  double dtaudx   = rec   -> dtaudx_of_x(x);
  double c        = Constants.c;
  double ckHp     = c * k / Hp;

  // Initial conditions for metric perturbations in tight coupling
  double Psi_init = - 2.0 / 3.0;
  Phi             = - Psi_init;
  
  // Initial conditions for matter perturbations in tight coupling
  delta_cdm       = - 3.0 / 2.0 * Psi_init;
  delta_b         = delta_cdm;
  v_cdm           = - ckHp / 2.0 * Psi_init;
  v_b             = v_cdm;

  // Initial conditions of the fort two moments of radiation perturbations in tight coupling
  Theta[0] = - 0.5 * Psi_init;
  Theta[1] = ckHp / 6.0 * Psi_init;
  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm];
  double &delta_b         =  y[Constants.ind_deltab];
  double &v_cdm           =  y[Constants.ind_vcdm];
  double &v_b             =  y[Constants.ind_vb];
  double &Phi             =  y[Constants.ind_Phi];
  double *Theta           = &y[Constants.ind_start_theta];

  // Defining convinient short-hand notation for some quantities
  double Hp     = cosmo -> Hp_of_x(x);
  double dtaudx = rec   -> dtaudx_of_x(x);
  double c      = Constants.c;
  double ckHp   = c * k / Hp;
  
  // Initial condition of metric perturbation after tight coupling
  Phi       = Phi_tc;

  // Initial condions for matter perturbations after tight coupling
  delta_cdm = delta_cdm_tc;
  delta_b   = delta_b_tc;
  v_cdm     = v_cdm_tc;
  v_b       = v_b_tc;
  
  // Initial conditions for the multipole moments of radiation solved for after tight coupling
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = - 20.0 * ckHp / (45.0 * dtaudx) * Theta[1];
  
  
  for (int ell = 3; ell < Constants.n_ell_theta; ell++){
    Theta[ell] = - ell / (2.0 * ell + 1.0) * ckHp /  dtaudx * Theta[ell - 1];
  }
  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  // Defining variables for some quantities used in determining the time of 
  // transition from tight coupling into full equation system:

  double x_tight_coupling_end;  // Ending time of tight coupling
  double Xe;                    
  double dtaudx;
  double Hp;
  double c   = Constants.c;
  double ckHp;
  int points = 5e3; // Grid size used in determining the transition time
  Vector x   = Utils::linspace(x_start, x_end, points); // Vector of times

  // Looping through time and testing when transition rules are triggered to find
  // transition time.
  for (int i = 0; i < points; i++){
    Xe      = rec   -> Xe_of_x(x[i]);
    dtaudx  = rec   -> dtaudx_of_x(x[i]);
    Hp      = cosmo -> Hp_of_x(x[i]);
    ckHp    = c * k / Hp;

    // Fist two conditions
    if (std::fabs(dtaudx) < 10 * std::min(1.0, ckHp)){
      x_tight_coupling_end = x[i];
      break;
    }
    // Second condition
    if (Xe < 0.999){
      x_tight_coupling_end = x[i];
      break; 
    }
  }
  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");
  // Sourse functionfor later use when computing the CMB angular power spectrum.
  Vector k_array(n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x); // Linearly spaced x values

  // Log-spaced wave numbers
  const double dk = (log10(Constants.k_max) - log10(Constants.k_min)) / (n_k - 1.0);

  for(int ik = 0; ik < n_k; ik++){
    k_array[ik] = log10(Constants.k_min) + ik * dk;
    k_array[ik] = pow(10, k_array[ik]);
  }
  
  const double c = Constants.c; // Light speed

  // Defining quantities needed from cosmo and rec object
  double Hp;
  double dHpdx;
  double ddHpddx;
  double g_tilde;
  double dg_tildedx;
  double ddg_tildeddx;
  double tau;
  
  // Defining  quantetie needed from current script
  double vb;
  double dvbdx;
  double Phi;
  double dPhidx;
  double Psi;
  double dPsidx;
  double Theta0;
  double Theta2;
  double dTheta2dx;
  double ddTheta2ddx;

  // The four terms of the temperature source function
  double term1;
  double term2;
  double term3;
  double term4;

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    Hp           = cosmo -> Hp_of_x(x);
    dHpdx        = cosmo -> dHpdx_of_x(x);
    ddHpddx      = cosmo -> ddHpddx_of_x(x);
    g_tilde      = rec   -> g_tilde_of_x(x);
    dg_tildedx   = rec   -> dgdx_tilde_of_x(x);
    ddg_tildeddx = rec   -> ddgddx_tilde_of_x(x);
    tau          = rec   -> tau_of_x(x);

    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];
      vb = v_b_spline(x, k);
      dvbdx = v_b_spline.deriv_x(x, k);
      Phi = Phi_spline(x, k);
      dPhidx = Phi_spline.deriv_x(x, k);
      Psi = Psi_spline(x, k);
      dPsidx = Psi_spline.deriv_x(x, k); 
      Theta0 = Theta_spline[0](x, k);
      Theta2 = Theta_spline[2](x, k);
      dTheta2dx = Theta_spline[2].deriv_x(x, k);
      ddTheta2ddx = Theta_spline[2].deriv_x(x, k);

      // Index of current sourse function value to save
      const int index = ix + n_x * ik;

      // Computing individual terms
      term1 = g_tilde * (Theta0 + Psi + 0.25 * Theta2);
      term2 = exp(-tau) * (dPsidx - dPhidx);
      term3 = - 1.0 / (c * k) * ( dHpdx * g_tilde * vb
                                 + Hp * dg_tildedx * vb
                                 + Hp * g_tilde * dvbdx);
      term4 = Theta2 * g_tilde * (dHpdx * dHpdx + Hp * ddHpddx)
            + 3.0 * Hp * dHpdx * (dg_tildedx * Theta2 + g_tilde * dTheta2dx)
            + Hp * Hp * (ddg_tildeddx * Theta2 + 2 * dg_tildedx * dTheta2dx + g_tilde * ddTheta2ddx);
      
      // Temperature source function
      ST_array[index] = term1 + term2 + term3 + term4;                      

    }
  }
  // Spline the temperature source function
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];

  // Constants and quantities used for derivatives
  const double H0             = cosmo -> get_H0();
  const double Hp             = cosmo -> Hp_of_x(x);
  const double dHpdx          = cosmo -> dHpdx_of_x(x);
  const double c              = Constants.c; 
  const double ckHp           = c * k / Hp;
  const double OmegaR0        = cosmo -> get_OmegaR(0);
  const double OmegaCDM0      = cosmo -> get_OmegaCDM(0);
  const double OmegaB0        = cosmo -> get_OmegaB(0);
  const double dtaudx         = rec -> dtaudx_of_x(x);
  const double ddtauddx       = rec -> ddtauddx_of_x(x);
  const double R              = 4.0 * OmegaR0 / (3.0 * OmegaB0) * exp(- x);
  double q;

  // Computing derivatived of the ODE system in tightly coupled regime
  double Theta2 = - 20.0 * ckHp / (45.0 * dtaudx) * Theta[1];
  double Psi = - Phi - 12.0 * H0 * H0 / (c * c * k * k * exp(2 * x)) 
                     * OmegaR0 * Theta2;
  dPhidx = Psi - c * c * k * k / (3.0 * Hp * Hp) * Phi
               + H0 * H0 / (2.0 * Hp * Hp) 
               * (OmegaCDM0 * exp(-x) * delta_cdm
               +  OmegaB0 * exp(-x) * delta_b
               + 4 * OmegaR0 * exp(- 2 * x) * Theta[0]);
  
  ddelta_cdmdx = ckHp * v_cdm - 3 * dPhidx;
  dv_cdmdx     = - v_cdm - ckHp * Psi;
  ddelta_bdx   = ckHp * v_b - 3 * dPhidx;
  dThetadx[0]  = - ckHp * Theta[1] - dPhidx;
  q            = - ((1 - R) * dtaudx + (1 + R) * ddtauddx) * (3 * Theta[1] + v_b)
                - ckHp * Psi
                + ckHp * (1 - dHpdx / Hp) * (-Theta[0] + 2 * Theta2) 
                - ckHp * dThetadx[0];
  q           /= (1 + R) * dtaudx + dHpdx / Hp - 1;
  dv_bdx       = (- v_b - ckHp * Psi + R * (q + ckHp * (-Theta[0] + 2 * Theta2) - ckHp * Psi));
  dv_bdx      /= (1 + R);
  dThetadx[1]  = (q - dv_bdx) / 3.0;

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];

  // Constants and quanties needed to compute the derivatives of the full system
  
  // Cosmo quantities
  const double H0             = cosmo -> get_H0();
  const double Hp             = cosmo -> Hp_of_x(x);
  const double dHpdx          = cosmo -> dHpdx_of_x(x);
  const double OmegaR0        = cosmo -> get_OmegaR(0);
  const double OmegaCDM0      = cosmo -> get_OmegaCDM(0);
  const double OmegaB0        = cosmo -> get_OmegaB(0);
  const double eta            = cosmo -> eta_of_x(x);
  //Recombination quantites
  const double dtaudx         = rec   -> dtaudx_of_x(x);
  // Constants and short-hand notations
  const double R              = 4.0 * OmegaR0 / (3.0 * OmegaB0) * exp(- x);
  const double c              = Constants.c; 
  const double ckHp           = c * k / Hp;
  const int ell_max           = Constants.n_ell_theta - 1;
  double q;


  // Metric time component perturbaton
  double Psi = - Phi - 12.0 * H0 * H0 / (c * c * k * k * exp(2 * x)) 
                     * OmegaR0 * Theta[2];
  
  // Derivative of metric space component perturbation
  dPhidx = Psi - c * c * k * k / (3 * Hp * Hp) * Phi
               + H0 * H0 / (2 * Hp * Hp) 
               * (OmegaCDM0 * exp(-x) * delta_cdm
               +  OmegaB0 * exp(-x) * delta_b
               + 4 * OmegaR0 * exp(- 2 * x) * Theta[0]);
  
  // Derivatives of density and velocity of matter perturbation
  ddelta_cdmdx = ckHp * v_cdm - 3 * dPhidx;
  dv_cdmdx     = - v_cdm - ckHp * Psi;
  ddelta_bdx   = ckHp * v_b - 3 * dPhidx;
  dv_bdx       = - v_b - ckHp * Psi + R * dtaudx * (3.0 * Theta[1] + v_b);

  // Derivatives of multipole moments of radiation perturbations
  dThetadx[0]  = - ckHp * Theta[1] - dPhidx;
  dThetadx[1]  =  ckHp / 3.0 * Theta[0] 
                  - 2.0 / 3.0 * ckHp * Theta[2] 
                  + ckHp / 3.0 * Psi 
                  + dtaudx * (Theta[1] + v_b / 3.0);

  for (int ell = 2; ell < ell_max; ell++){   
    dThetadx[ell] =   ell / (2.0 * ell + 1.0) * ckHp * Theta[ell - 1]
                    - (ell + 1.0) / (2.0 * ell + 1.0) * ckHp * Theta[ell + 1]
                    + dtaudx * Theta[ell]; 
    if (ell == 2){
      dThetadx[ell] -= 1.0 / 10.0 * Theta[2] * dtaudx;
    }
  }
  
  dThetadx[ell_max] = ckHp * Theta[ell_max - 1]
                    - (ell_max + 1.0) / (Hp * eta) * c * Theta[ell_max]
                    + dtaudx * Theta[ell_max];

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_Delta_cdm(const double x, const double k) const{
  // Computes co-moving CDM density contrast 
  return Delta_cdm_spline(x, k);
}

double Perturbations::get_delta_cdm(const double x, const double k) const{
  // Computes CDM density contrast
  return delta_cdm_spline(x, k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  // Computes baryon density contrast
  return delta_b_spline(x, k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  // Computes CDM velocity perturbation
  return v_cdm_spline(x, k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  // Computes baryon density perturbartion
  return v_b_spline(x, k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  // Computes space-space metric perturbation
  return Phi_spline(x, k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  // Computes time-time metric perturbation
  return Psi_spline(x, k);
}

double Perturbations::get_Pi(const double x, const double k) const{
  // Computes polarization tensor
  return Pi_spline(x, k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  // Computes source function
  return ST_spline(x, k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  // Computes polarization source function
  return SE_spline(x, k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  // Computes photon multipole moments
  return Theta_spline[ell](x, k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  // Computes photon polarization multipole moments 
  return Theta_p_spline[ell](x, k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  // Computes neutrino multipole moments
  return Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 1e4;

  auto x_array = Utils::linspace(x_start, x_end, npts);
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo -> eta_of_x(0.0) - cosmo -> eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x, k, 0)   << " ";
    fp << get_Theta(x, k, 1)   << " ";
    fp << get_Theta(x, k, 2)   << " ";
    fp << get_Phi(x, k)       << " ";
    fp << get_Psi(x, k)       << " ";
    fp << get_delta_cdm(x, k)       << " ";
    fp << get_delta_b(x, k)       << " ";
    fp << get_v_cdm(x, k)       << " ";
    fp << get_v_b(x, k)       << " ";
    fp << get_Source_T(x, k)       << " ";
    fp << get_Delta_cdm(x, k)       << " ";

    //fp << get_Source_T(x,k)  << " ";
    //fp << get_Source_T(x, k) * Utils::j_ell(5,   arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    //fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";
    fp << Utils::j_ell(225, arg)           << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

