#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
double h, 
double OmegaB, 
double OmegaCDM,
double Neff, 
double TCMB) :
h(h),
OmegaB(OmegaB),
OmegaCDM(OmegaCDM),
Neff(Neff), 
TCMB(TCMB)
{

	H0 		= Constants.H0_over_h * h; // s^-1 Hubble constant 
	OmegaR 	= pow(M_PI, 3) / 15 * pow(Constants.k_b * TCMB, 4) / (pow(Constants.hbar, 3)
    	       * pow(Constants.c, 5)) * 8 * M_PI * Constants.G / (3 * H0 * H0);       // Radiation density parameter today
	OmegaLambda = 1 - OmegaB - OmegaR - OmegaCDM;
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
Utils::StartTiming("Eta");

Vector x_array;
x_array = Utils::linspace(x_start, x_end, npts); // Log-scale factors 

// The ODE for deta/dx
ODEFunction detadx = [&](double x, const double *eta, double *detadx){
	detadx[0] = Constants.c / Hp_of_x(x);
return GSL_SUCCESS;
};

// Solving ODE for eta(x)
double eta_init = 0;
Vector eta_ic{eta_init};

ODESolver ode;
ode.solve(detadx, x_array, eta_ic);

auto all_data = ode.get_data();
Vector eta_values(npts);

// Filling array with eta values
for (int i = 0; i < all_data.size(); i++){
  	eta_values[i] = all_data[i][0];
}


// Createing a spline of eta(x)
eta_of_x_spline.create(x_array, eta_values, "Spline of eta(x)"); 

Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  	// Hubble parameter as a function of the log-scale factor x
  	double H = H0 * sqrt( (OmegaB + OmegaCDM) * exp(-3 * x)
                      + OmegaR * exp(-4 * x)
                      + OmegaLambda );
  	return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
	// Expansion rate as a function of the log-scale factor
  	double H = BackgroundCosmology::H_of_x(x) * exp(x);
  	return H;
}

double BackgroundCosmology::alpha(double x) const{
  	// Factor used for dHpdx_of_x and ddHpddx_of_x
  
  	return ( - (OmegaB + OmegaCDM) * exp(-x)
          	- 2 * OmegaR * exp(-2 * x)
          	+ 2 * OmegaLambda * exp(2 * x) );
}

double BackgroundCosmology::beta(double x) const{
  	// Factor used for ddHpddx_of_x
  
  	return ( (OmegaB + OmegaCDM) * exp(-x)
          	+ 4 * OmegaR * exp(-2 * x)
          	+ 4 * OmegaLambda * exp(2 * x) );
}

double BackgroundCosmology::dHpdx_of_x(double x) const{
	// Derivative of Hubble parameter wrt and as a function of the log-scale factor
  	return 0.5 * H0 * H0 * alpha(x) / Hp_of_x(x);
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{
	// Double derivative of Hubble parameter wrt and as a function of the log
	// scale factor
  	return 0.5 * H0 * H0 * beta(x) / Hp_of_x(x)
         	- 0.25 * pow(H0, 4) * alpha(x) * alpha(x) / pow(Hp_of_x(x), 3);
}

double BackgroundCosmology::get_OmegaB(double x) const{ 
	// Baryon density parameter as a function of the log-scale factor
  	double Hubble_frac_sq = H_of_x(x) * H_of_x(x) / (H0 * H0);
	if(x == 0.0) return OmegaB;

  	else{
    	return OmegaB * exp(-3 * x) / Hubble_frac_sq;
  	}
}

double BackgroundCosmology::get_OmegaR(double x) const{ 
  	// Radiation density parameter as a function of the log-scale factor
  	double Hubble_frac_sq = H_of_x(x) * H_of_x(x) / (H0 * H0);

  	if(x == 0.0) return OmegaR;

  	else{
		return OmegaR * exp(-4 * x) / Hubble_frac_sq;
  	}
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  	// CDM density parameter as a function of the log-scale factor
  	double Hubble_frac_sq = H_of_x(x) * H_of_x(x) / (H0 * H0);

  	if(x == 0.0) return OmegaCDM;

  	else{
    	return OmegaCDM * exp(-3 * x) / Hubble_frac_sq;
  	}
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  	// Dark energy density parameter as a function of the log-scale factor
  	double Hubble_frac_sq = H_of_x(x) * H_of_x(x) / (H0 * H0);

  	if(x == 0.0) return OmegaLambda;
	
	return OmegaLambda / Hubble_frac_sq;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB() const{ 
  return TCMB; 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  //std::cout << "OmegaK:      " << OmegaK      << "\n";
  //std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = Constants.x_start;
  const double x_max =  Constants.x_end;
  const int    n_pts =  npts;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);
 
  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << ddHpddx_of_x(x)      << " ";
    fp <<"\n";
    //fp << get_OmegaNu(x)     << " ";
    //fp << get_OmegaK(x)      << " ";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

