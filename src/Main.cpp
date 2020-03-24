#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.7;
  double OmegaB      = 0.046;
  double OmegaCDM    = 0.224;
  double OmegaLambda = 0.72995;
  double Neff        = 3.046;
  double TCMB        = 2.7255;
  /*
  // Debug Background parameters
  double h           = 0.7;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.45;
  double OmegaLambda = 0.5;
  double Neff        = 3.046;
  double TCMB        = 2.7255;
  */
  

  // Recombination parameters
  double Yp          = 0.0;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  
  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  
  //=========================================================================
  // Module III
  //=========================================================================
  /*
  std::pair<double,double> xrange (Constants.x_start, Constants.x_end);
  double tol = (Constants.x_end - Constants.x_start) / 5000; 
  double tau1 = Utils::binary_search_for_value(cosmo::tau_of_x_spline, 1.0, xrange, tol);
  std::cout << tau1 << std::endl;
  */
  return 0;
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 * Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
