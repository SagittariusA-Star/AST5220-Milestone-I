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
  double OmegaLambda; // = 0.72995;
  double Neff        = 3.046;
  double TCMB        = 2.7255;
  /*
  // Debug Background parameters
  double h           = 0.7;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.45;
  double OmegaLambda; //= 0.5;
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
  
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue;
  kvalue = 0.3 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.1.txt");
  /*
  */
  /*
  kvalue = 5e-5 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k5e-5.txt");
  
  kvalue = 2.85e-4 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k2.85e-4.txt");
  
  kvalue = 1.623e-3 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k1.623e-3.txt");

  kvalue = 9.244e-3 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k9.244e-3.txt");

  kvalue = 5.266e-2 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k5.266e-2.txt");

  kvalue = 3e-1 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k3e-1.txt");
  // Remove when module is completed
  
  */  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert);
  power.solve();
  power.output("cells.txt");
  
  Utils::EndTiming("Everything");
}
