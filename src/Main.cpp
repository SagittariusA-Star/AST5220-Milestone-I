#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"



void example_make_2D_spline(){
 
  // The range and number of points
  const double xmin = 0.0;
  const double xmax = 2.0*M_PI;
  const int    nx   = 49;
  const double ymin = 0.0;
  const double ymax = 2.0*M_PI;
  const int    ny   = 101;

  // A test function to generate some data with
  // and its derivatives (to test the spline derivative routines)
  auto function = [&](double x, double y){
    return x + y + x*y;
  };

  // Make an x-array, an y-array and a z = z(x,y) array
  Vector x_array = Utils::linspace(xmin, xmax, nx);
  Vector y_array = Utils::linspace(ymin, ymax, ny);
  Vector z_array = Vector(nx*ny);
  
  // Fill the z array
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++)
      z_array[ix + nx * iy] = function(x_array[ix], y_array[iy]);

  // Make a spline
  Spline2D z_spline(x_array, y_array, z_array, "Test 2D spline");

  // Test that it gives good results
  std::cout << "Example getting function values from 2D spline\n";
  std::cout << "# x   y  z_spline(x,y)   z_exact(x,y)\n";
  for(int ix = 0; ix < 4; ix++){
    for(int iy = 0; iy < 4; iy++){
      double x = (ix+1.0)/4.0;
      double y = (iy+1.0)/4.0;
      std::cout << std::setw(8) << x             << " ";
      std::cout << std::setw(8) << y             << " ";
      std::cout << std::setw(8) << z_spline(x,y) << " ";
      std::cout << std::setw(8) << function(x,y) << "\n";
    }
  }
  
  //=================================================================
  //=================================================================

  // We have the same functionality as for the 1D spline
  // and we can compute up to second order derivatives using
  // z_spline.deriv_x(x,y)
  // z_spline.deriv_y(x,y)
  // z_spline.deriv_xx(x,y)
  // z_spline.deriv_yy(x,y)
  // z_spline.deriv_xy(x,y)
  
  // Test that it gives good results
  std::cout << "Example getting derivatives from 2D spline\n";
  std::cout << "# x   y  ddzdxdy_spline(x,y)   ddzdxdy_exact(x,y)\n";
  for(int ix = 0; ix < 4; ix++){
    for(int iy = 0; iy < 4; iy++){
      double x = (ix+1.0)/4.0;
      double y = (iy+1.0)/4.0;
      std::cout << std::setw(8) << x             << " ";
      std::cout << std::setw(8) << y             << " ";
      std::cout << std::setw(8) << z_spline.deriv_xy(x,y) << " ";
      std::cout << std::setw(8) << 1.0           << "\n";
    }
  }
  
  //=================================================================
  //=================================================================
  
  // Alternative way of making a spline: from a 2D vector z[ix][iy]
  Vector2D zz_array(nx, Vector(ny));
  for(int ix = 0; ix < nx; ix++)
    for(int iy = 0; iy < ny; iy++){
      zz_array[ix][iy] = function(x_array[ix], y_array[iy]);
    }
  Spline2D zz_spline(x_array, y_array, zz_array, "Test 2D spline");

  //=================================================================
  //=================================================================
  Vector x_data = Utils::linspace(xmin, xmax, 1e3);
  for (int i = 0; i < 1e3; i++){
    std::cout << z_spline(x_data[i], M_PI) << std::endl;
  }

}


int main(int argc, char **argv){
  Utils::StartTiming("Everything");
  
  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  /*
  double h           = 0.7;
  double OmegaB      = 0.046;
  double OmegaCDM    = 0.224;
  double OmegaLambda = 0.72995;
  double Neff        = 3.046;
  double TCMB        = 2.7255;
  */
  // Debug Background parameters
  double h           = 0.7;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.45;
  double OmegaLambda = 0.5;
  double Neff        = 3.046;
  double TCMB        = 2.7255;
  

  // Recombination parameters
  double Yp          = 0.0;

  //=========================================================================
  // Module I
  //=========================================================================
  example_make_2D_spline();
  return 0;
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
  double kvalue = 0.01 / Constants.Mpc;
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
