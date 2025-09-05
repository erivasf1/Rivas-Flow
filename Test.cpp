// Unit testing file that uses the Catch2 header library
#define CATCH_CONFIG_MAIN
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdarg.h>
#include <catch2/catch.hpp> //for unit testing

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "EulerOperator.h"
#include "DataManager.h"
#include "Output.h"
#include "TimeIntegrator.h"
#include "MeshAccess.hpp"

using namespace std;


//Start using Testcases as cases to test whole classes or segments of the Main file; not individual functions

TEST_CASE(" MeshGen fcns. "){

  //TODO: use a simple cartesian grid for testing
  const char* mesh_file = "Grids/InletGrids/Inlet.53x17.grd";

  MeshGenBASE* mesh = new MeshGen2D(mesh_file);



  //Clean-up of dynamic objects
  delete mesh;

}
//EulerOperator Fcns.

/*TEST_CASE(" EulerOperator " ){

  //Mesh Specifications
  int cellnum = 10.0;
  double xmin = -1.0; double xmax = 1.0;
  vector<double> xcoords;

  //Fluid properties
  double P0 = 10.0 * 1000.0;
  double Pb = 10.0 * 1000.0;
  double T0 = 10.0;
  double gamma = 1.4; //specific heat ratio
  bool cond{false}; //flow condition (false = supersonic & true = subsonic)

  // Flux Specifications
  const bool flux_scheme{false}; //true for JST Damping & false for Upwind
  const bool upwind_scheme{false}; //true for Van Leer & false for Rhoe
  const bool flux_accuracy{true}; //true for 1st order & false for 2nd order

  //Field variables
  vector<array<double,3>> Field(cellnum); //stores primitive variable sols.
  vector<array<double,3>> ExpectedField(cellnum); //stores primitive variable sols.

  //Pointers to Field variables
  vector<array<double,3>>* field = &Field; //pointer to Field solutions
  vector<array<double,3>>* expected_field = &ExpectedField; //pointer to Field solutions

  //Object Initializations
  MeshGen1D Mesh(xmin,xmax,cellnum); //mesh

  SpaceVariables1D Sols; //for operating on Field variables

  Tools tool; //utilities object

  Euler1D Euler(cellnum,P0,Pb,T0,gamma); //for performing Euler Eq. operations 

  //Pointers to Objects
  MeshGen1D* mesh = &Mesh;
  Euler1D* euler = &Euler;

  //Generating 1D Mesh
  Mesh.GenerateMesh(xcoords);


  double expected_density,expected_velocity,expected_pressure;

  //Initialize Field with trivial solutions
  for (int i=0;i<(int)Field.size();i++){
    for (int j=0;j<3;j++)
      Field[i][j] = 10.0*j + 1; 
  }  

  //Testing location
  int loc = cellnum/2; int nbor = (cellnum/2) + 1;

  SECTION( "Roe-Avg. Computation" ){ //Verified
    
    //Testing fcn.
    double abar_test;
    array<double,3> tested_roeavg = euler->ComputeRoeAvgVars(field,loc,nbor,abar_test);
      
    //Expected fcn.
    array<double,3> expected_roeavg;
    double rho_loc = (*field)[loc][0];double rho_nbor = (*field)[nbor][0];
    double u_loc = (*field)[loc][1];double u_nbor = (*field)[nbor][1];
    double p_loc = (*field)[loc][2];double p_nbor = (*field)[nbor][2];
    double ht_loc = (gamma/(gamma-1.0))*(p_loc/rho_loc); //pressure work
    ht_loc += pow(u_loc,2.0) / 2.0; //kinetic energy
    double ht_nbor = (gamma/(gamma-1.0))*(p_nbor/rho_nbor); //pressure work
    ht_nbor += pow(u_nbor,2.0) / 2.0; //kinetic energy

    double R_ihalf = sqrt(rho_nbor/rho_loc);

    expected_roeavg[0] = sqrt(rho_loc*rho_nbor);
    //expected_roeavg[0] = R_ihalf * rho_loc;
    expected_roeavg[1] = ((R_ihalf*u_nbor) + u_loc) / (R_ihalf + 1.0);
    expected_roeavg[2] = ((R_ihalf*ht_nbor) + ht_loc) / (R_ihalf + 1.0);
    
    for (int i=0;i<3;i++) {
      REQUIRE(tested_roeavg[i] == Approx(expected_roeavg[i]));
    }
    

  }

  double abar;
  array<double,3> roe_avg = euler->ComputeRoeAvgVars(field,loc,nbor,abar);

  SECTION( "Roe Eigen-values" ){ //Verified

    //Tested fcn.
    array<double,3> test_eigenvals = euler->ComputeRoeEigenVals(roe_avg,abar);

    //Expected fcn.
    array<double,3> expected_eigenvals;
    double ubar = roe_avg[1];
    expected_eigenvals[0] = ubar;
    expected_eigenvals[1] = ubar + abar;
    expected_eigenvals[2] = ubar - abar;

    for (int i=0;i<3;i++) {
      REQUIRE(test_eigenvals[i] == Approx(expected_eigenvals[i]));
    }
    


  }


  SECTION( "Roe Eigen-vectors" ){  //Verified

    //Tested fcn.
    array<array<double,3>,3> test_eigenvecs = euler->ComputeRoeEigenVecs(roe_avg,abar);

    //Expected fcn.
    array<array<double,3>,3> expected_eigenvectors;

    double rho_bar = roe_avg[0];
    double u_bar = roe_avg[1];
    double ht_bar = roe_avg[2];

    double r2_scale = rho_bar / (2.0*abar);
    double r3_scale = -rho_bar / (2.0*abar);
    //1st eigenvector
    expected_eigenvectors[0][0] = 1.0; 
    expected_eigenvectors[0][1] = u_bar; 
    expected_eigenvectors[0][2] = pow(u_bar,2.0)*0.5; 

    //2nd eigenvector
    expected_eigenvectors[1][0] = 1.0; 
    expected_eigenvectors[1][1] = u_bar + abar; 
    expected_eigenvectors[1][2] = ht_bar + (u_bar*abar);

    //3rd eigenvector
    expected_eigenvectors[2][0] = 1.0; 
    expected_eigenvectors[2][1] = u_bar - abar; 
    expected_eigenvectors[2][2] = ht_bar - (u_bar*abar);

    //Scaling 2nd and 3rd eigenvectors
    for (int n=0;n<3;n++)
      expected_eigenvectors[1][n] *= r2_scale;

    for (int n=0;n<3;n++)
      expected_eigenvectors[2][n] *= r3_scale;
    
    for (int i=0;i<3;i++){
      for (int j=0;j<3;j++)
      REQUIRE(test_eigenvecs[i][j] == Approx(expected_eigenvectors[i][j]));
    }

  }

  SECTION( "Roe Wave Amps" ){ //Verified

    //Tested fcn.
    array<double,3> test_waveamps = euler->ComputeRoeWaveAmps(roe_avg,field,abar,loc,nbor);

    //Expected fcn.
    array<double,3> expected_waveamps;
    double drho = (*field)[nbor][0] - (*field)[loc][0];
    double du = (*field)[nbor][1] - (*field)[loc][1];
    double dP = (*field)[nbor][2] - (*field)[loc][2];

    double rho_bar = roe_avg[0];
    double u_bar = roe_avg[1];
    double ht_bar = roe_avg[2];

    expected_waveamps[0] = drho - (dP/pow(abar,2.0));
    expected_waveamps[1] = du + (dP/(rho_bar*abar));
    expected_waveamps[2] = du - (dP/(rho_bar*abar));

    for (int i=0;i<3;i++)
    REQUIRE(test_waveamps[i] == Approx(expected_waveamps[i]));


  }

  array<double,3> eigenvals = euler->ComputeRoeEigenVals(roe_avg,abar);
  array<array<double,3>,3> eigenvecs = euler->ComputeRoeEigenVecs(roe_avg,abar);
  array<double,3> wave_amps = euler->ComputeRoeWaveAmps(roe_avg,field,abar,loc,nbor);

  SECTION( "Roe Flux" ){ 
    //Test fcn.
    array<double,3> test_roeflux = euler->ComputeRoeFlux(field,loc,nbor);

    //Expected fcn.
    array<double,3> expected_roeflux;
    array<double,3> flux_loc = euler->ComputeSpatialFlux_CELL(field,loc);
    array<double,3> flux_nbor = euler->ComputeSpatialFlux_CELL(field,nbor);
    array<double,3> shock_tube;
   
    //shock tube portion
    for (int i=0;i<3;i++){
      for (int j=0;j<3;j++) 
        shock_tube[i] = 0.5*abs(eigenvals[i])*wave_amps[i]*eigenvecs[i][j];
    }

    for (int i=0;i<3;i++)
      expected_roeflux[i] = 0.5*(flux_loc[i]+flux_nbor[i]) - shock_tube[i];



    for (int i=0;i<3;i++)
      REQUIRE(test_roeflux[i] == Approx(expected_roeflux[i]));

  }

}*/


/*TEST_CASE(" 2D Mesh Functions " ){

  SECTION(" Mesh Indexing ") {
  
    int Nx = 3;
    int Ny = 3;

    vector<array<double,4>> Domain(Nx*Ny); 
    vector<array<double,4>>* domain = &Domain; 
 
    array<double,4> exact_val{1.0,2.0,3.0,4.0};
    (*domain)[4] = exact_val; //!< setting cell id 4 to test_val
 
    array<double,4> test_val = fieldij(domain,1,1,Nx);
    for (int i=0;i<4;i++)
      REQUIRE(test_val[i] == Approx(exact_val[i]));
 
  }
}
*/

