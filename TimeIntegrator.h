// Responsible for integrating in time
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include "ExactNozzle.h"
#include "EulerOperator.h"
#include "DataManager.h"

class EulerExplicit {
  //double xmin,xmax;
  MeshGenBASE* mesh;
  EulerBASE* euler;
  const double CFL;
  double Density_min = 1.0e-7;
  double Velocity_min = 1.0e-7;
  double Pressure_min = 1.0e-7;

  double Density_max = 1.0e6;
  double Velocity_max = 1.0e8;
  double Pressure_max = 1.0e10;

  public:
  EulerExplicit(MeshGenBASE* &m,EulerBASE* &e,const double &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,4>>* &field); // Outputting the local time step for every cell in the domain
  vector<double> ComputeGlobalTimeStep(vector<array<double,4>>* &field); // Outputting the local time step for every cell in the domain

  void FWDEulerAdvance(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega); //Computing the new solution at the next time step

  void SolutionLimiter(vector<array<double,4>>* &field); //Re-assigning primitive variables to specified max and min limits

  void UnderRelaxationCheck(array<double,4> ResidPrevNorm,array<double,4> ResidNorm,double C,array<bool,4> &check); //looks for residual norm increase to mark needed under-relaxation for each equation 

  bool CheckStallResids(int &count,array<double,4> &ResidNorms,array<double,4> &Prev_ResidualNorms,SpaceVariables2D* &sol); //checks if residuals are stalled

  ~EulerExplicit();


};

//TODO: Add RUNGE-KUTTA Classes (2nd and 4th stage)



#endif
