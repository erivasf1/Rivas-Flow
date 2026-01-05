// Responsible for integrating in time
#ifndef _TIMEINTEGRATOR_H_
#define _TIMEINTEGRATOR_H_
#include "ExactNozzle.h"
#include "EulerOperator.h"
#include "DataManager.h"

class EulerExplicitBASE {
  //double xmin,xmax;
  public:
  MeshGenBASE* mesh;
  EulerBASE* euler;
  const double CFL;
  double Density_min = 1.0e-7;
  //double Velocity_min = 1.0e-7;
  double Pressure_min = 1.0e-7;

  double Density_max = 1.0e6;
  double Velocity_max = 1.0e8;
  double Pressure_max = 1.0e10;

  EulerExplicitBASE(MeshGenBASE* &m,EulerBASE* &e,const double &c);
 
  vector<double> ComputeLocalTimeStep(vector<array<double,4>>* &field); // Outputting the local time step for every cell in the domain
  vector<double> ComputeGlobalTimeStep(vector<array<double,4>>* &field); // Outputting the local time step for every cell in the domain

  void FWDEulerAdvance(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega); //Computing the new solution at the next time step

  void SolutionLimiter(vector<array<double,4>>* &field); //Re-assigning primitive variables to specified max and min limits

  void UnderRelaxationCheck(array<double,4> ResidPrevNorm,array<double,4> ResidNorm,double C,array<bool,4> &check); //looks for residual norm increase to mark needed under-relaxation for each equation 

  bool CheckStallResids(int &count,array<double,4> &ResidNorms,array<double,4> &Prev_ResidualNorms,SpaceVariables2D* &sol); //checks if residuals are stalled

  virtual void ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall); 

  virtual ~EulerExplicitBASE();


};

class EulerExplicit : public EulerExplicitBASE {


  public:
  EulerExplicit(MeshGenBASE* &m,EulerBASE* &e,const double &c);

  void ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall); //calls FWDEulerAdvance fcn.

  ~EulerExplicit();

};

//TODO: Add RUNGE-KUTTA Classes (2nd and 4th stage)

class RungeKutta2 : public EulerExplicit {

  double alpha1 = 1.0 / 2.0;
  double alpha2 = 1.0;

  vector<array<double,4>> Field_intermediate_cons; //in conserved variables form
  vector<array<double,4>> Field_intermediate_prim; //in primitive variables form
  vector<array<double,4>> Resid_intermediate;

  vector<array<double,4>>* field_interm_cons;
  vector<array<double,4>>* field_interm_prim;
  vector<array<double,4>>* resid_interm;

  public:
  
  RungeKutta2(MeshGenBASE* &m,EulerBASE* &e,const double &c);
  
  void ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall) override;


  ~RungeKutta2();

};

class RungeKutta4 : public EulerExplicit {

  double alpha1 = 1.0 / 4.0;
  double alpha2 = 1.0 / 3.0;
  double alpha3 = 1.0 / 2.0;
  double alpha4 = 1.0;

  vector<array<double,4>> Field_intermediate_cons; //in conserved variables form
  vector<array<double,4>> Field_intermediate_prim; //in primitive variables form
  vector<array<double,4>> Resid_intermediate;

  vector<array<double,4>>* field_interm_cons;
  vector<array<double,4>>* field_interm_prim;
  vector<array<double,4>>* resid_interm;

  public:
  
  RungeKutta4(MeshGenBASE* &m,EulerBASE* &e,const double &c);
  
  void ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall) override;


  ~RungeKutta4();

};



#endif
