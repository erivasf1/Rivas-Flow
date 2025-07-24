// Responsible for Euler Eqs. Computations
#ifndef _EULEROPERATOR_H_
#define _EULEROPERATOR_H_
#include "ExactNozzle.h"
#include "MeshGen.h"
#include "DataManager.h"

//Forward Declarations
class SpaceVariables2D;

//TODO: Make BASE Class
class EulerBASE {
  
  int top_cond,btm_cond,left_cond,right_cond; //boundary conds. of domain

  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double gamma = 1.4; //specific heat ratio
  const double MolMass = 28.96; // kg/kmol

  public:
  MeshGenBASE* mesh;

  double R = Ru / MolMass; //specific gas constant

  vector<array<double,4>>* mms_source; //pointer to source terms

  int flux_scheme; //flux scheme id (0=JST, 1=VanLeer, 2=Roe)
  int flux_accuracy; //flux accuracy id (1=1st order, 2 = 2nd order)

  int cell_imax; //NUMBER of cells in the i dir!
  int cell_jmax; //NUMBER of cells in the j dir!

  EulerBASE(int &cell_inum,int &cell_jnum,int &scheme,int &accuracy,int &top,int &btm,int &left,int &right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source);

  double ComputeMachNumber(array<double,4> &sols);
  double GetGamma();

  //Compute Conserved & Primitive Variables
  array<double,4> ComputeConserved(vector<array<double,4>>* &field,int &i,int &j);
  //Initial Conditions
  virtual void SetInitialConditions(vector<array<double,4>>* &field); //Complete (tested)
  //TODO:Boundary Conditions -- refer to this paper:https://arc-aiaa-org.ezproxy.lib.vt.edu/doi/10.2514/3.11983  
  void SetupBoundaryConditions();
  //ComputeLeftBoundaryCondition -- inflow
  //ComputeRightBoundaryCondition -- outflow
  //ComputeTopBoundaryCondition -- freestream
  //ComputeBottomBoundaryCondition -- slip wall
  //Spatial Fluxes
  array<double,4> ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,4>>* field,int loci,int locj,int nbori,int nborj); //1st order upwind schemes
  array<double,4> ComputeSpatialFlux_UPWIND2ndOrder(vector<array<double,4>>* field,int loci,int locj,int nbori,int nborj); //2nd order upwind schemes
  //VanLeer
  array<double,4> VanLeerCompute(array<double,4> &field_state,bool sign);
  double GetC(double M,bool sign); //c value
  double GetAlpha(double M,bool sign); //alpha value
  double GetBeta(double M); //beta value
  double GetVanLeerM(double M,bool sign); //Van Leer MachNumber
  double GetD(double M,bool sign); //D value
  double GetP2Bar(double M,bool sign); //Pressure double bar
  //Roe
  //MUSCL Extrapolation
  //Artificial Dissipation (JST Damping Only)
  //Source Term
  //Residual
  virtual void ComputeResidual(vector<array<double,4>>* &resid,vector<array<double,4>>* &field);
  //MMS
  virtual void ManufacturedPrimitiveSols(vector<array<double,4>>* &field,SpaceVariables2D* &sols);
  virtual void EvalSourceTerms(SpaceVariables2D* &sols); //source terms for all governing equations
  virtual ~EulerBASE();
};

class Euler1D { //TODO: Make this into inherit class
  //vector<double> &xcoords;
  // parameters for 1D nozzle
  double stag_pressure; //stagnation pressure
  double back_pressure; //back pressure
  double stag_temperature; //stagnation temperature
  double gamma;
  const double Ru = 8314.0; // J/(kmol*K) -- universal gas constant   
  const double MolMass = 28.96; // kg/kmol
 
  double Density_max = 1.0e9;
  double Velocity_max = 1.0e9;
  double Pressure_max = 1.0e9;

  public:
  //double dx; //cell thickness
  int interior_cellnum; //holds the cell num (interior only)
  int total_cellnum; //holds the total cell num (including boundary conditions)
  double R = Ru / MolMass; //specific gas constant

  Euler1D(); //empty constructor for unit testing

  Euler1D(int &cellnum,double &P0,double &BP,double &T0,double &g); //constructor for Main file

  // Primitive & Conserved variables fcns.
  array<double,3> ComputeConserved(vector<array<double,3>>* &field,int loc);

  void ComputePrimitive(vector<array<double,3>>* &field,array<double,3> &Conserved,int loc);

  // BOUNDARY + INITIAL CONDITIONS FCNS.
  void SetInitialConditions(vector<array<double,3>>* &field,vector<double> &xcoords); //Complete (tested)
  void SetBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //ADDS ghost cells nodes and computes their values
  void ComputeTotalBoundaryConditions(vector<array<double,3>>* &field,bool &cond); //only computes their values
  void ComputeInflowBoundaryConditions(vector<array<double,3>>* &field); //only computes their values
  void ComputeOutflowBoundaryConditions(vector<array<double,3>>* &field,bool& cond); //only computes their values

  // RESIDUAL FCNS.
  void ComputeResidual(vector<array<double,3>>* &resid,vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,vector<double> &xcoords,double &dx,bool flux_scheme,bool flux_accuracy,bool upwind_scheme,double epsilon,bool &resid_stall);

  // SPATIAL FLUXES FCNS. (INCLUDING SOURCE TERM)
  // Central Difference using Central quadrature fcn.
  array<double,3> ComputeSpatialFlux_CELL(array<double,3> &field_state); //at cell
  array<double,3> ComputeSpatialFlux_BASE(vector<array<double,3>>* &field,int loc,int nbor); //at cell interface
  // Upwind Schemes
  array<double,3> ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,3>>* &field,bool method,int loc,int nbor); //1st order upwind schemes
  array<double,3> ComputeSpatialFlux_UPWIND2ndOrder(vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,bool method,int loc, int nbor,double epsilon,bool &resid_stall); //2nd order upwind schemes
  //VanLeer Fcns.
  array<double,3> VanLeerCompute(array<double,3> &field_state,bool sign); //returns convective+pressure flux of specified state (either right or left)
  double GetC(double M,bool sign); //c value
  double GetAlpha(double M,bool sign); //alpha value
  double GetBeta(double M); //beta value
  double GetVanLeerM(double M,bool sign); //Van Leer MachNumber
  double GetD(double M,bool sign); //D value
  double GetP2Bar(double M,bool sign); //Pressure double bar
  //Roe Fcns.
  array<double,3> ComputeRoeFlux(array<double,3> &field_ltstate,array<double,3> &field_rtstate); //rho-avg eigenvalues
  array<double,3> ComputeRoeWaveAmps(array<double,3> &roe_vars,array<double,3> &field_ltstate,array<double,3> &field_rtstate,double abar); //rho-avg eigenvalues
  array<double,3> ComputeRoeEigenVals(array<double,3> &rho_vars,double abar); //rho-avg eigenvalues
  array<array<double,3>,3> ComputeRoeEigenVecs(array<double,3> &roe_vars,double abar); //rho-avg eigenvectors
  array<double,3> ComputeRoeAvgVars(array<double,3> &field_ltstate,array<double,3> &field_rtstate,double &abar); //rho-avg. vars

  double ComputeSourceTerm(vector<array<double,3>>* &field,int loc,vector<double> &xcoords,double dx);

  // MUSCL extrapolation + Flux Limiters
  array<array<double,3>,2> MUSCLApprox(vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,int loc,int nbor,double epsilon,bool &resid_stall); //outputs the left and right state primitive variables
  array<array<double,3>,2> ComputeBetaLimiter(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor,double beta); //computes limiter using the beta limiter method
  array<double,3> ComputeVanLeerLimiter(array<double,3> &r_vec); //computes limiter using the Van Leer method
  array<array<double,3>,2> ComputeRVariation(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor); //consecutive variation

  array<double,3> ComputeRPlusVariation(vector<array<double,3>>* &field,int loc,int r_nbor,int nbor); //plus part of consecutive variation
  array<double,3> ComputeRMinusVariation(vector<array<double,3>>* &field,int loc,int l_nbor,int nbor); //minus part of consecutive variation

  // ARTIFICIAL DISSIPATON FCNS. (USING JST DAMPENING)
  array<double,3> Compute2ndOrderDamping(vector<array<double,3>>* &field,int loc); // viscous term for shocks (c(2))
  array<double,3> Compute4thOrderDamping(vector<array<double,3>>* &field,int loc); // prevents odd-even decoupling (c(4))
  
  double GetEpsilon2(vector<array<double,3>>* &field,int loc); //sensor that detects "sharp" gradients
  double GetEpsilon4(vector<array<double,3>>* &field,int loc); 
 
  double GetLambda(vector<array<double,3>>* &field,int loc);
  double GetNu(vector<array<double,3>>* &field,int loc); //switching fcn.
  double GetMachNumber(vector<array<double,3>>* &field,int loc); //used for GetLambda fcn.
  
  // SUPPLEMENTAL FCNS. (MAY BE USED FOR OTHER FCNS. OF OTHER CLASSES)
  double GetLambdaMax(vector<array<double,3>>* &field,int loc); //extracts largest eigenvalue for a given cell
  static double GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3> &sol_right); //testing x-velocity for now
  double ComputeMachNumber(array<double,3> &sols); //computes Mach number for solution vars. that are not part of the field -- same formula as GetMachNumber

  ~Euler1D();


};


// EulerOperator Class for 2D Problems
class Euler2D : public EulerBASE {
  

  public:
  double Mach_bc,T_stag,P_stag,alpha; //free-stream and initial conditions, depending on case
  
  Euler2D(int case_2d,int scheme,int accuracy,int cell_inum,int cell_jnum,int top,int btm,int left,int right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source); //constructor determines val. of const. parameters (e.g. freestream Mach #, angle-of-attack)

  void InitSolutions(vector<array<double,4>>* &field,int cellnum);
  void SetInitialConditions(vector<array<double,4>>* &field) override;


  void ComputeResidual(vector<array<double,4>>* &resid,vector<array<double,4>>* &field) override;

  ~Euler2D();

};

// EulerOperator Class for 2D Problems w/ MMS
class Euler2DMMS : public EulerBASE {

  //Source Terms Constants
  //Note: Manufactured Sol. from Mathematica!
  // Using Roy AIAA 2002 paper for constants -- Supersonic Condition
  double L = 1.0;
  double rho0 = 1.0;
  double press0 = 1.0e5;
  double uvel0 = 800.0;
  double vvel0 = 800.0;

  double Pi = M_PI;
  double rhox = 0.15;double rhoy = -0.1;
  double uvelx = 50.0;double uvely = -30.0;
  double vvelx = -75.0;double vvely = 40.0;
  double pressx = 0.2e5;double pressy = 0.5e5;
  double wvel0 = 0.0;

  public:
  Euler2DMMS(int cell_inum,int cell_jnum,int scheme,int accuracy,int top,int btm, int left,int right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source);

  void SetInitialConditions(vector<array<double,4>>* &field) override; 

  void ManufacturedPrimitiveSols(vector<array<double,4>>* &field,SpaceVariables2D* &Sols) override;

  double ContinuitySourceTerm(double x,double y);
  double XMomentumSourceTerm(double x,double y);
  double YMomentumSourceTerm(double x,double y);
  double EnergySourceTerm(double x,double y);
  void EvalSourceTerms(SpaceVariables2D* &sols) override; //source terms for all governing equations
  void ComputeResidual(vector<array<double,4>>* &resid,vector<array<double,4>>* &field) override;
  ~Euler2DMMS();

};

#endif
