// 2D FVM Euler Eq. Solver - Erick Rivas
//Quasi-1D Nozzle, Euler Eqs. of Converging-to-Diverging Nozzle - Erick Rivas
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdarg.h>

#include "ExactNozzle.h"
#include "MeshAccess.hpp"
#include "MeshGen.h"
#include "EulerOperator.h"
#include "DataManager.h"
#include "Output.h"
#include "TimeIntegrator.h" 

using namespace std;

//For ease of defining parameters
enum CASE_2D { 
  INLET,AIRFOIL1,AIRFOIL2};

enum BOUNDARY_COND {
  INFLOW,OUTFLOW,SLIP_WALL};


int main() {

  //Timing Purposes
  double start_time,stop_time;
  start_time = MPI_Wtime();

  //! INITIALIZATION
  // Scenario
  int scenario = 2; //1 = 1D, 2 = 2D, 3 = 2D MMS
  CASE_2D case_2d = INLET;

  // Constants for 1D case
  double xmin = -1.0;
  double xmax = 1.0;

  // Boundary Conditions Specification
  BOUNDARY_COND top_cond = SLIP_WALL; //TEMP
  BOUNDARY_COND btm_cond = SLIP_WALL;
  BOUNDARY_COND left_cond = INFLOW;
  BOUNDARY_COND right_cond = OUTFLOW;

  [[maybe_unused]]bool cond_loc{false}; //true for subsonic & false for supersonic (FOR EXACT SOL.)
  [[maybe_unused]]bool cond_bc{true}; //true for subsonic & false for supersonic (FOR OUTFLOW BC)

  // Mesh Specifications
  const char* meshfile = "Grids/InletGrids/Inlet.53x17.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  [[maybe_unused]]int cellnum = 100; //recommending an even number for cell face at the throat of nozzle (NOTE: will get reassigned val. if mesh is provided)

  // Temporal Specifications
  const int iter_max = 1;
  int iterout = 1; //number of iterations per solution output
  const double CFL = 0.1; //CFL number (must <= 1 for Euler Explicit integration)
  //const double CFL = 2.9e-4; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{false}; //true = local time stepping; false = global time stepping

  // Flux Specifications
  int flux_scheme{1}; //0=JST, 1=Van Leer, 2 = Roe (this will eventually be used!)
  int flux_accuracy{1}; //1=1st order, 2=2nd order
  [[maybe_unused]] const double ramp_stop = 1.0e-7; //stopping criteria for ramping fcn. of transitioning from 1st to 2nd
  double epsilon = 1.0; //ramping value used to transition from 1st to 2nd order
  bool resid_stall{false};//for detecting if residuals have stalled
  [[maybe_unused]]int stall_count = 0;

  // Under-Relaxation Parameters
  double C = 1.2; //residual norm check
  array<double,4> Omega{1.0,1.0,1.0,1.0}; //FWD Advance Limiter
  int subiter_max = 1e2; //max number of relaxation sub-iterations
  array<bool,4> check{false,false,false,false}; //false by default to check if under-relaxation is needed

  // Governing Eq. Residuals
  double cont_tol = 1.0e-8;
  double xmom_tol = 1.0e-8;
  double ymom_tol = 1.0e-8;
  double energy_tol = 1.0e-8;

  //! GENERATING MESH 
  MeshGenBASE* mesh; 
  if (scenario == 1){ //1D Nozzle Mesh Case
    mesh = new MeshGenNozzle(xmin,xmax,mesh->cellnumber);
    mesh->GenerateMesh();
  }
  else if ((scenario == 2) || (scenario == 3)) {
    mesh = new MeshGen2D(meshfile);
    mesh->OutputMesh(); //!< outputs mesh file for visualization
  }
  else{
    cerr<<"Unknown scenario number!"<<endl;
    return 0;
  }
  
  //! DATA ALLOCATION
  //Field variables
  vector<array<double,4>> Field(mesh->cellnumber); //stores primitive variable sols.
  vector<array<double,4>> FieldStar(mesh->cellnumber); //stores intermediate primitive variable sols.
  vector<array<double,4>> FieldStall(mesh->cellnumber); //stores primitive variable sols. before stall (if detected)
  vector<array<double,4>> FieldMS(mesh->cellnumber); //stores manufactured sol.
  vector<array<double,4>> FieldMS_Source(mesh->cellnumber); //stores manufactured source term for all cells

  vector<array<double,4>> ExactField(mesh->cellnumber); //stores exact cell-averaged primitve variable sols.
  vector<array<double,4>> ExactFaces(mesh->cellnumber+1); //stores exact primitve variable sols. at cell faces

  vector<array<double,4>> Residual(mesh->cellnumber); //stores the local residuals per eq.
  vector<array<double,4>> ResidualStar(mesh->cellnumber); //stores the intermediate stage of primtive variables
  vector<array<double,4>> InitResidual(mesh->cellnumber); //stores the initial residual

  vector<double> TimeSteps(mesh->cellnumber); //for storing the time step (delta_t) for each cell
  array<double,4> ResidualNorms; //for storing the global residual norms
  array<double,4> ResidualStarNorms; //stores the intermediate global residual norms
  array<double,4> Prev_ResidualNorms; //for storing the previous global residual norms

  //Pointers to Field variables
  vector<array<double,4>>* field = &Field; //pointer to Field solutions
  vector<array<double,4>>* field_star = &FieldStar; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_stall = &FieldStall; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_ms = &FieldMS; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_ms_source = &FieldMS_Source; //pointer to intermediate Field solutions
  [[maybe_unused]] vector<array<double,4>>* exact_sols = &ExactField; //pointer to exact solution field values
  [[maybe_unused]] vector<array<double,4>>* exact_faces = &ExactFaces; //pointer to exact solution field values
  vector<array<double,4>>* resid = &Residual; //pointer to residual field values per cell
  vector<array<double,4>>* resid_star = &ResidualStar; //pointer to intermediate residual field values per cell
  vector<array<double,4>>* init_resid = &InitResidual; //pointer to residual field values per cell
  vector<double>* time_steps = &TimeSteps;

  //!OBJECT INITIALIZATIONS

  //Euler Operator spec.
  // Specifying EulerOperator
  EulerBASE* euler;
  //Temp -- will add scenario == 1 once 1D section is fixed!
  if (scenario == 2) 
    euler = new Euler2D(case_2d,mesh->cell_imax,mesh->cell_jmax,flux_scheme,flux_accuracy,top_cond,btm_cond,left_cond,right_cond,mesh,field_ms_source);
  else if (scenario == 3)
    euler = new Euler2DMMS(mesh->cell_imax,mesh->cell_jmax,flux_scheme,flux_accuracy,top_cond,btm_cond,left_cond,right_cond,mesh,field_ms_source);
  else{
    cerr<<"Error: scenario # not recognized!"<<endl;
    return 0;
  }
  
  //Time Integrator spec.
  //TODO: Specifying time integrator via Polymorphism
  EulerExplicit Time(mesh,euler,CFL); //for computing time steps

  SpaceVariables2D Sols; //for operating on Field variables

  Output Error; //for discretization error operations

  //Pointers to Objects
  SpaceVariables2D* sols = &Sols;

  EulerExplicit* time = &Time;
  Output* error = &Error;

  //! PRINTING OUT SIMULATION INFO
  // Title
  if (meshfile)
    Tools::print("2D EULER EQ. SOLVER\n");
  else
    Tools::print("1D EULER EQ. SOLVER\n");
  // Case Spec
  if (meshfile){
    Tools::print("-Mesh Selected: ");
    Tools::print("%s\n",meshfile);
  }
  else{
    Tools::print("-Case Selected: ");
    (cond_bc == true) ? Tools::print("Shock Wave Case\n") : Tools::print("Isentropic Case\n");
  }
  // Spatial Stats
  Tools::print("-Spatial Statistics:\n");
  Tools::print("--Cell Number: %d\n",mesh->cellnumber);
  //Tools::print("--Delta x: %f\n",dx);
  // Temporal Stats
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("--Time-stepping method: ");
  (timestep == true) ? Tools::print("Local time-stepping\n") : Tools::print("Global time-stepping\n");  
  // Temporal Stats
  Tools::print("-Flux Statistics:");
  if (flux_scheme == 1) //Van Leer Flux 
    (flux_accuracy == 1) ? Tools::print(" 1st Order Van Leer Scheme\n") : Tools::print(" 2nd Order Van Leer Scheme\n");
  else if (flux_scheme == 2) //Roe Flux
    (flux_accuracy == 1) ? Tools::print(" 1st Order Roe's Scheme\n") : Tools::print(" 2nd Order Roe's Scheme\n");
  
  else
    Tools::print(" JST Damping\n");


  //! COMPUTING MANUFACTURED SOLUTION AND SOURCE TERMS (MMS ONLY)
  if (scenario == 3){
    string mms_sol_filename = "ManufacturedSols.dat";
    string mms_source_filename = "SourceTerms.dat";
    euler->ManufacturedPrimitiveSols(field_ms,sols); //!< computing manufactured sol.
    euler->EvalSourceTerms(sols); //!< computing manufactured source terms
    error->OutputPrimitiveVariables(field_ms,mms_sol_filename,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);
    error->OutputManufacturedSourceTerms(field_ms_source,mms_source_filename,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  }

  //! SETTING INITIAL CONDITIONS
  //Tools::print("At initial conditions\n");
  euler->SetInitialConditions(field);
  //Field = FieldMS; //for now, setting field to manufactured sol.

  //!TODO: COMPUTING EXACT SOLUTION -- ONLY FOR 1D QUASI-STEADY NOZZLE
  /*
  if ((cond_bc == false) && (!meshfile)){ //Compute Exact Solution if isentropic case is selected
    array<double,3> sol;
    double area;
    for (int i=0;i<(int)exact_faces->size();i++) {
      area = tool.AreaVal(xcoords[i]);
      cond_loc = (xcoords[i] < 0) ? true:false; 
      SuperSonicNozzle Nozzle(area,area_star,stag_pressure,stag_temp,cond_loc);
      Nozzle.ComputeExactSol(sol);
 
      (*exact_faces)[i] = sol; //assigning to exact faces vector
    
    }
    // Computing cell-average sol. for all cells
    sols->ComputeCellAveragedSol(exact_faces,exact_sols,xcoords);
  }
  */


  // SETTING BOUNDARY CONDITIONS
  //generates ghost cells here too
  euler->Setup2DBoundaryConditions(field,error); 

  //time->SolutionLimiter(field_test);


  //!< Outputting initial solutions with BC's
  //const char* filename2 = "InitSolutionswBCs.txt";
  //sols->OutputPrimitiveVariables(field,euler,filename2);

  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,4> InitNorms;

  //if (upwind_scheme == false && flux_accuracy == false) //reverting to 1st order if 2nd order Roe is selected
    //flux_accuracy = true;

  euler->ComputeResidual(init_resid,field);

  //debug: Residual
  const char* resid_file = "InitialResiduals.dat"; 
  error->OutputPrimitiveVariables(init_resid,resid_file,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  InitNorms = sols->ComputeSolutionNorms(init_resid); //computing L2 norm of residuals
  //InitNorms = sols->ComputeSolutionNorms(init_resid); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Y-Momentum:%e\n",InitNorms[2]);
  Tools::print("--Energy:%e\n\n",InitNorms[3]);


  (*resid) = (*init_resid);//!< setting initial residual to intermediate
  ResidualNorms = InitNorms;

  string it,name; //used for outputting file name
  int iter; //iteration number

  //Opening file that stores residuals
  ofstream myresids;
  myresids.open("SolResids.dat");
  myresids<<"variables= \"Iteration num.\" \"Continuity\" \"X-Momentum\" \"Y-Momentum\" \"Energy\""<<endl;
  myresids<<"zone T= "<<"\""<<0<<"\""<<endl;
  myresids<<"DATAPACKING=POINT"<<endl;
  myresids<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  //myresids<<"Iteration"<<"  "<<"Contintuity"<<"  "<<"X-Momentum"<<"  "<<"Energy"<<endl;

  myresids<<0<<"  "<<InitNorms[0]<<"  "<<InitNorms[1]<<"  "<<InitNorms[2]<<"  "<<InitNorms[3]<<endl; //printing out the initial residuals first

  //Printing Primitive vars. for TECPLOT visualization
  std::string filename_totalsols = "AllSolutions.dat";
  error->OutputPrimitiveVariables(field,filename_totalsols,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

  //Assigning Intermediate Field to Initial Field (including residuals)
  (*field_star) = (*field);
  (*resid_star) = (*resid);

  (*resid) = (*init_resid);//!< setting initial residual to intermediate
  ResidualNorms = InitNorms;

  //! BEGIN OF MAIN LOOP
  for (iter=1;iter<iter_max;iter++){

    // For Switching from 1st to 2nd order accurate
/*    if ((upwind_scheme == false) && (flux_accuracy == true) && (iter == 2e5)){ //resetting back to 2nd Order Roe - if Roe 2nd order selected
      flux_accuracy = false;
      iterout = 10;
      subiter_max = 5;
    }
 */   
  
    //! COMPUTE TIME STEP
    // if global time step, chosen then create a vector<double> of the smallest time step
    (*time_steps) = (timestep == true) ? time->ComputeLocalTimeStep(field) : time->ComputeGlobalTimeStep(field);

    //! COMPUTE NEW SOL. VALUES 
    time->FWDEulerAdvance(field_star,resid_star,time_steps,Omega);//TESTING
    time->SolutionLimiter(field_star); //applies solution limiter to all cells (including ghost cells)

    //! ENFORCE BOUNDARY CONDITIONS
    euler->Enforce2DBoundaryConditions(field_star,false);
    //time->SolutionLimiter(field_star); //temporarily reapplying the limiter


    euler->ComputeResidual(resid_star,field_star);
    ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
    time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);


    //! UNDER-RELAXATION CHECK
    if (check[0]==true || check[1] == true || check[2] == true || check[3] == true){ //perform under-relaxation if any of these are true
      for (int j=0;j<subiter_max;j++){
        for (int i=0;i<4;i++) //!< reassigns omega to half of current value if under-relaxation detected
          Omega[i] = (check[i] == true) ?  Omega[i] /= 2.0 : Omega[i] = 1.0;

        (*field_star) = (*field); //resetting primitive variables to previous time step values
        time->FWDEulerAdvance(field_star,resid_star,time_steps,Omega); //advancing intermediate solution w/ under-relaxation factor 
        time->SolutionLimiter(field_star); //temporarily reapplying the limiter
        euler->ComputeResidual(resid_star,field_star);
        ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
        time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);

        //checking if new residuals do not need under-relaxation
        if (check[0]==false && check[1] == false && check[2] == false && check[3] == false){ 
          for (int i=0;i<4;i++) //!< resetting omega to 1
            Omega[i] = 1.0;
          break;
        }

        if (j == subiter_max)
          Tools::print("Under-relaxation Failed!\n");
      }
    }

    // Checking if Residuals are stalled
    /*if (resid_stall == false){ //only checking if haven't already been marked as stalled
      Prev_ResidualNorms = ResidualNorms;
      resid_stall = time->CheckStallResids(stall_count,ResidualNorms,Prev_ResidualNorms,sols);
      if (resid_stall == true)
        FieldStall = Field; //setting previous field sols. before stall
    }
    */


    //Assinging New Time Step Values to Intermediate Values
    (*field) = (*field_star);
    (*resid) = (*resid_star); ResidualNorms = ResidualStarNorms;

    //Computing Ramping Value
    //epsilon = sols->ComputeRampValue(ResidualNorms,InitNorms,ramp_stop);


    //! OUTPUT SOL. IN .DAT FILE EVERY "ITEROUT" STEPS
    if (iter % iterout == 0) {
      it = to_string(iter);
      name = "SolResults/Iteration";
      name += it;
      name += ".dat";
      const char* filename_iter = name.c_str();
      //writing all time-steps in unique files
      error->OutputPrimitiveVariables(field,filename_iter,false,iter,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);
      //writing in 1 file
      error->OutputPrimitiveVariables(field,filename_totalsols,true,iter,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

      //Printing to TECPLOT
      //error->OutputPrimitiveVariables(field,filename_totalsols,true,iter,xcoords);

      //Printing Residual Norms to Screen
      Tools::print("------Iteration #: %d----------\n",iter);
      Tools::print("Continuity:%e\nX-Momentum:%e\nY-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2],ResidualNorms[3]);
      //debug:
      Tools::print("Epsilon: %e\n",epsilon);
      if (resid_stall == true)
        Tools::print("Residuals are detected to be stalled!\n"); //printing message is temp. for now

      // Writing Residuals history to "SolResids.txt" file
      myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

    }


    //! CHECK FOR CONVERGENCE (w/ respect to the intial residual norms)
    if (ResidualNorms[0]/InitNorms[0] <= cont_tol && ResidualNorms[1]/InitNorms[1] <= xmom_tol && 
        ResidualNorms[2]/InitNorms[2] <= ymom_tol && ResidualNorms[3]/InitNorms[3] <= energy_tol)
      break;

  }

  //! FINAL OUTPUT OF SOLUTION

  if (iter==iter_max)
    Tools::print("Failed to converge!\n");

  else {
    Tools::print("\n");
    Tools::print("------------------------------------------------------------\n");
    Tools::print("CONGRATS you converged!\n");
    Tools::print("Continuity: %e\nX-Momentum: %e\nEnergy: %e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);
    const char* filename_final = "ConvergedSolution.dat" ;
    error->OutputPrimitiveVariables(field,filename_final,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);
  }

  //Closing Residuals file
  myresids.close();


  //TODO:! EVALUATE DISCRETIZATION NORMS FOR GRID CONVERGENCE AND PRINT OUT TO FILE
  /*if (cond_bc == false){
    field->erase(field->begin()); field->erase(field->begin()); //!< erasing ghost cells
    field->erase(field->end()); field->erase(field->end());

    vector<array<double,3>> Errors(Field);
    vector<array<double,3>>* errors = &Errors;

    error->DiscretizationErrorNorms(field,exact_sols,errors,sols);
  }*/

  stop_time = MPI_Wtime();
  Tools::print("Elapsed time: %fs\n",stop_time-start_time);

  //! CLEANUP
  if (meshfile){
    delete euler;
    delete mesh;
  }

  return 0;
}
