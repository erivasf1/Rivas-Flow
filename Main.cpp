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

enum CASE_2D { //For ease of defining parameters
  INLET,AIRFOIL1,AIRFOIL2
};

enum BOUNDARY_COND {
  INFLOW,OUTFLOW,SLIP_WALL
};


int main() {

  //Timing Purposes
  double start_time,stop_time;
  start_time = MPI_Wtime();

  //! INITIALIZATION
  // Scenario
  int scenario = 3; //1 = 1D, 2 = 2D, 3 = 2D MMS
  CASE_2D case_2d = AIRFOIL1;
  // Constants
  double xmin = -1.0;
  double xmax = 1.0;
  double stag_pressure = 300.0 * 1000.0; //kPa -> Pa
  double back_pressure = 120.0 * 1000.0; //kPa -> Pa (for subsonic outflow cond.)
  double stag_temp = 600.0; //K
  double gamma = 1.4; //specific heat ratio
  double area_star = Tools::AreaVal(0.5*(xmin+xmax)); //area at throat

  // Boundary Conditions Specification
  BOUNDARY_COND top_cond = OUTFLOW;
  BOUNDARY_COND btm_cond = INFLOW;
  BOUNDARY_COND left_cond = INFLOW;
  BOUNDARY_COND right_cond = OUTFLOW;

  bool cond_loc{false}; //true for subsonic & false for supersonic (FOR EXACT SOL.)
  bool cond_bc{true}; //true for subsonic & false for supersonic (FOR OUTFLOW BC)

  // Mesh Specifications
  int cellnum = 100; //recommending an even number for cell face at the throat of nozzle (NOTE: will get reassigned val. if mesh is provided)
  vector<double> xcoords; //stores the coords of the cell NODES!!! (i.e. size of xcoords is cellnum+1)!
  vector<double> ycoords; //stores the coords of the cell NODES!!! (i.e. size of xcoords is cellnum+1)!
  double dx;
  const char* meshfile = "Grids/CurvilinearGrids/curv2d9.grd"; //name of 2D file -- Note: set to NULL if 1D case is to be ran
  //const char* meshfile = NULL;

  // Temporal Specifications
  const int iter_max = 1e6;
  int iterout = 500; //number of iterations per solution output
  const double CFL = 0.1; //CFL number (must <= 1 for Euler Explicit integration)
  //const double CFL = 2.9e-4; //CFL number (must <= 1 for Euler Explicit integration)
  bool timestep{false}; //true = local time stepping; false = global time stepping

  // Flux Specifications
  bool flux_scheme{false}; 
  int flux_scheme1{1}; //0=JST, 1=Van Leer, 2 = Roe (this will eventually be used!)
  bool upwind_scheme{false}; //true for Van Leer & false for Roe
  bool flux_accuracy{false}; 
  int flux_accuracy1{1}; //1=1st order, 2=2nd order
  [[maybe_unused]] const double ramp_stop = 1.0e-7; //stopping criteria for ramping fcn. of transitioning from 1st to 2nd
  double epsilon = 1.0; //ramping value used to transition from 1st to 2nd order
  bool resid_stall{false};//for detecting if residuals have stalled
  int stall_count = 0;

  // Under-Relaxation Parameters
  double C = 1.2; //residual norm check
  array<double,3> Omega{1.0,1.0,1.0}; //FWD Advance Limiter
  //int subiter_max = 0; //max number of relaxation sub-iterations
  int subiter_max = 1e2; //max number of relaxation sub-iterations

  // Governing Eq. Residuals
  double cont_tol = 1.0e-10;
  double xmom_tol = 1.0e-10;
  double energy_tol = 1.0e-10;

  //! GENERATING MESH -- TODO: Employ polymorphism here too!
  MeshGenBASE* mesh; 
  if (scenario == 1){ //1D Nozzle Mesh Case
    mesh = new MeshGenNozzle(xmin,xmax,cellnum);
    mesh->GenerateMesh();
  }
  else if ((scenario == 2) || (scenario == 3)) {
    mesh = new MeshGen2D(meshfile);
    mesh->OutputMesh(); //!< outputs mesh file for Tecplot visualization
  }
  else{
    cerr<<"Unknown scenario number!"<<endl;
    return 0;
  }
  

  //MeshGen1D* mesh = &Mesh; //1D Mesh

  /*if (meshfile){ //2D Mesh Case -- read from file
    xcoords = mesh_2d->xcoords;
    ycoords = mesh_2d->ycoords;
    //zcoords = mesh_2d->zcoords;
    cellnum = mesh_2d->cellnumber;
    mesh_2d -> OutputMesh(); //outputs mesh for tecplot visualization
  }
  else{ //1D Mesh Case -- Uniform
    mesh->GenerateMesh(xcoords); //stores all coords in xcoords list
    dx = abs(xcoords[0]-xcoords[1]); //delta x distance
  }
  //debug:
  //Tools::print("Size of xcoords: %d\n",xcoords.size());
  //Tools::print("Size of ycoords: %d\n",ycoords.size());
  //Tools::print("imax: %d & jmax: %d\n",mesh_2d->imax,mesh_2d->jmax);
  */

  //! DATA ALLOCATION
  //TODO: Change Field variable array size to 4
  //Field variables
  vector<array<double,3>> Field(cellnum); //stores primitive variable sols.
  vector<array<double,3>> FieldStar(cellnum); //stores intermediate primitive variable sols.
  vector<array<double,3>> FieldStall(cellnum); //stores primitive variable sols. before stall (if detected)
  vector<array<double,4>> FieldTest(mesh->cellnumber); 
  vector<array<double,4>> FieldMS(mesh->cellnumber); //stores manufactured sol.
  vector<array<double,4>> FieldMS_Source(mesh->cellnumber); //stores manufactured source term for all cells

  vector<array<double,3>> ExactField(cellnum); //stores exact cell-averaged primitve variable sols.
  vector<array<double,3>> ExactFaces(cellnum+1); //stores exact primitve variable sols. at cell faces

  vector<array<double,3>> Residual(cellnum); //stores the local residuals per eq.
  vector<array<double,3>> ResidualStar(cellnum); //stores the intermediate stage of primtive variables
  vector<array<double,3>> InitResidual(cellnum); //stores the initial residual

  vector<double> TimeSteps; //for storing the time step (delta_t) for each cell
  array<double,3> ResidualNorms; //for storing the global residual norms
  array<double,3> Prev_ResidualNorms; //for storing the previous global residual norms

  //Pointers to Field variables
  vector<array<double,3>>* field = &Field; //pointer to Field solutions
  vector<array<double,3>>* field_star = &FieldStar; //pointer to intermediate Field solutions
  vector<array<double,3>>* field_stall = &FieldStall; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_test = &FieldTest; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_ms = &FieldMS; //pointer to intermediate Field solutions
  vector<array<double,4>>* field_ms_source = &FieldMS_Source; //pointer to intermediate Field solutions
  vector<array<double,3>>* exact_sols = &ExactField; //pointer to exact solution field values
  vector<array<double,3>>* exact_faces = &ExactFaces; //pointer to exact solution field values
  vector<array<double,3>>* resid = &Residual; //pointer to residual field values per cell
  vector<array<double,3>>* resid_star = &ResidualStar; //pointer to intermediate residual field values per cell
  vector<array<double,3>>* init_resid = &InitResidual; //pointer to residual field values per cell
  vector<double>* time_steps = &TimeSteps;

  //!OBJECT INITIALIZATIONS

  SpaceVariables1D Sols; //for operating on Field variables

  SpaceVariables2D Sols_TEST; //for operating on Field variables

  Tools tool; //utilities object

  Euler1D Euler(cellnum,stag_pressure,back_pressure,stag_temp,gamma); //for performing Euler Eq. operations 

  EulerExplicit Time(cellnum); //for computing time steps

  Output Error; //for discretization error operations

  //Pointers to Objects

  SpaceVariables1D* sols = NULL;
  SpaceVariables2D* sols_test = &Sols_TEST;
  Euler1D* euler = &Euler;

  EulerBASE* euler_test;
  //Temp -- will add scenario == 1 once 1D section is fixed!
  if (scenario == 2) 
    euler_test = new Euler2D(case_2d,mesh->Nx-1,mesh->Ny-1,flux_scheme1,flux_accuracy1,top_cond,btm_cond,left_cond,right_cond,mesh,field_ms_source);
  else if (scenario == 3)
    euler_test = new Euler2DMMS(mesh->Nx-1,mesh->Ny-1,flux_scheme1,flux_accuracy1,top_cond,btm_cond,left_cond,right_cond,mesh,field_ms_source);
  else{
    cerr<<"Error: scenario # not recognized!"<<endl;
    return 0;
  }
  



  //if (mesh_2d) //reassigns pointer to newly created 2D SpaceVariable
    //sols = new SpaceVariables2D();

  /*if (mesh_2d){ //reassigns pointer to newly created 2D Objects
    euler = new Euler2D();
    sols = new SpaceVariables2D();
  }
  else
    sols = new SpaceVariables1D();
  */

  EulerExplicit* time = &Time;
  Output* error = &Error;

  //Intermediate variables for under-relaxation
  array<double,3> ResidualStarNorms; //stores the intermediate global residual norms
  array<bool,3> check{false,false,false}; //false by default to check if under-relaxation is needed


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
  Tools::print("--Cell Number: %d\n",cellnum);
  Tools::print("--Delta x: %f\n",dx);
  // Temporal Stats
  Tools::print("-Temporal Statistics:\n");
  Tools::print("--CFL: %f\n",CFL);
  Tools::print("--Time-stepping method: ");
  (timestep == true) ? Tools::print("Local time-stepping\n") : Tools::print("Global time-stepping\n");  
  // Temporal Stats
  Tools::print("-Flux Statistics:");
  if (flux_scheme == false){ //Upwind Schemes
    if (flux_accuracy == true)//1st order
      (upwind_scheme == true) ? Tools::print(" 1st Order Van Leer Scheme\n") : Tools::print(" 1st Order Roe's Scheme\n");
    else if (flux_accuracy == false)//2nd order
      (upwind_scheme == true) ? Tools::print(" 2nd Order Van Leer Scheme\n") : Tools::print(" 2nd Order Roe's Scheme\n");
  }
  else
    Tools::print(" JST Damping\n");

  //debug: for visualizing manufactured sol.
  /* 
  vector<array<double,4>> FieldTest(cellnum); 
  vector<array<double,4>>* field_test = &FieldTest;
  Euler2D Erick; SpaceVariables2D Emma;
  string file = "2DSols.dat";
  Erick.ManufacturedPrimitiveSols(field_test,mesh_2d->imax,mesh_2d->jmax,xcoords,ycoords,Emma,cellnum);
  Emma.AllOutputPrimitiveVariables(field_test,file,false,0,xcoords,ycoords,cellnum,mesh_2d->imax,mesh_2d->jmax);
  return 0;
  */ 

  //! SETTING INITIAL CONDITIONS
  //Tools::print("At initial conditions\n");
  euler_test->SetInitialConditions(field_test);

  //! COMPUTING EXACT SOLUTION -- ONLY FOR 1D QUASI-STEADY NOZZLE
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

  //Debug: Temporarily set initial conditions to exact solutions
  //Field = ExactField;
  //Debug: printing initial conditions w/ no BCs
  const char* filename = "InitialSolutions.txt";
  sols->OutputPrimitiveVariables(field,euler,filename);

  // SETTING BOUNDARY CONDITIONS
  /*
  mesh->GenerateGhostCells(1,1,1,1);
  const char* right_ghost_coords = "RightGhostCoords.dat";
  const char* left_ghost_coords = "LeftGhostCoords.dat";
  const char* top_ghost_coords = "TopGhostCoords.dat";
  const char* btm_ghost_coords = "BtmGhostCoords.dat";
  error->OutputGhostCells(right_ghost_coords,mesh->right_xcoords,mesh->right_ycoords,mesh->Ny,2);
  error->OutputGhostCells(left_ghost_coords,mesh->left_xcoords,mesh->left_ycoords,mesh->Ny,2);
  error->OutputGhostCells(top_ghost_coords,mesh->top_xcoords,mesh->top_ycoords,mesh->Nx,2);
  error->OutputGhostCells(btm_ghost_coords,mesh->btm_xcoords,mesh->btm_ycoords,mesh->Nx,2);
  //euler_test->Setup2DBoundaryConditions(field_test); //TODO: Get rid of this & use GenerateGhostCells + Enforce2DBoundaryConditions
  */
  //debug: MMS enforcing boundary conditions
  //TODO: make a fcn. to visualize all the ghost cells (not just the last layer)
  const char* filename_ghost = "Ghost_Cell_Solutions.dat";
  vector<array<double,4>>* btm_ghost_cells = &mesh->btm_cells; //creating a pointer to use in output primitive variables fcn.
  euler_test->Setup2DBoundaryConditions(field_test,error);
  error->OutputPrimitiveVariables(btm_ghost_cells,filename_ghost,false,0,mesh->btm_xcoords,mesh->btm_ycoords,(int)mesh->btm_cells.size(),mesh->cell_imax+1,2);
  delete mesh; delete euler_test;
  return 0;

  time->SolutionLimiter(field);

  //! COMPUTING MANUFACTURED SOLUTION AND SOURCE TERMS (MMS ONLY)
  if (scenario == 3){
    string mms_sol_filename = "ManufacturedSols.dat";
    string mms_source_filename = "SourceTerms.dat";
    euler_test->ManufacturedPrimitiveSols(field_ms,sols_test); //!< computing manufactured sol.
    euler_test->EvalSourceTerms(sols_test); //!< computing manufactured source terms
    error->OutputPrimitiveVariables(field_ms,mms_sol_filename,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);
    error->OutputManufacturedSourceTerms(field_ms_source,mms_source_filename,false,0,mesh->xcoords,mesh->ycoords,mesh->cellnumber,mesh->Nx,mesh->Ny);

    delete euler_test;
    delete mesh;
  }
  return 0;

  //!< Outputting initial solutions with BC's
  const char* filename2 = "InitSolutionswBCs.txt";
  sols->OutputPrimitiveVariables(field,euler,filename2);

  // COMPUTING INITIAL RESIDUAL NORMS
  // using ResidSols spacevariable
  array<double,3> InitNorms;

  //if (upwind_scheme == false && flux_accuracy == false) //reverting to 1st order if 2nd order Roe is selected
    //flux_accuracy = true;

  euler->ComputeResidual(init_resid,field,field_stall,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme,epsilon,resid_stall); //computing upwind flux 1st order initially
  //euler_test->ComputeResidual(init_resid);
  InitNorms = sols->ComputeSolutionNorms(init_resid); //computing L2 norm of residuals
  Tools::print("-Initial Residual Norms\n");
  Tools::print("--Continuity:%e\n",InitNorms[0]);
  Tools::print("--X-Momentum:%e\n",InitNorms[1]);
  Tools::print("--Energy:%e\n\n",InitNorms[2]);


  (*resid) = (*init_resid);//!< setting initial residual to intermediate
  ResidualNorms = InitNorms;

  string it,name; //used for outputting file name
  int iter; //iteration number

  //Opening file that stores residuals
  ofstream myresids;
  myresids.open("SolResids.dat");
  myresids<<"variables= \"Iteration num.\" \"Continuity\" \"Momentum\"  \"Energy\""<<endl;
  myresids<<"zone T= "<<"\""<<0<<"\""<<endl;
  myresids<<"DATAPACKING=POINT"<<endl;
  myresids<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  //myresids<<"Iteration"<<"  "<<"Contintuity"<<"  "<<"X-Momentum"<<"  "<<"Energy"<<endl;

  myresids<<0<<"  "<<InitNorms[0]<<"  "<<InitNorms[1]<<"  "<<InitNorms[2]<<endl; //printing out the initial residuals first

  //Printing to TECPLOT
  std::string filename_totalsols = "AllSolutions.dat";
  sols->AllOutputPrimitiveVariables(field,euler,filename_totalsols,false,0,xcoords);

  //Assigning Intermediate Field to Initial Field (including residuals)
  (*field_star) = (*field);
  (*resid_star) = (*resid);

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
    (*time_steps) = (timestep == true) ? time->ComputeLocalTimeStep(field,euler,CFL,dx) : time->ComputeGlobalTimeStep(field,euler,CFL,dx);

    //! COMPUTE NEW SOL. VALUES 
    time->FWDEulerAdvance(field_star,resid_star,euler,time_steps,Omega,mesh);//TESTING
    time->SolutionLimiter(field_star); //applies solution limiter to all cells (including ghost cells)

    //! COMPUTE BOUNDARY CONDITIONS
    euler->ComputeTotalBoundaryConditions(field_star,cond_bc);
    time->SolutionLimiter(field_star); //temporarily reapplying the limiter


    euler->ComputeResidual(resid_star,field_star,field_stall,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme,epsilon,resid_stall);
    ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
    time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);


    //! UNDER-RELAXATION CHECK
    if (check[0]==true || check[1] == true || check[2] == true){ //perform under-relaxation if any of these are true
      for (int j=0;j<subiter_max;j++){
        for (int i=0;i<3;i++) //!< reassigns omega to half of current value if under-relaxation detected
          Omega[i] = (check[i] == true) ?  Omega[i] /= 2.0 : Omega[i] = 1.0;

        (*field_star) = (*field); //resetting primitive variables to previous time step values
        time->FWDEulerAdvance(field_star,resid_star,euler,time_steps,Omega,mesh); //advancing intermediate solution w/ under-relaxation factor 
        euler->ComputeResidual(resid_star,field_star,field_stall,xcoords,dx,flux_scheme,flux_accuracy,upwind_scheme,epsilon,resid_stall); //compute under-relaxed residual
        ResidualStarNorms = sols->ComputeSolutionNorms(resid_star);
        time->UnderRelaxationCheck(ResidualNorms,ResidualStarNorms,C,check);

        if (check[0]==false && check[1] == false && check[2] == false){ //checking if new residuals now do not need under-relaxation
        for (int i=0;i<3;i++) //!< resetting omega to 1
          Omega[i] = 1.0;
        break;
        }

        if (j == subiter_max)
          Tools::print("Under-relaxation Failed!\n");
      }
    }

    // Checking if Residuals are stalled
    if (resid_stall == false){ //only checking if haven't already been marked as stalled
      Prev_ResidualNorms = ResidualNorms;
      resid_stall = time->CheckStallResids(stall_count,ResidualNorms,Prev_ResidualNorms,sols);
      if (resid_stall == true)
        FieldStall = Field; //setting previous field sols. before stall
    }


    

    //Assinging New Time Step Values to Intermediate Values
    (*field) = (*field_star);
    (*resid) = (*resid_star); ResidualNorms = ResidualStarNorms;

    //Computing Ramping Value
    //epsilon = sols->ComputeRampValue(ResidualNorms,InitNorms,ramp_stop);


    //! OUTPUT SOL. IN TEXT FILE EVERY "ITEROUT" STEPS
    if (iter % iterout == 0) {
      it = to_string(iter);
      name = "SolResults/Iteration";
      name += it;
      name += ".txt";
      const char* filename_iter = name.c_str();
      sols->OutputPrimitiveVariables(field,euler,filename_iter);

      //Printing to TECPLOT
      Sols.AllOutputPrimitiveVariables(field,euler,filename_totalsols,true,iter,xcoords);

      //Printing Residual Norms to Screen
      Tools::print("------Iteration #: %d----------\n",iter);
      Tools::print("Continuity:%e\nX-Momentum:%e\nEnergy:%e\n",ResidualNorms[0],ResidualNorms[1],ResidualNorms[2]);
      //debug:
      Tools::print("Epsilon: %e\n",epsilon);
      if (resid_stall == true)
        Tools::print("Residuals are detected to be stalled!\n"); //printing message is temp. for now

      // Writing Residuals history to "SolResids.txt" file
      myresids<<iter<<"  "<<ResidualNorms[0]<<"  "<<ResidualNorms[1]<<"  "<<ResidualNorms[2]<<endl;

    }


    //! CHECK FOR CONVERGENCE (w/ respect to the intial residual norms)
    if (ResidualNorms[0]/InitNorms[0] <= cont_tol && ResidualNorms[1]/InitNorms[1] <= xmom_tol && ResidualNorms[2]/InitNorms[2] <= energy_tol)
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
    const char* filename_final = "ConvergedSolution.txt" ;
    sols->OutputPrimitiveVariables(field,euler,filename_final);
  }

  //Closing Residuals file
  myresids.close();


  //! EVALUATE DISCRETIZATION NORMS FOR GRID CONVERGENCE AND PRINT OUT TO FILE
  if (cond_bc == false){
    field->erase(field->begin()); field->erase(field->begin()); //!< erasing ghost cells
    field->erase(field->end()); field->erase(field->end());

    vector<array<double,3>> Errors(Field);
    vector<array<double,3>>* errors = &Errors;

    error->DiscretizationErrorNorms(field,exact_sols,errors,sols);
  }

  stop_time = MPI_Wtime();
  Tools::print("Elapsed time: %fs\n",stop_time-start_time);

  //! CLEANUP
  if (meshfile){
    delete euler;
  }

  return 0;
}
