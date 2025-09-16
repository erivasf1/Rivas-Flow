//User-defined functions
#include "TimeIntegrator.h" 

// EULEREXPLICIT DEFINITIONS

//-----------------------------------------------------------
EulerExplicit::EulerExplicit(MeshGenBASE* &m,EulerBASE* &e,const double &c)
: mesh(m),euler(e) ,CFL(c)
{}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeLocalTimeStep(vector<array<double,4>>* &field){

  //refer to Sec7-2D-Slide:36
  array<double,2> lambda_max; //speed of sound + velocity
  vector<double> time_steps(mesh->cellnumber);
  double area_internal_i,area_internal_j,vol;
  int index;

  for (int j=0;j<mesh->cell_jmax;j++){ //looping through all interior cells
    for (int i=0;i<mesh->cell_imax;i++){ //looping through all interior cells

    //retrieving cell volume and interior cell areas
    area_internal_i = mesh->GetInteriorCellFaceArea(i,j,0); 
    area_internal_j = mesh->GetInteriorCellFaceArea(i,j,1); 
    vol = mesh->GetCellVolume(i,j);

    index = i + (j*mesh->cell_imax);
    lambda_max = euler->GetLambdaMax(field,i,j); //note: [0]:i & [1]:j

    //error handling
    if (std::isnan(lambda_max[0]) || isnan(lambda_max[1]) ){
      Tools::print("Infinitiy detected!\n");
      Tools::print("Velocity at loc(%d): %f\n",i,(*field)[i][0]);
      //Tools::print("Mach # at loc(%d): %f\n",i,euler->ComputeMachNumber(field,i));
    }

    time_steps[index] = CFL * vol / 
                 ( (lambda_max[0]*area_internal_i)+(lambda_max[1]*area_internal_j) );
    //Tools::print("CFL:%f,dx: %f,lambda_max:%f\n",CFL,dx,lambda_max);
     
    }
  }

  return time_steps;
}

//-----------------------------------------------------------
vector<double> EulerExplicit::ComputeGlobalTimeStep(vector<array<double,4>>* &field){

  //extracting smallest local time step of all cells
  vector<double> time_steps = ComputeLocalTimeStep(field); //local time steps list for all cells
  double min_time_step = 1.0e5; //temp. value for min time_step

  for (int n=0;n<(int)time_steps.size();n++){
    if (time_steps[n] < min_time_step)
      min_time_step = time_steps[n];
  }

  vector<double> global_time_steps((int)time_steps.size(),min_time_step);

  return global_time_steps;

}
//-----------------------------------------------------------
void EulerExplicit::FWDEulerAdvance(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega){

  double vol;
  array<double,4> conserve;
  int cell_index;
  //use indexing of interior cells for Resid!
  // Field still has ghost cells
  //First, convert to conservative to compute conservative values at new time step
  //Second, extract primitive variables from newly calculated conservative variables
  for (int j=0;j<mesh->cell_jmax;j++){ 
    for (int i=0;i<mesh->cell_imax;i++){ 

      cell_index = i + (j*mesh->cell_imax);
      vol = mesh->GetCellVolume(i,j); //acquiring cell vol
      conserve = euler->ComputeConserved(field,i,j); //!< computing conservative values

      for (int n=0;n<4;n++) // advancing to new timestep of conservative variable
        conserve[n] -= Omega[n]*((*time_steps)[n] / vol) * (*resid)[cell_index][n];

      euler->ComputePrimitive(field,conserve,i,j); //!< extracting new primitive variables & reassigning field w/ new values

    }
  }
  
  return;

}

//-----------------------------------------------------------
void EulerExplicit::SolutionLimiter(vector<array<double,4>>* &field){

  for (int n=0;n<(int)field->size();n++){
    //Density
    (*field)[n][0] = std::min(Density_max,std::max(Density_min,(*field)[n][0]));

    //X-Velocity
    (*field)[n][1] = std::min(Velocity_max,std::max(Velocity_min,(*field)[n][1]));

    //Y-Velocity
    (*field)[n][2] = std::min(Velocity_max,std::max(Velocity_min,(*field)[n][2]));

    //Pressure
    (*field)[n][3] = std::min(Pressure_max,std::max(Pressure_min,(*field)[n][3]));

    //Printing out message if limiter kicks in
    
    if ((*field)[n][0] == Density_max || (*field)[n][0] == Density_min)
      Tools::print("Limiter was hit for density at cell %d | val is now:%e\n",n,(*field)[n][0]);

    if ((*field)[n][1] == Velocity_max || (*field)[n][1] == Velocity_min)
      Tools::print("Limiter was hit for x-velocity at cell %d | val is now:%e\n",n,(*field)[n][1]);

    if ((*field)[n][2] == Velocity_max || (*field)[n][2] == Velocity_min)
      Tools::print("Limiter was hit for y-velocity at cell %d | val is now:%e\n",n,(*field)[n][2]);

    if ((*field)[n][3] == Pressure_max || (*field)[n][3] == Pressure_min)
      Tools::print("Limiter was hit for pressure at cell %d | val is now:%e\n",n,(*field)[n][3]);
   
  }

  return;

}

//-----------------------------------------------------------
void EulerExplicit::UnderRelaxationCheck(array<double,4> ResidPrevNorm,array<double,4> ResidNorm,double C,array<bool,4> &check){

  for (int i=0;i<4;i++){
    if (ResidNorm[i] > C*ResidPrevNorm[i])
      check[i] = true; //assigning true to corresponding equation
  }


  return;
}

//-----------------------------------------------------------
bool EulerExplicit::CheckStallResids(int &count,array<double,4> &ResidNorms,array<double,4> &Prev_ResidualNorms,SpaceVariables2D* &sol){

  int count_tol = 1e5; double diff_tol = 1.0e-2; 
  double resid_avg = sol->ComputeNormAvg(ResidNorms); 
  double prev_resid_avg = sol->ComputeNormAvg(Prev_ResidualNorms); 

  count = (abs(resid_avg-prev_resid_avg) <= diff_tol) ? count + 1 : count;

  bool stall = (count > count_tol) ? true : false;

  return stall;

}

//-----------------------------------------------------------
EulerExplicit::~EulerExplicit(){}
