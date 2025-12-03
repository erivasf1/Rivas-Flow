//User-defined functions
#include "TimeIntegrator.h" 

// EULEREXPLICITBASE DEFINITIONS

//-----------------------------------------------------------
EulerExplicitBASE::EulerExplicitBASE(MeshGenBASE* &m,EulerBASE* &e,const double &c)
: mesh(m),euler(e) ,CFL(c)
{}

//-----------------------------------------------------------
vector<double> EulerExplicitBASE::ComputeLocalTimeStep(vector<array<double,4>>* &field){

  //refer to Sec7-2D-Slide:37
  array<double,2> lambda_max; //speed of sound + velocity
  vector<double> time_steps(mesh->cellnumber);
  double area_iavg, area_javg;

  double vol;
  int index;


  for (int j=0;j<mesh->cell_jmax;j++){ //looping through all interior cells
    for (int i=0;i<mesh->cell_imax;i++){ //looping through all interior cells

    //retrieving cell volume and interior cell areas
    vol = mesh->GetCellVolume(i,j);
    area_iavg = mesh->GetInteriorCellArea(i,j,0);
    area_javg = mesh->GetInteriorCellArea(i,j,1); 

    index = i + (j*mesh->cell_imax);

    //retrieving outward normal unit vectors
    array<double,2> n_iplushalf = mesh->ComputeOutwardUnitVector(i,j,3);
    array<double,2> n_iminushalf = mesh->ComputeOutwardUnitVector(i,j,2);
    array<double,2> n_jplushalf = mesh->ComputeOutwardUnitVector(i,j,0);
    array<double,2> n_jminushalf = mesh->ComputeOutwardUnitVector(i,j,1);

    //averaging outward normal unit vectors to get NORMAL UNIT VECTOR
    //NOTE: vectors are subtracted because avg. unit vector is not outward!
    array<double,2> n_iavg, n_javg;
    for (int n=0;n<2;n++){
      n_iavg[n] = (n_iplushalf[n] - n_iminushalf[n]) / 2.0;
      n_javg[n] = (n_jplushalf[n] - n_jminushalf[n]) / 2.0;
    }

    //computing lambda max for both in i and j dir.
    lambda_max = euler->GetLambdaMax(field,i,j);

    //error handling
    if (std::isnan(lambda_max[0]) || isnan(lambda_max[1]) ){
      Tools::print("Infinitiy detected!\n");
      Tools::print("Velocity at loc(%d): %f\n",i,(*field)[i][0]);
      //Tools::print("Mach # at loc(%d): %f\n",i,euler->ComputeMachNumber(field,i));
    }

    //computing time step!
    time_steps[index] = CFL * vol / 
                        ( (lambda_max[0]*area_iavg)+(lambda_max[1]*area_javg) );
     
    }
  }

  return time_steps;
}

//-----------------------------------------------------------
vector<double> EulerExplicitBASE::ComputeGlobalTimeStep(vector<array<double,4>>* &field){

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
void EulerExplicitBASE::FWDEulerAdvance(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega){

  double vol;
  array<double,4> conserve;
  array<double,4> primitive;
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

      primitive = euler->ComputePrimitive(conserve); //!< extracting new primitive variables & reassigning field w/ new values

      fieldij(field,i,j,mesh->cell_imax) = primitive;

    }
  }
  
  return;

}

//-----------------------------------------------------------
void EulerExplicitBASE::SolutionLimiter(vector<array<double,4>>* &field){

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
void EulerExplicitBASE::UnderRelaxationCheck(array<double,4> ResidPrevNorm,array<double,4> ResidNorm,double C,array<bool,4> &check){

  for (int i=0;i<4;i++){
    if (ResidNorm[i] > C*ResidPrevNorm[i])
      check[i] = true; //assigning true to corresponding equation
  }


  return;
}

//-----------------------------------------------------------
bool EulerExplicitBASE::CheckStallResids(int &count,array<double,4> &ResidNorms,array<double,4> &Prev_ResidualNorms,SpaceVariables2D* &sol){

  int count_tol = 1e5; double diff_tol = 1.0e-2; 
  double resid_avg = sol->ComputeNormAvg(ResidNorms); 
  double prev_resid_avg = sol->ComputeNormAvg(Prev_ResidualNorms); 

  count = (abs(resid_avg-prev_resid_avg) <= diff_tol) ? count + 1 : count;

  bool stall = (count > count_tol) ? true : false;

  return stall;

}

//-----------------------------------------------------------
void EulerExplicitBASE::ComputeNewSolution(vector<array<double,4>>* &,vector<array<double,4>>* &,vector<double>* &,array<double,4> &,vector<array<double,4>>* &,bool &) {
  return;
}
//-----------------------------------------------------------
EulerExplicitBASE::~EulerExplicitBASE(){}

// EULEREXPLICIT DEFINITIONS
//-----------------------------------------------------------
EulerExplicit::EulerExplicit(MeshGenBASE* &m,EulerBASE* &e,const double &c) : EulerExplicitBASE(m,e,c) {}
//-----------------------------------------------------------
void EulerExplicit::ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &,bool &){

  FWDEulerAdvance(field,resid,time_steps,Omega);

  return;

}
//-----------------------------------------------------------
EulerExplicit::~EulerExplicit(){}
//-----------------------------------------------------------

// RUNGEKUTTA2 DEFINITIONS
//-----------------------------------------------------------
RungeKutta2::RungeKutta2(MeshGenBASE* &m,EulerBASE* &e,const double &c) : EulerExplicit(m,e,c) {

  Field_intermediate_cons.resize(mesh->cellnumber); 
  Field_intermediate_prim.resize(mesh->cellnumber); 
  Resid_intermediate.resize(mesh->cellnumber); 

  field_interm_cons = &Field_intermediate_cons; //assigning pointers
  field_interm_prim = &Field_intermediate_prim; 
  resid_interm = &Resid_intermediate;
}

//-----------------------------------------------------------
void RungeKutta2::ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall) {

  //Reference: Class notes section 8, page 5
  
  double dt_over_V;
  array<double,2> alpha{alpha1,alpha2}; 

  array<double,4> conserve;
  array<double,4> primitive;
  int index;

  //! Converting to conserved variabls
  for (int j=0;j<mesh->cell_jmax;j++){
    for (int i=0;i<mesh->cell_imax;i++){

      conserve = euler->ComputeConserved(field,i,j);
      fieldij(field_interm_cons,i,j,mesh->cell_imax) = conserve;

    }
  }

  //! Intermediate stage
  (*resid_interm) = (*resid); //setting resid. interm. to original resid.
  for (int n=0;n<(int)alpha.size();n++){

    double alpha_n = alpha[n];

    // compute U(n) vals.
    for (int j=0;j<mesh->cell_jmax;j++){
      for (int i=0;i<mesh->cell_imax;i++){
        index = i + (j*mesh->cell_imax);
        dt_over_V = (*time_steps)[index] / mesh->GetCellVolume(i,j);
        for (int q=0;q<4;q++){
          (*field_interm_cons)[index][q] = (*field_interm_cons)[index][q] - (alpha_n * dt_over_V * (*resid_interm)[index][q]);
        }
      }
    }

    if (n < (int)alpha.size()-1){ //only evaluating interm. resid. once!
      // eval resid. of U(n) -- need to convert back to primitive to eval. resid.
      for (int j=0;j<mesh->cell_jmax;j++){
        for (int i=0;i<mesh->cell_imax;i++){
          conserve = fieldij(field_interm_cons,i,j,mesh->cell_imax);
          primitive = euler->ComputePrimitive(conserve);
          fieldij(field_interm_prim,i,j,mesh->cell_imax) = primitive;
        }
      }
      euler->ComputeResidual(resid_interm,field_interm_prim,field_stall,resid_stall);
    }


  }



  //! Updating old time step vals. to new time step
  for (int j=0;j<mesh->cell_jmax;j++){
    for (int i=0;i<mesh->cell_imax;i++){
      conserve = fieldij(field_interm_cons,i,j,mesh->cell_imax);
      primitive = euler->ComputePrimitive(conserve);
      fieldij(field,i,j,mesh->cell_imax) = primitive;
    }
  }

  return;

}
//-----------------------------------------------------------
RungeKutta2::~RungeKutta2(){}
//-----------------------------------------------------------

// RUNGEKUTTA4 DEFINITIONS
//-----------------------------------------------------------
RungeKutta4::RungeKutta4(MeshGenBASE* &m,EulerBASE* &e,const double &c) : EulerExplicit(m,e,c) {

  Field_intermediate_cons.resize(mesh->cellnumber); 
  Field_intermediate_prim.resize(mesh->cellnumber); 
  Resid_intermediate.resize(mesh->cellnumber); 

  field_interm_cons = &Field_intermediate_cons; //assigning pointers
  field_interm_prim = &Field_intermediate_prim; 
  resid_interm = &Resid_intermediate;
}

//-----------------------------------------------------------
void RungeKutta4::ComputeNewSolution(vector<array<double,4>>* &field,vector<array<double,4>>* &resid,vector<double>* &time_steps,array<double,4> &Omega,vector<array<double,4>>* &field_stall,bool &resid_stall) {

  //Reference: Class notes section 8, page 5
  
  double dt_over_V;
  array<double,4> alpha{alpha1,alpha2,alpha3,alpha4}; 

  array<double,4> conserve;
  array<double,4> primitive;
  int index;

  //! Converting to conserved variables
  for (int j=0;j<mesh->cell_jmax;j++){
    for (int i=0;i<mesh->cell_imax;i++){

      conserve = euler->ComputeConserved(field,i,j);
      fieldij(field_interm_cons,i,j,mesh->cell_imax) = conserve;

    }
  }

  //! Intermediate stage
  (*resid_interm) = (*resid); //setting resid. interm. to original resid.
  for (int n=0;n<(int)alpha.size();n++){

    double alpha_n = alpha[n];

    // compute U(n) vals.
    for (int j=0;j<mesh->cell_jmax;j++){
      for (int i=0;i<mesh->cell_imax;i++){
        index = i + (j*mesh->cell_imax);
        dt_over_V = (*time_steps)[index] / mesh->GetCellVolume(i,j);
        for (int q=0;q<4;q++){
          (*field_interm_cons)[index][q] = (*field_interm_cons)[index][q] - (alpha_n * dt_over_V * (*resid_interm)[index][q]);
        }
      }
    }

    if (n < (int)alpha.size()-1){ //only evaluating interm. resid. once!
      // eval resid. of U(n) -- need to convert back to primitive to eval. resid.
      for (int j=0;j<mesh->cell_jmax;j++){
        for (int i=0;i<mesh->cell_imax;i++){
          conserve = fieldij(field_interm_cons,i,j,mesh->cell_imax);
          primitive = euler->ComputePrimitive(conserve);
          fieldij(field_interm_prim,i,j,mesh->cell_imax) = primitive;
        }
      }
      euler->ComputeResidual(resid_interm,field_interm_prim,field_stall,resid_stall);
    }


  }



  //! Updating old time step vals. to new time step
  for (int j=0;j<mesh->cell_jmax;j++){
    for (int i=0;i<mesh->cell_imax;i++){
      conserve = fieldij(field_interm_cons,i,j,mesh->cell_imax);
      primitive = euler->ComputePrimitive(conserve);
      fieldij(field,i,j,mesh->cell_imax) = primitive;
    }
  }

  return;

}
//-----------------------------------------------------------
RungeKutta4::~RungeKutta4(){}
//-----------------------------------------------------------
