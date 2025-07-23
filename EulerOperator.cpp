//User-defined functions
#include "EulerOperator.h" 
#include "MeshAccess.hpp"

// EULERBASE DEFINITIONS

//-----------------------------------------------------------
EulerBASE::EulerBASE(int &cell_inum,int &cell_jnum,int &scheme,int &accuracy,int &top,int &btm,int &left,int &right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source)
  : cell_imax(cell_inum), cell_jmax(cell_jnum),flux_scheme(scheme),flux_accuracy(accuracy),top_cond(top),btm_cond(btm),left_cond(left), right_cond(right), mesh(mesh_ptr), mms_source(source) {}

//-----------------------------------------------------------
array<double,4> EulerBASE::ComputeConserved(vector<array<double,4>>* &field,int &i,int &j){

  array<double,4> res; //to store conserved variables

  //compute cell id
  int cell_id = i + (j*cell_imax);

  //density
  res[0] = (*field)[cell_id][0]; 

  //x-mom.
  res[1] = (*field)[cell_id][0] * (*field)[cell_id][1]; 

  //y-mom.
  res[2] = (*field)[cell_id][0] * (*field)[cell_id][2]; 

  //energy
  res[3] = (*field)[cell_id][2]/(gamma-1.0) + (0.5)* (*field)[cell_id][0]*pow((*field)[cell_id][1],2.0); //TODO: Fix this!

  return res;
  
}
//-----------------------------------------------------------
void EulerBASE::SetInitialConditions([[maybe_unused]] vector<array<double,4>>* &field){
  return;
}

//-----------------------------------------------------------
void EulerBASE::EvalSourceTerms(SpaceVariables2D* &){
//void EulerBASE::EvalSourceTerms(vector<array<double,4>>* &MMS_Source,SpaceVariables2D* &sols,MeshGenBASE* &mesh){
  return;
}

//-----------------------------------------------------------
void EulerBASE::ManufacturedPrimitiveSols(vector<array<double,4>>* &,SpaceVariables2D* &){
//void EulerBASE::ManufacturedPrimitiveSols(vector<array<double,4>>* &field,SpaceVariables2D* &sols,MeshGenBASE* &mesh){
  return;
}

//-----------------------------------------------------------
void EulerBASE::SetupBoundaryConditions(){
  //TODO: Should only include the GenerateGhostCells fcn. of MeshGen
  // & then initialize the ghost cells with a specified val.
  //cond: 0=Inflow, 1 = Outflow, 2 = Slip Wall
  //id: 0 = reflect, 1 = extend
  int btm_id = (btm_cond==0 || btm_cond==1) ? 1:0;
  int top_id = (top_cond==0 || top_cond==1) ? 1:0;
  int left_id = (left_cond==0 || left_cond==1) ? 1:0;
  int right_id = (right_cond==0 || right_cond==1) ? 1:0;
  
  mesh->GenerateGhostCells(left_id,right_id,btm_id,top_id);

  //initializing ghost cells
/*
  else{
    cerr<<"Unknown Boundary Condition applied to Right Section!"<<endl;
    return; 
  }
*/
  
  return;

}
//-----------------------------------------------------------
void EulerBASE::ComputeResidual(vector<array<double,4>>* &,vector<array<double,4>>* &){
  return;
}

//-----------------------------------------------------------
EulerBASE::~EulerBASE(){}
//-----------------------------------------------------------
// EULER1D DEFINITIONS

Euler1D::Euler1D(int& cellnum,double &P0,double &BP,double &T0,double &g)
  : interior_cellnum(cellnum), stag_pressure(P0) , back_pressure(BP), stag_temperature(T0) 
  , gamma(g) {}


//-----------------------------------------------------------
Euler1D::Euler1D() {}


//-----------------------------------------------------------
array<double,3> Euler1D::ComputeConserved(vector<array<double,3>>* &field,int loc){

  array<double,3> res; //to store conserved variables

  //density
  res[0] = (*field)[loc][0]; 

  //x-mom.
  res[1] = (*field)[loc][0] * (*field)[loc][1]; 

  //energy
  res[2] = (*field)[loc][2]/(gamma-1.0) + (0.5)* (*field)[loc][0]*pow((*field)[loc][1],2.0);

  return res;

}

//-----------------------------------------------------------
void Euler1D::ComputePrimitive(vector<array<double,3>>* &field,array<double,3> &Conserved,int loc) {

  //Computing Primitive variables given the conserved variables

  //density
  (*field)[loc][0] = Conserved[0];

  //x-velocity
  (*field)[loc][1] = Conserved[1] / Conserved[0];

  //pressure
  (*field)[loc][2] = (gamma-1.0) * (Conserved[2]-0.5*(pow(Conserved[1],2.0)/Conserved[0]));

  return;
}



//-----------------------------------------------------------
void Euler1D::SetInitialConditions(vector<array<double,3>>* &field,vector<double> &xcoords){

  // Mach number is set to vary linearly via M(x) = 9/10(x) + 1
  // And flow quantities are calculated using isentropic conditions
  // ASSUMPTION: Mach number at left face is equal to the cell averaged Mach number of a given cell (may be fine as an initial condition)
  double M,psi,T,a;

  //Tools::print("SetInitialConditions\n");
  for (int i=0;i<interior_cellnum;i++){
   // Tools::print("cell index: %d\n",i);

    M = (9.0/10.0)*xcoords[i] + 1.0; //local Mach number
    psi = 1.0+(gamma-1.0)/2.0 * pow(M,2.0);
    
    //pressure calc.
    (*field)[i][2] = pow(psi,gamma/(gamma-1.0));
    (*field)[i][2] = stag_pressure / (*field)[i][2]; 
    //Tools::print("pressure: %f\n",field[i][2]);
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    (*field)[i][0] = (*field)[i][2] / (R*T); 
    //Tools::print("density: %f\n",field[i][0]);

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    (*field)[i][1] = abs(M*a);
    //Tools::print("velocity: %f\n",field[i][1]);
     
  }    
  
  return;

}
//-----------------------------------------------------------
void Euler1D::SetBoundaryConditions(vector<array<double,3>>* &field,bool &cond){

  //Inserting 4 ghost cells for inflow and outflow locations (2 for each)
  //iterator it = Field.begin();
  array<double,3> empty{0.0,0.0,0.0};
  //Inserting empty solution vectors into inflow and outflow ghost cells
  //TODO: Check if 2 empty lists were added to the start of Field domain instead of overwriting the 1st 2 elements
  //Inflow
  field->insert(field->begin(),empty); //!< temporarily setting ghost cells to empty arrays 
  field->insert(field->begin(),empty); 

  //Outflow
  field->push_back(empty); 
  field->push_back(empty); 

  total_cellnum = (int)field->size(); //!< saving the new total num. of cells (w/ ghost cells)

  //Calculating Boundary Condition values
  ComputeTotalBoundaryConditions(field,cond);

  return;
  
}

//-----------------------------------------------------------
void Euler1D::ComputeTotalBoundaryConditions(vector<array<double,3>>* &field,bool &cond){

  //Combines setting the inflow and outflow boundary conditions
  ComputeInflowBoundaryConditions(field);

  ComputeOutflowBoundaryConditions(field,cond);


  return;
}
//-----------------------------------------------------------
void Euler1D::ComputeInflowBoundaryConditions(vector<array<double,3>>* &field){

  // Domain:
  // [G1,G2,I1,...,Imax,G3,G4] --> READ THIS!!!
  //Inflow -- linear extrapolation of Mach Number (refer to class notes section 3 slide 34)
  // note: use the absolute velocity when computing Mach number to prevent negative velocities at the inflow
  double M0,M1,M2; //temp. values
  double psi;
  double T,a;
  double M_limit = 1.0e-3; //Mach num. limiter
  for (int i=1;i>-1;i--){ //reverse for loop for ease of indexing!
    // i=0 (G1) & i=1 (G2)
    //acquiring Mach number from neighboring nodes
    M1 = GetMachNumber(field,i+1);
    M2 = GetMachNumber(field,i+2);
    M0 = 2.0*M1 - M2;

    //applying limiter
    M0 = (M0<M_limit) ? M_limit : M0;

    //USING ISENTROPIC CONDITIONS for thermodynamic variables

    psi = 1.0+(gamma-1.0)/2.0 * pow(M0,2.0);
    //pressure calc.
    (*field)[i][2] = pow(psi,gamma/(gamma-1.0));
    (*field)[i][2] = stag_pressure / (*field)[i][2]; 
    
    //density calc.
    T = stag_temperature / psi; // local temperature
    (*field)[i][0] = (*field)[i][2] / (R*T); 

    //velocity calc.
    a = sqrt(gamma*R*T); //local speed of sound
    (*field)[i][1] = abs(M0*a);
    

    //Tools::print("B.C.\n");
    //Tools::print("Cell index:%d\n",i);
    //Tools::print("Density: %f, Velocity:%f,Pressure:%f\n",field[i][0],field[i][1],field[i][2]);
   }

  return;
}

//-----------------------------------------------------------
void Euler1D::ComputeOutflowBoundaryConditions(vector<array<double,3>>* &field,bool& cond){

  //TODO: Figure out why [i-1] node is not retrieving the node of interior node

  if (cond == false){ //supersonic case 
    //using simple extrapolation from class notes section 3 slide 36
    int c = (int)field->size() - 2; //index of the 1st outflow ghost cell
    for (int i=0;i<2;i++){ //for both outflow ghost cells
      (*field)[i+c][0] = 2.0* (*field)[(i+c)-1][0] - (*field)[(i+c)-2][0]; //density
      (*field)[i+c][1] = 2.0* (*field)[(i+c)-1][1] - (*field)[(i+c)-2][1]; //velocity
      (*field)[i+c][2] = 2.0* (*field)[(i+c)-1][2] - (*field)[(i+c)-2][2]; //pressure

    }

  }

  else { //subsonic
    //fixing back pressure at boundary, not at ghost cell
    int c = (int)field->size() - 2; //index of the 1st outflow ghost cell
    (*field)[c][2] = 2.0*back_pressure - (*field)[c-1][2]; //extrapolated pressure with back pressure enforced
    (*field)[c+1][2] = 2.0* (*field)[c][2] - (*field)[c-1][2]; //regular extrapolated pressure for last ghost cell

    for (int n=c;n<c+2;n++){ //reg. extrapolation
      for (int i=0;i<2;i++) //only iterating for density & velocity
        (*field)[n][i] = 2.0* (*field)[n-1][i] - (*field)[n-2][i]; //density
    }

  }


  return;
}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_CELL(array<double,3> &field_state){

  array<double,3> flux;

  //Continuity Flux
  flux[0] = field_state[0] * field_state[1]; //rho*u

  //X-Momentum Flux
  flux[1] = field_state[0]*pow(field_state[1],2.0) + field_state[2]; //rho*u^2 + P

  //Energy Flux
  flux[2] =(gamma/(gamma-1.0))*field_state[2]*field_state[1] + (field_state[0]*pow(field_state[1],3.0))*0.5; //(gamma/gamma-1)*P*u + (rho*u^3)/2

  return flux; 

}
//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_BASE(vector<array<double,3>>* &field,int loc,int nbor){

  array<double,3> V_face; //primitive variable vector

  //RECONSTRUCTION: approximating primitive variables at face using central quadrature
  for (int i=0;i<3;i++)
    V_face[i] = 0.5*((*field)[loc][i]+(*field)[nbor][i]);

  array<double,3> flux;

  //Continuity Flux
  flux[0] = V_face[0]*V_face[1]; //rho*u

  //X-Momentum Flux
  flux[1] = V_face[0]*pow(V_face[1],2.0) + V_face[2]; //rho*u^2 + P

  //Energy Flux
  flux[2] =(gamma/(gamma-1.0))*V_face[2]*V_face[1] + (V_face[0]*pow(V_face[1],3.0))*0.5; //(gamma/gamma-1)*P*u + (rho*u^3)/2


  return flux; 

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,3>>* &field,bool method,int loc,int nbor){

  //NOTE: assumes loc is left state and nbor is right state (for 1st order accuracy)
  array<double,3> left_state = (*field)[loc];
  array<double,3> right_state = (*field)[nbor];

  array<double,3> flux; //total flux
  [[maybe_unused]] array<double,3> flux_rtstate; //right state flux (for upwinding in +c wave speed)
  [[maybe_unused]] array<double,3> flux_ltstate; //left state flux (for upwinding in -c wave speed)


  if (method == true){ //Van Leer Method
    flux_rtstate = VanLeerCompute(right_state,false); //false for negative c case
    flux_ltstate = VanLeerCompute(left_state,true); //true for positive c case
   
    for (int n=0;n<3;n++) //summing up left and right state fluxes
      flux[n] = flux_rtstate[n] + flux_ltstate[n];

  }
  else if (method == false) //Roe's Method
    flux = ComputeRoeFlux(left_state,right_state);


  return flux;
}
//-----------------------------------------------------------
array<double,3> Euler1D::ComputeSpatialFlux_UPWIND2ndOrder(vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,bool method,int loc,int nbor,double epsilon,bool &resid_stall){

  //NOTE: using MUSCL extrapolation to compute left and right states
  array<array<double,3>,2> field_states = MUSCLApprox(field,field_stall,loc,nbor,epsilon,resid_stall);
  array<double,3> left_state = field_states[0];
  array<double,3> right_state = field_states[1];
  array<double,3> flux; //total flux
  array<double,3> flux_rtstate; //right state flux (for upwinding in +c wave speed)
  array<double,3> flux_ltstate; //left state flux (for upwinding in -c wave speed)

  if (method == true){ //Van Leer Method
    flux_rtstate = VanLeerCompute(right_state,false); //false for negative c case
    flux_ltstate = VanLeerCompute(left_state,true); //true for positive c case
   
    for (int n=0;n<3;n++) //summing up left and right state fluxes
      flux[n] = flux_rtstate[n] + flux_ltstate[n];

  }
  else if (method == false) //Roe's Method
    flux = ComputeRoeFlux(left_state,right_state);


  return flux;
}


//-----------------------------------------------------------
array<double,3> Euler1D::VanLeerCompute(array<double,3> &field_state,bool sign){

  array<double,3> flux;
  //scalars
  double M = ComputeMachNumber(field_state); //local Mach Number
  double a = field_state[1] / M; //local speed of sound
 
  //total energy term(h_t)
  double ht = (gamma/(gamma-1.0))*(field_state[2]/field_state[0]); //pressure work
  ht += pow(field_state[1],2.0) / 2.0; //kinetic energy

  //vectors
  array<double,3> convect_vec{1.0,field_state[1],ht}; //vector multiple for convective flux
  array<double,3> pressure_vec{0.0,field_state[2],0.0}; //vector multiple for pressure flux

  //Convective Flux
  double C = GetC(M,sign);
  for (int n=0;n<3;n++)
    flux[n] = field_state[0]*a*C*convect_vec[n]; 

  //Pressure Flux
  double D = GetD(M,sign);
  for (int n=0;n<3;n++)
    flux[n] += D*pressure_vec[n];
  

  return flux;
}

//-----------------------------------------------------------
double Euler1D::GetC(double M,bool sign){

  double alpha = GetAlpha(M,sign);
  double beta = GetBeta(M);
  double M_vl = GetVanLeerM(M,sign);

  double C = alpha*(1.0+beta)*M - beta*M_vl;
 
  return C;

}
//-----------------------------------------------------------
double Euler1D::GetAlpha(double M,bool sign){

  //sign = true for positive and false for negative
  double alpha = (sign == true) ? 0.5*(1.0+std::copysign(1.0,M)) : 0.5*(1.0-std::copysign(1.0,M));

  return alpha;

}
//-----------------------------------------------------------
double Euler1D::GetBeta(double M){

  double beta = -std::max(0.0,1.0-(int)abs(M)); //(int) -- explicit type conversion

  return beta;
}
//-----------------------------------------------------------
double Euler1D::GetVanLeerM(double M,bool sign){
  
  double M_VL = (sign == true) ? 0.25*pow((M+1.0),2.0) : -0.25*pow((M-1.0),2.0);
  return M_VL;
}
//-----------------------------------------------------------
double Euler1D::GetD(double M,bool sign){

  double alpha = GetAlpha(M,sign);
  double beta = GetBeta(M);
  double p2bar = GetP2Bar(M,sign);

  double D = alpha*(1.0+beta) - beta*p2bar;

  return D;
}
//-----------------------------------------------------------
double Euler1D::GetP2Bar(double M,bool sign){

  double M_vl = GetVanLeerM(M,sign);

  double p2bar = (sign == true) ? M_vl*(-M+2.0) : M_vl*(-M-2.0);

  return p2bar;
}
//-------------------------------------------------------------------------------------------------------
array<double,3> Euler1D::ComputeRoeFlux(array<double,3> &field_ltstate,array<double,3> &field_rtstate){ 

  //Note: Refer to Lecture Notes Section 6,Page 46
  array<double,3> flux;
  array<double,3> flux_left = ComputeSpatialFlux_CELL(field_ltstate);
  array<double,3> flux_right = ComputeSpatialFlux_CELL(field_rtstate);
  
  double abar; //roe-avg speed of sound -> assigned in subsequent fcn.
  array<double,3> roe_vars = ComputeRoeAvgVars(field_ltstate,field_rtstate,abar); //roe-avg. variables

  array<double,3> roe_eigenvals = ComputeRoeEigenVals(roe_vars,abar); 
  array<array<double,3>,3> roe_eigenvecs = ComputeRoeEigenVecs(roe_vars,abar);
  array<double,3> wave_amps = ComputeRoeWaveAmps(roe_vars,field_ltstate,field_rtstate,abar); 

  //shock-tube problem waves (assuming Riemann problem)
  array<double,3> shock_waves;
  array<array<double,3>,3> j_vec;
  for (int j=0;j<3;j++){ //computing j vectors
    for (int n=0;n<3;n++)
    j_vec[j][n] = abs(roe_eigenvals[j])*wave_amps[j]*roe_eigenvecs[j][n];
  }
  for (int j=0;j<3;j++) //summing up j vectors
    shock_waves[j] = 0.5*(j_vec[0][j]+j_vec[1][j]+j_vec[2][j]);
  
  for (int n=0;n<3;n++) 
    flux[n] = 0.5*(flux_left[n]+flux_right[n]) - shock_waves[n];

  return flux;
}

//---------------------------------------------------------------------------------------------------------------
array<double,3> Euler1D::ComputeRoeWaveAmps(array<double,3> &roe_vars,array<double,3> &field_ltstate,array<double,3> &field_rtstate,double abar){

  array<double,3> wave_amps;
  double delta_rho = field_rtstate[0] - field_ltstate[0];
  double delta_u = field_rtstate[1] - field_ltstate[1];
  double delta_p = field_rtstate[2] - field_ltstate[2];

  wave_amps[0] = delta_rho - (delta_p / pow(abar,2.0));
  wave_amps[1] = delta_u + (delta_p / (roe_vars[0]*abar));
  wave_amps[2] = delta_u - (delta_p / (roe_vars[0]*abar));

  return wave_amps;

}

//-----------------------------------------------------------
array<array<double,3>,3> Euler1D::ComputeRoeEigenVecs(array<double,3> &roe_vars,double abar){
 
  //Refer to notes: Section 6, Page 44

  array<double,3> evec1{1.0,roe_vars[1],pow(roe_vars[1],2.0)*0.5};
  array<double,3> evec2{1.0,roe_vars[1]+abar,roe_vars[2]+(roe_vars[1]*abar)};
  array<double,3> evec3{1.0,roe_vars[1]-abar,roe_vars[2]-(roe_vars[1]*abar)};

  for (int n=0;n<3;n++) //multiplying by scaling factor
    evec2[n] *= roe_vars[0]/(2.0*abar);

  for (int n=0;n<3;n++) //multiplying by scaling factor
    evec3[n] *= -roe_vars[0]/(2.0*abar);

  array<array<double,3>,3> FullEigenVecs{evec1,evec2,evec3};

  return FullEigenVecs;
}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeRoeEigenVals(array<double,3> &roe_vars,double abar){

  array<double,3> roe_eigenvals{roe_vars[1],roe_vars[1]+abar,roe_vars[1]-abar};

  //Modification to eigenvalues by Harten,1983 to handle rarefaction fans by smoothing them
  double epsilon = 0.5;

  for (int i=0;i<3;i++)
    roe_eigenvals[i] = (abs(roe_eigenvals[i]) < 2.0*epsilon*abar) ? (pow(roe_eigenvals[i],2.0)/(4.0*epsilon*abar)) + (epsilon*abar) : abs(roe_eigenvals[i]);

  return roe_eigenvals;

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeRoeAvgVars(array<double,3> &field_ltstate,array<double,3> &field_rtstate,double &abar){

  //preliminary terms
  array<double,3> roe_avg;
  double R_ihalf = sqrt(field_rtstate[0]/field_ltstate[0]);
  //double R_ihalf = sqrt((*field)[nbor][0]/(*field)[loc][0]);

  double ht_left = (gamma/(gamma-1.0))*(field_ltstate[2]/field_ltstate[0]); //pressure work
  ht_left += pow(field_ltstate[1],2.0) / 2.0; //kinetic energy

  double ht_right = (gamma/(gamma-1.0))*(field_rtstate[2]/field_rtstate[0]); //pressure work
  ht_right += pow(field_rtstate[1],2.0) / 2.0; //kinetic energy


  //Rho-avg.
  roe_avg[0] = R_ihalf * field_ltstate[0];
  //roe_avg[0] = R_ihalf * (*field)[loc][0];

  //U-avg.
  roe_avg[1] = (R_ihalf * field_rtstate[1] + field_ltstate[1]) / (R_ihalf + 1.0); 

  //ht-avg.
  roe_avg[2] = (R_ihalf * ht_right + ht_left) / (R_ihalf + 1.0);

  //a-avg
  abar = (gamma-1.0) * (roe_avg[2] - (pow(roe_avg[1],2.0)*0.5));
  abar = sqrt(abar);

  return roe_avg;

}
//-----------------------------------------------------------
array<array<double,3>,2> Euler1D::MUSCLApprox(vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,int loc,int nbor,double epsilon,bool &resid_stall){

  array<double,3> left_state,right_state;
  array<array<double,3>,2> total_state;
  //if (loc == interior_cellnum - 1)
    //Tools::print("I am here\n");  

  //kappa scheme: -1 = fully upwinded; 0 = upwind biased
  // 1/3 = 3rd order upwind scheme;1/2 = QUICK scheme;1 = central scheme
  double kappa = 1.0;

  //computing flux limiters if not frozen
  // Options: Beta & Van Leer
  //NOTE: Refer to Lecture Notes: Section 7 Page 24
  //psi+- of i-1/2, psi-of i+3/2, psi+ of i-1/2
  [[maybe_unused]] double beta = 1.5;
  int r_nbor = nbor+1; int l_nbor = loc - 1;

  //evaluating limiters based on if resid. are stalled or not
  vector<array<double,3>> *field_selected = (resid_stall == false) ? field:field_stall;

  array<double,3> r_centralplus = ComputeRPlusVariation(field_selected,loc,r_nbor,nbor); //i+1/2
  array<double,3> r_centralminus = ComputeRMinusVariation(field_selected,loc,l_nbor,nbor); //i+1/2
  array<double,3> r_upwindminus = ComputeRMinusVariation(field_selected,loc+1,l_nbor+1,nbor+1); //i+3/2
  array<double,3> r_upwindplus = ComputeRPlusVariation(field_selected,loc-1,r_nbor-1,nbor-1); //i-1/2

  //evaluating limiters (psi vectors)
  array<double,3> psi_centralplus = ComputeVanLeerLimiter(r_centralplus); //i+1/2
  array<double,3> psi_centralminus = ComputeVanLeerLimiter(r_centralminus); //i+1/2
  array<double,3> psi_upwind_ltstate = ComputeVanLeerLimiter(r_upwindplus); //i-1/2
  array<double,3> psi_upwind_rtstate = ComputeVanLeerLimiter(r_upwindminus); //i+3/2
  
  //computing left state
  for (int n=0;n<3;n++){
    left_state[n] = (1.0-kappa)*psi_upwind_ltstate[n]*((*field)[loc][n]-(*field)[loc-1][n]) + (1.0+kappa)*psi_centralminus[n]*((*field)[nbor][n] - (*field)[loc][n]);
    left_state[n] = (*field)[loc][n] + ((epsilon/4.0)*left_state[n]);
  }
  //computing right state
  for (int n=0;n<3;n++){
    right_state[n] = (1.0-kappa)*psi_upwind_rtstate[n]*((*field)[nbor+1][n]-(*field)[nbor][n]) + (1.0+kappa)*psi_centralplus[n]*((*field)[nbor][n] - (*field)[loc][n]);
    right_state[n] = (*field)[nbor][n] - ((epsilon/4.0)*right_state[n]);
  }

  //combining left and right state
  total_state[0] = left_state;
  total_state[1] = right_state;

  //TODO: Check left and right state vals. in GDB
  return total_state;
}

//-----------------------------------------------------------
array<array<double,3>,2> Euler1D::ComputeBetaLimiter(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor,double beta){

  //Note: This is really a state limiter!
  //Using Beta Limiter and VanLeer Limiter

  array<array<double,3>,2> r_vec = ComputeRVariation(field,loc,nbor,r_nbor,l_nbor); //consecutive variation ratios

  //psi+ limiter
  array<double,3> psi_plus;
  for (int n=0;n<3;n++)
    psi_plus[n] = std::max(0.0,std::max(std::min(beta*r_vec[0][n],1.0),std::min(r_vec[0][n],beta)));

  //psi- limiter
  array<double,3> psi_minus;
  for (int n=0;n<3;n++){
    psi_minus[n] = std::max(0.0,std::max(std::min(beta*r_vec[1][n],1.0),std::min(r_vec[1][n],beta))); //need 2 max fcns. b/c std::max only takes 2 arguments by default
 
  }

  //combining + and - psi limiter
  array<array<double,3>,2> psi_total{psi_plus,psi_minus};

  return psi_total;

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeVanLeerLimiter(array<double,3> &r_vec){

  //Note: Refer to Class Notes Section 7, Page 20
  
  //array<array<double,3>,2> r_vec = ComputeRVariation(field,loc,nbor,r_nbor,l_nbor); //consecutive variation ratios
  //array<double,3> r_plus = ComputeRPlusVariation(field,loc,r_nbor,nbor);
  //array<double,3> r_minus = ComputeRMinusVariation(field,loc,l_nbor,nbor);

  //psi limiter
  array<double,3> psi;
  for (int n=0;n<3;n++)
    psi[n] = std::max(0.0,(r_vec[n] + abs(r_vec[n])) / (1.0+r_vec[n]));


  return psi;

}

//-----------------------------------------------------------
array<array<double,3>,2> Euler1D::ComputeRVariation(vector<array<double,3>>* &field,int loc,int nbor,int r_nbor,int l_nbor){

  double num,denom;
  double delta = 1.0e-6; //to prevent 0/0
  //r+ vector
  array<double,3> r_plus;
  if (r_nbor == total_cellnum){ //boundary conditions case for r+ - outflow 
    for (int n=0;n<3;n++) 
      r_plus[n] = 1.0; //TODO: assumes the consecutive variation is 1 due to the outflow BC
  }
  else{ //normal case
    for (int n=0;n<3;n++){ 
      num = (*field)[r_nbor][n] - (*field)[nbor][n];
      denom = (*field)[nbor][n] - (*field)[loc][n];
      denom = std::copysign(std::max(abs(denom),delta),denom); 
      r_plus[n] = num / denom;
    }
  }

  //r- vector
  array<double,3> r_minus;
  if (loc==0){ //boundary condition case for r- - inflow
    for (int n=0;n<3;n++) 
      r_minus[n] = 1.0; //TODO: assumes the consecutive variation is 1 due to the outflow BC
  }
  else{ //normal case
    for (int n=0;n<3;n++) {
      num = (*field)[loc][n] - (*field)[l_nbor][n];
      denom = (*field)[nbor][n] - (*field)[loc][n];
      denom = std::copysign(std::max(abs(denom),delta),denom); 
      r_minus[n] = num / denom;

    }
  }

  //total r vectors
  array<array<double,3>,2> total_r{r_plus,r_minus};

  return total_r;

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeRPlusVariation(vector<array<double,3>>* &field,int loc,int r_nbor,int nbor){

  array<double,3> r_plus;
  double num,denom;
  double delta = 1.0e-6; //to prevent 0/0

  for (int n=0;n<3;n++){ 
    num = (*field)[r_nbor][n] - (*field)[nbor][n];
    denom = (*field)[nbor][n] - (*field)[loc][n];
    denom = std::copysign(std::max(abs(denom),delta),denom); 
    r_plus[n] = num / denom;
  }

  return r_plus;

}

//-----------------------------------------------------------
array<double,3> Euler1D::ComputeRMinusVariation(vector<array<double,3>>* &field,int loc,int l_nbor,int nbor){

  array<double,3> r_minus;
  double num,denom;
  double delta = 1.0e-6; //to prevent 0/0

  for (int n=0;n<3;n++) {
    num = (*field)[loc][n] - (*field)[l_nbor][n];
    denom = (*field)[nbor][n] - (*field)[loc][n];
    denom = std::copysign(std::max(abs(denom),delta),denom); 
    r_minus[n] = num / denom;
  }

  return r_minus;

}

//-----------------------------------------------------------
double Euler1D::ComputeSourceTerm(vector<array<double,3>>* &field,int loc,vector<double> &xcoords,double dx) {

  //Note: loc is index of cell in Field -> xcoords index is i-2 of left face
  //Get area for i+1/2 and i-1/2 locations (refer to read me for data indexing)
  //dx = abs(xcoords[0]-xcoords[1]); //assigning dx val. here
  double A_rface = Tools::AreaVal(xcoords[loc-1]);
  double A_lface = Tools::AreaVal(xcoords[loc-2]);
  // pressure at i = loc
  double p = (*field)[loc][2];
  double dAdx = (A_rface-A_lface) / dx;
  // source = P_i * (A_i+1/2 - A_i-1/2) / dx * dx
  double res = p * dAdx;

  return res;

}

//-----------------------------------------------------------
double Euler1D::GetEpsilon2(vector<array<double,3>>* &field,int loc) {

  double Nu = GetNu(field,loc);
  double Nuleft = GetNu(field,loc-1);
  double Nuright = GetNu(field,loc+1);
  double Nuright2 = GetNu(field,loc+2);

  double kappa2 = 1.0/2.0; //typically from 1/4<kappa2<1/2
  
  double max = std::max({Nu,Nuleft,Nuright,Nuright2}); //acquring max nu
  double res = kappa2 * max;
  return res;

}


//-----------------------------------------------------------
double Euler1D::GetNu(vector<array<double,3>>* &field,int loc){

  int total_cells = (int)field->size();
  double res;
  //Tools::print("loc:%d\n",loc);

  if (loc == total_cells-1){ //check for evaluating nu at last outflow ghost cell
    res = (*field)[loc][2] - (2.0* (*field)[loc][2]) +  (*field)[loc-1][2]; //assuming Pressure is same at cell loc
    res /= (*field)[loc][2] + (2.0* (*field)[loc][2]) +  (*field)[loc-1][2]; //assuming Pressure is same at cell loc
    res = 0.0;
    return res;
  }

  if (loc == 0){ //check for evaluating new at 1st inflow ghost cell
    res = (*field)[loc+1][2] - (2.0* (*field)[loc][2]) +  (*field)[loc][2]; //assuming Pressure is same at cell loc
    res /= (*field)[loc+1][2] + (2.0* (*field)[loc][2]) +  (*field)[loc][2]; //assuming Pressure is same at cell loc
    return res;
  }

  res = (*field)[loc+1][2] - (2.0* (*field)[loc][2]) +  (*field)[loc-1][2];
  res /= (*field)[loc+1][2] + (2.0* (*field)[loc][2]) +  (*field)[loc-1][2];
  res = abs(res);

  return res;

}


//-----------------------------------------------------------
double Euler1D::GetMachNumber(vector<array<double,3>>* &field,int loc){

  //Using insentropic conditions only for thermodynamic variables
  double T = (*field)[loc][2] / ((*field)[loc][0]*R); //if T is negative than M will be -nan!!!

  //using actual x-vel. value now
  double M = (*field)[loc][1] / sqrt(gamma * R * T); //M = u/a

  return M;

}


//-----------------------------------------------------------
double Euler1D::GetLambda(vector<array<double,3>>* &field,int loc){

  double lambda_i = GetLambdaMax(field,loc); //Lambda bar

  double lambda_iright = GetLambdaMax(field,loc+1); //Lambda_right bar

  double res = (lambda_i + lambda_iright) / 2.0;
  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute2ndOrderDamping(vector<array<double,3>>* &field,int loc){

  //following Roy's class notes nomenclature (section 3 slide 31)
  //returns solely \arrow{d^2} vector!
  //look into applying damping terms to boundary?
  array<double,3> conserved = ComputeConserved(field,loc);
  array<double,3> conserved_nbor = ComputeConserved(field,loc+1);
  
  double lambda = GetLambda(field,loc); //at cell face (i+1/2)
  double epsilon = GetEpsilon2(field,loc); //sensor for detecting shocks (will have to tweak the constant later)

  array<double,3> res;
   for (int i=0;i<3;i++) //computing 2nd order damping here
    res[i] = lambda*epsilon*(conserved_nbor[i]-conserved[i]);
  
  return res;

}



//-----------------------------------------------------------
double Euler1D::GetEpsilon4(vector<array<double,3>>* &field,int loc){

  double epsilon2 = GetEpsilon2(field,loc);
  double kappa4 = 1.0/30.0; //typically ranges from: 1/64<kappa4<1/32
  double res = std::max(0.0,(kappa4 - epsilon2));

  return res;

}


//-----------------------------------------------------------
array<double,3> Euler1D::Compute4thOrderDamping(vector<array<double,3>>* &field,int loc){

  double lambda = GetLambda(field,loc);
  double epsilon = GetEpsilon4(field,loc); //looks for gradients that cause odd-even decoupling (will have to tweak the constant later)
  //TODO: Convert to conservative variables HERE ONLY!
  array<double,3> conserved = ComputeConserved(field,loc);
  array<double,3> conserved_leftnbor = ComputeConserved(field,loc-1);
  array<double,3> conserved_rightnbor = ComputeConserved(field,loc+1);
  array<double,3> conserved_right2nbor = ComputeConserved(field,loc+2);

  array<double,3> res;
  for (int i=0;i<3;i++) //computing 4th order damping term here
    res[i] = lambda*epsilon*(conserved_right2nbor[i] - 3.0*conserved_rightnbor[i] + 3.0*conserved[i] - conserved_leftnbor[i]);


  return res;
}


//-----------------------------------------------------------
void Euler1D::ComputeResidual(vector<array<double,3>>* &resid,vector<array<double,3>>* &field,vector<array<double,3>>* &field_stall,vector<double> &xcoords,double &dx,bool flux_scheme,bool flux_accuracy,bool upwind_scheme,double epsilon,bool &resid_stall){ 

  //following nomenclature from class notes
  [[maybe_unused]] array<double,3> F_right,F_left; //left and right face spatial fluxes 
  [[maybe_unused]] array<double,3> D2_right,D2_left,D4_right,D4_left; //left and right face damping terms
  array<double,3> TotalF_right,TotalF_left; //left and right total fluxes (spatial + damping)
  array<double,3> Source{0.0,0.0,0.0}; //source term vector
  double S; //source term (only for x-mom. eq.)
  double A_left,A_right; //Area of corresponding faces

  for (int n=0;n<total_cellnum;n++){ //looping through all interior nodes
    if ((n==0) | (n==1) | (n==total_cellnum-2) | (n==total_cellnum-1)) //skipping the ghost cell nodes
      continue;

  //if (n == interior_cellnum - 1) //debug
    //Tools::print("I am here\n");  

    //Fluxes Evaluation
    //note: \arrow{F}_(i-1/2) is same as \arrow{F}_(i+1/2) of cell to the left!
    if (flux_scheme == true){ //JST Damping Case -- central quadrature
      // Spatial Fluxes
       F_right = ComputeSpatialFlux_BASE(field,n,n+1);
       F_left = ComputeSpatialFlux_BASE(field,n-1,n);

      //JST Damping Terms (need a D2_left flux and D2_right flux vector; similar for D4)
      //note: \arrow{D}_(i-1/2) is same as \arrow{D}_(i+1/2) of cell to the left!
      // right face
      D2_right = Compute2ndOrderDamping(field,n);
      D4_right = Compute4thOrderDamping(field,n);
      // left face
      D2_left = Compute2ndOrderDamping(field,n-1);
      D4_left = Compute4thOrderDamping(field,n-1);

      //Total Flux Terms
      for (int i=0;i<3;i++){
        TotalF_right[i] = F_right[i] - (D2_right[i]-D4_right[i]);
        TotalF_left[i] = F_left[i] - (D2_left[i]-D4_left[i]);
      }
    }

    //TODO: Upwind Cases
    else{
      if (flux_accuracy == true){ //1st order accurate case
        TotalF_right = ComputeSpatialFlux_UPWIND1stOrder(field,upwind_scheme,n,n+1);
        TotalF_left = ComputeSpatialFlux_UPWIND1stOrder(field,upwind_scheme,n-1,n);

      }
      else if (flux_accuracy == false){ //2nd order accurate case
        TotalF_right = ComputeSpatialFlux_UPWIND2ndOrder(field,field_stall,upwind_scheme,n,n+1,epsilon,resid_stall);
        TotalF_left = ComputeSpatialFlux_UPWIND2ndOrder(field,field_stall,upwind_scheme,n-1,n,epsilon,resid_stall);
      }
    }

    //Source Term (external pressure) ONLY for x-mom. eq.
    // also, area already evaluated, but may need to be multiplied dx?
    S = ComputeSourceTerm(field,n,xcoords,dx);
    Source[1] = S; //adding scalar to vector


    //Area Evaluations
    A_left = Tools::AreaVal(xcoords[n-2]);
    A_right = Tools::AreaVal(xcoords[n-1]);

    //Residual cal.
    for (int i=0;i<3;i++)
      (*resid)[n-2][i] = (TotalF_right[i]*A_right - TotalF_left[i]*A_left) - Source[i]*dx;
    

  }
  return;

}

//-----------------------------------------------------------
double Euler1D::GetLambdaMax(vector<array<double,3>>* &field,int loc){

  double M = GetMachNumber(field,loc); //cell-averaged Mach number
  double a = abs((*field)[loc][1]) / M; //cell-averaged speed of sound
  double lambda_max = abs((*field)[loc][1]) + a;

  return lambda_max;

}

//-----------------------------------------------------------
double Euler1D::GetCellAverageSol(double &A_left,double &A_right,double &dx,array<double,3> &sol_left,array<double,3> &sol_right){

  //using weighted geometric cell-average formula
  double num = 0.5*(A_right + A_left) * (sol_right[1]+sol_left[1]);
  double vol = 0.5*(A_right + A_left) * dx;

  double sol = num / vol;

  return sol;

}

//-----------------------------------------------------------
double Euler1D::ComputeMachNumber(array<double,3> &sols){

  //Using insentropic conditions only for thermodynamic variables
  double T = sols[2] / (sols[0]*R); //if T is negative than M will be -nan!!!

  //using actual x-vel. value now
  double M = sols[1] / sqrt(gamma * R * T); //M = u/a

  return M;
}

//-----------------------------------------------------------
Euler1D::~Euler1D(){}
//-----------------------------------------------------------

// EULER2D DEFINITIONS
//-----------------------------------------------------------
Euler2D::Euler2D(int case_2d,int cell_inum,int cell_jnum,int scheme,int accuracy,int top,int btm,int left,int right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source) : EulerBASE(cell_inum,cell_jnum,scheme,accuracy,top,btm,left,right,mesh_ptr,source){
  //Case assignments
  if (case_2d == 0){ //30 deg inlet case
    Mach_bc = 4.0;
    P_stag = 12270.0;
    T_stag = 217.0;
    alpha = 0.0;

  }
  else if (case_2d == 1){ // 0 deg AOA airfoil
    Mach_bc = 0.84;
    P_stag = 65855.8;
    T_stag = 300.0;
    alpha = 0.0;

  }
  else if (case_2d == 2){ // 8 deg AOA airfoil
    Mach_bc = 0.75;
    P_stag = 67243.5;
    T_stag = 300.0;
    alpha = 8.0;

  }
  else {
    cerr<<"Unknown 2D Case!"<<endl;
  }

  return;
}
//-----------------------------------------------------------
void Euler2D::InitSolutions(vector<array<double,4>>* &field,int cellnum){
  
  //Setting w/ arb. sin function
  for (int i=0;i<cellnum;i++) {
    (*field)[i][0] = 1.0 + sin(i);
    (*field)[i][1] = 2.0 + sin(i);
    (*field)[i][2] = 3.0 + sin(i);
    (*field)[i][3] = 4.0 + sin(i);
  }

  return;
}
//-----------------------------------------------------------
void Euler2D::SetInitialConditions(vector<array<double,4>>* &field){

  int cellnum = cell_imax * cell_jmax;
  //Setting w/ arb. sin function
  for (int i=0;i<cellnum;i++) {
    (*field)[i][0] = 1.0 + sin(i);
    (*field)[i][1] = 2.0 + sin(i);
    (*field)[i][2] = 3.0 + sin(i);
    (*field)[i][3] = 4.0 + sin(i);
  }

  return;

}

//-----------------------------------------------------------
array<double,4> Euler2D::ComputeSpatialFlux_UPWIND1stOrder(vector<array<double,4>>* field,int loci,int locj,int nbori,int nborj){
  
  //NOTE: assumes loc is left state and nbor is right state (for 1st order accuracy)
  array<double,4> loc_state = fieldij(field,loci,locj,cell_imax);
  array<double,4> nbor_state = fieldij(field,nbori,nborj,cell_imax);

  array<double,4> flux; //total flux
  [[maybe_unused]] array<double,4> flux_rtstate; //right state flux (for upwinding in +c wave speed)
  [[maybe_unused]] array<double,4> flux_ltstate; //left state flux (for upwinding in -c wave speed)


  if (flux_scheme == 1){ //Van Leer Method
    //flux_rtstate = VanLeerCompute(right_state,false); //false for negative c case
    //flux_ltstate = VanLeerCompute(left_state,true); //true for positive c case
   
    for (int n=0;n<3;n++) //summing up left and right state fluxes
      flux[n] = flux_rtstate[n] + flux_ltstate[n];

  }
  else if (flux_scheme == 2){ //Roe's Method
    //flux = ComputeRoeFlux(left_state,right_state);
  }
  else { //error handling
    cerr<<"Error: Unknown flux specification!"<<endl;
  }

  return flux;

}

//-----------------------------------------------------------
array<double,4> Euler2D::ComputeSpatialFlux_UPWIND2ndOrder(vector<array<double,4>>* field,int loci,int locj,int nbori,int nborj){}



//-----------------------------------------------------------
void Euler2D::ComputeResidual(vector<array<double,4>>* &resid,vector<array<double,4>>* &field){

  array<double,4> flux_right,flux_left,flux_top,flux_btm;
  double area_right,area_left,area_btm,area_top;
  array<double,4> res;
  array<double,4> field_loc;
  array<double,4> field_rnbor,field_lnbor,field_topnbor,field_btmnbor;

  //normal residual computation (no source term)
  for (int j=0;j<cell_jmax;j++){
    for (int i=0;i<cell_imax;i++){
  
      //TODO:fluxes evaluation
      if (flux_scheme==0){ //JST Damping
      }
      else if ((flux_scheme!=0) && (flux_accuracy==1)){ //1st order Upwind Schemes
        flux_right = ComputeSpatialFlux_UPWIND1stOrder(field,i,j,i+1,j);
        flux_left = ComputeSpatialFlux_UPWIND1stOrder(field,i-1,j,i,j);
        flux_top = ComputeSpatialFlux_UPWIND1stOrder(field,i,j,i,j+1);
        flux_btm = ComputeSpatialFlux_UPWIND1stOrder(field,i,j-1,i,j);
      }
      else if ((flux_scheme!=0) && (flux_accuracy==2)){ //2nd ordr Upwind Schemes w/ MUSCL extrapolation
        flux_right = ComputeSpatialFlux_UPWIND2ndOrder(field,i,j,i+1,j);
        flux_left = ComputeSpatialFlux_UPWIND2ndOrder(field,i-1,j,i,j);
        flux_top = ComputeSpatialFlux_UPWIND2ndOrder(field,i,j,i,j+1);
        flux_btm = ComputeSpatialFlux_UPWIND2ndOrder(field,i,j-1,i,j);
      }
      else {} //error handling

      //area evaluation
      area_right = mesh->GetInteriorCellArea(i,j,3);
      area_left = mesh->GetInteriorCellArea(i,j,2);
      area_top = mesh->GetInteriorCellArea(i,j,0);
      area_btm = mesh->GetInteriorCellArea(i,j,1);

      //residual calc.
      for (int n=0;n<cell_imax;n++)
        res[n] = (flux_right[n]*area_right-flux_left[n]*area_left) + (flux_top[n]*area_top - flux_btm[n]*area_btm);

      fieldij(resid,i,j,cell_imax) = res;
    }
  }

  return;
}

//-----------------------------------------------------------
Euler2D::~Euler2D(){}
//-----------------------------------------------------------

// EULER2DMMS DEFINITIONS
//-----------------------------------------------------------
Euler2DMMS::Euler2DMMS(int cell_inum,int cell_jnum,int scheme,int accuracy,int top,int btm, int left,int right,MeshGenBASE* &mesh_ptr,vector<array<double,4>>* &source) : EulerBASE(cell_inum,cell_jnum,scheme,accuracy,top,btm,left,right,mesh_ptr,source){}
//-----------------------------------------------------------
void Euler2DMMS::SetInitialConditions(vector<array<double,4>>* &field){

  int cellnum = cell_imax * cell_jmax;

  for (int i=0;i<cellnum;i++) {
    (*field)[i][0] = 0.0;
    (*field)[i][1] = 0.0;
    (*field)[i][2] = 0.0;
    (*field)[i][3] = 0.0;
  }

  return;
}

//-----------------------------------------------------------
void Euler2DMMS::ManufacturedPrimitiveSols(vector<array<double,4>>* &field,SpaceVariables2D* &sols){

  int cellid;
  double x,y; //x and y coords

  cell_imax = mesh->Nx - 1;
  cell_jmax = mesh->Ny - 1;

  //Evaluating Cell-Center Coords 
  vector<double> cell_center_xcoords(mesh->cellnumber);
  vector<double> cell_center_ycoords(mesh->cellnumber);

  sols->ComputeCellCenteredCoordinate(mesh->xcoords,mesh->ycoords,cell_center_xcoords,cell_center_ycoords,cell_imax);

  //int cell_imax = imax - 1; //last cell index in i-row
  //int cell_jmax = jmax - 1; //last cell index in j-row
  //int Nx = cell_imax + 1;
 


  // Solving Manufactured Solution
  
  for (int j=0;j<cell_jmax;j++){
    for (int i=0;i<cell_imax;i++){

      cellid = i + j*cell_imax; //accesses index in 1D flat vectors to retrieve x&y coords of pt
      x = cell_center_xcoords[cellid]; y = cell_center_ycoords[cellid];

      //Rho
      fieldij(field,i,j,cell_imax)[0] = rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L);

      //U
      fieldij(field,i,j,cell_imax)[1] = uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L));

      //V
      fieldij(field,i,j,cell_imax)[2] = vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L));

      //P
      fieldij(field,i,j,cell_imax)[3] = press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L);

    }
  }
  
  return;
}
//-----------------------------------------------------------
double Euler2DMMS::ContinuitySourceTerm(double x,double y){

  double source_term = (3.0 * Pi * uvelx * cos((3.0 * Pi * x) / (2.0 * L))) *
    (rho0 + rhoy * cos((Pi * y) / (2.0 * L)) + rhox * sin((Pi * x) / L)) / (2.0 * L)
  + (2.0 * Pi * vvely * cos((2.0 * Pi * y) / (3.0 * L))) *
    (rho0 + rhoy * cos((Pi * y) / (2.0 * L)) + rhox * sin((Pi * x) / L)) / (3.0 * L)
  + (Pi * rhox * cos((Pi * x) / L) *
     (uvel0 + uvely * cos((3.0 * Pi * y) / (5.0 * L)) + uvelx * sin((3.0 * Pi * x) / (2.0 * L)))) / L
  - (Pi * rhoy * sin((Pi * y) / (2.0 * L)) *
     (vvel0 + vvelx * cos((Pi * x) / (2.0 * L)) + vvely * sin((2.0 * Pi * y) / (3.0 * L)))) / (2.0 * L);


  return source_term;

}

//-----------------------------------------------------------
double Euler2DMMS::XMomentumSourceTerm(double x,double y){

  double source_term = (2.0 * Pi * vvely * cos((2.0 * Pi * y) / (3.0 * L))) *
    (rho0 + rhoy * cos((Pi * y) / (2.0 * L)) + rhox * sin((Pi * x) / L)) *
    (uvel0 + uvely * cos((3 * Pi * y) / (5.0 * L)) + uvelx * sin((3.0 * Pi * x) / (2.0 * L))) / (3.0 * L)
    
  - (Pi * rhoy * (uvel0 + uvely * cos((3.0 * Pi * y) / (5.0 * L)) + uvelx * sin((3.0 * Pi * x) / (2.0 * L))) *
     sin((Pi * y) / (2.0 * L)) * (vvel0 + vvelx * cos((Pi * x) / (2.0 * L)) + vvely * sin((2.0 * Pi * y) / (3.0 * L)))) / (2.0 * L)
    
  - (3.0 * Pi * uvely * (rho0 + rhoy * cos((Pi * y) / (2.0 * L)) + rhox * sin((Pi * x) / L)) *
     sin((3.0 * Pi * y) / (5.0 * L)) * (vvel0 + vvelx * cos((Pi * x) / (2.0 * L)) + vvely * sin((2.0 * Pi * y) / (3.0 * L)))) / (5.0 * L);

  return source_term;

}

//-----------------------------------------------------------
double Euler2DMMS::YMomentumSourceTerm(double x,double y){

  double source_term = (Pi * pressy * cos((Pi * y) / L)) / L
  + (4.0 * Pi * vvely * cos((2.0 * Pi * y) / (3.0 * L))) *
    (rho0 + rhoy * cos((Pi * y) / (2.0 * L)) + rhox * sin((Pi * x) / L)) *
    (vvel0 + vvelx * cos((Pi * x) / (2.0 * L)) + vvely * sin((2.0 * Pi * y) / (3.0 * L))) / (3.0 * L)
  - (Pi * rhoy * sin((Pi * y) / (2.0 * L)) *
     pow(vvel0 + vvelx * cos((Pi * x) / (2.0 * L)) + vvely * sin((2.0 * Pi * y) / (3.0 * L)), 2.0)) / (2.0 * L);

  return source_term;

}

//-----------------------------------------------------------
double Euler2DMMS::EnergySourceTerm(double x,double y){

  double source_term = 
(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)))*
    ((-2.0*Pi*pressx*sin((2.0*Pi*x)/L))/L + (rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))*
       ((-2.0*Pi*pressx*sin((2.0*Pi*x)/L))/((-1.0 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))) + 
         ((3.0*Pi*uvelx*cos((3.0*Pi*x)/(2.0*L))*(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L))))/L - 
            (Pi*vvelx*sin((Pi*x)/(2.0*L))*(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L))))/L)/2.0 - 
         (Pi*rhox*cos((Pi*x)/L)*(press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L)))/
          ((-1.0 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L),2.0))) + 
      (Pi*rhox*cos((Pi*x)/L)*((pow(wvel0,2.0) + pow(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)),2.0) + 
              pow(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L)),2.0))/2.0 + 
           (press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L))/
            ((-1.0 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L)))))/L) + 
   (3.0*Pi*uvelx*cos((3.0*Pi*x)/(2.0*L))*(press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L) + 
        (rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))*
         ((pow(wvel0,2.0) + pow(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)),2.0) + 
              pow(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L)),2.0))/2.0 + 
           (press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L))/
            ((-1.0 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))))))/(2.0*L) + 
   (2.0*Pi*vvely*cos((2.0*Pi*y)/(3.0*L))*(press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L) + 
        (rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))*
         ((pow(wvel0,2.0) + pow(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)),2.0) + 
              pow(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L)),2.0))/2.0 + 
           (press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L))/
            ((-1.0 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))))))/(3.0*L) + 
   (vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L)))*
    ((Pi*pressy*cos((Pi*y)/L))/L - (Pi*rhoy*sin((Pi*y)/(2.0*L))*
         ((pow(wvel0,2.0) + pow(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)),2.0) + 
              pow(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L)),2.0))/2.0 + 
           (press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L))/
            ((-1.0 + gamma)*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L)))))/(2.0*L) + 
      (rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))*
       ((Pi*pressy*cos((Pi*y)/L))/((-1.0 + gamma)*L*(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L))) + 
         ((-6.0*Pi*uvely*(uvel0 + uvely*cos((3.0*Pi*y)/(5.0*L)) + uvelx*sin((3.0*Pi*x)/(2.0*L)))*sin((3.0*Pi*y)/(5.0*L)))/(5.0*L) + 
            (4*Pi*vvely*cos((2.0*Pi*y)/(3.0*L))*(vvel0 + vvelx*cos((Pi*x)/(2.0*L)) + vvely*sin((2.0*Pi*y)/(3.0*L))))/(3.0*L))/2.0 + 
         (Pi*rhoy*sin((Pi*y)/(2.0*L))*(press0 + pressx*cos((2.0*Pi*x)/L) + pressy*sin((Pi*y)/L)))/
          (2.0*(-1.0 + gamma)*L*pow(rho0 + rhoy*cos((Pi*y)/(2.0*L)) + rhox*sin((Pi*x)/L),2.0))));


  return source_term;

}

//-----------------------------------------------------------
void Euler2DMMS::EvalSourceTerms(/*vector<array<double,4>>* &mms_source,*/SpaceVariables2D* &sols){

  double x,y;
  int cellid;
  cell_imax = mesh->Nx-1;
  cell_jmax = mesh->Ny-1;
  //Evaluating Cell-Center Coords 
  vector<double> cell_center_xcoords(mesh->cellnumber);
  vector<double> cell_center_ycoords(mesh->cellnumber);

  sols->ComputeCellCenteredCoordinate(mesh->xcoords,mesh->ycoords,cell_center_xcoords,cell_center_ycoords,cell_imax);

  //Evaluating Source Terms
  
  for (int j=0;j<cell_jmax;j++){
    for (int i=0;i<cell_imax;i++){

      cellid = i + j*cell_imax; //accesses index in 1D flat vectors to retrieve x&y coords of pt
      x = cell_center_xcoords[cellid]; y = cell_center_ycoords[cellid];

      //Mass Conservation
      fieldij(mms_source,i,j,cell_imax)[0] = ContinuitySourceTerm(x,y);

      //X-Momentum
      fieldij(mms_source,i,j,cell_imax)[1] = XMomentumSourceTerm(x,y);

      //Y-Momentum
      fieldij(mms_source,i,j,cell_imax)[2] = YMomentumSourceTerm(x,y);

      //Energy
      fieldij(mms_source,i,j,cell_imax)[3] = EnergySourceTerm(x,y);

    }
  }

}

//-----------------------------------------------------------
void Euler2DMMS::ComputeResidual(vector<array<double,4>>* &resid,vector<array<double,4>>* &field){


  array<double,4> flux_right,flux_left,flux_top,flux_btm;
  double area_right,area_left,area_btm,area_top;
  double vol;
  array<double,4> res;
  array<double,4> field_loc;
  array<double,4> field_rnbor,field_lnbor,field_topnbor,field_btmnbor;

  array<double,4> mms;
  

  for (int j=0;j<cell_jmax;j++){
    for (int i=0;i<cell_imax;i++){
      //primitive vars. of current cell and nbor cells
      field_loc = fieldij(field,i,j,cell_imax); 
      field_rnbor = fieldij(field,i+1,j,cell_imax); 
      field_lnbor = fieldij(field,i-1,j,cell_imax); 
      field_topnbor = fieldij(field,i,j+1,cell_imax); 
      field_btmnbor = fieldij(field,i,j-1,cell_imax); 
  
      //TODO:fluxes evaluation
      //flux_right = (scheme==1) ? VanLeer(accuracy,field_loc,field_rnbor,mesh_ptr) : Roe(accuracy,field_loc,field_rnbor,mesh_ptr);
      //flux_left = (scheme==1) ? VanLeer(accuracy,field_loc,field_lnbor) : Roe(accuracy,field_loc,field_lnbor);
      //flux_top = (scheme==1) ? VanLeer(accuracy,field_loc,field_topnbor) : Roe(accuracy,field_loc,field_topnbor);
      //flux_btm = (scheme==1) ? VanLeer(accuracy,field_loc,field_btmnbor) : Roe(accuracy,field_loc,field_btmnbor);
      
      //volume evaluation (for source term) + source term retrievel 
      vol = mesh->GetCellVolume(i,j);
      mms = fieldij(mms_source,i,j,cell_imax);

      //area evaluation
      area_right = mesh->GetInteriorCellArea(i,j,3);
      area_left = mesh->GetInteriorCellArea(i,j,2);
      area_top = mesh->GetInteriorCellArea(i,j,0);
      area_btm = mesh->GetInteriorCellArea(i,j,1);

      //residual calc.
      for (int n=0;n<4;n++)
        res[n] = (flux_right[n]*area_right-flux_left[n]*area_left) + (flux_top[n]*area_top - flux_btm[n]*area_btm) - mms[n]*vol;

      fieldij(resid,i,j,cell_imax) = res;
    }
  }

  return;
}

//-----------------------------------------------------------
Euler2DMMS::~Euler2DMMS(){}
