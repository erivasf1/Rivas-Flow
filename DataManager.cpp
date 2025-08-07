//User-defined functions
#include "DataManager.h" 

// SPACEVARIABLE1D DEFINITIONS

SpaceVariables1D::SpaceVariables1D()
{}

//---------------------------------------------------------
array<double,3> SpaceVariables1D::ComputeSolutionNorms(vector<array<double,3>>* &resid){

  array<double,3> norm{0.0,0.0,0.0};

  double imax = (double)resid->size(); //number of interior nodes 

  //using L2 norm  
  //Tools::print("Calculating Global Norm\n");
  for (int i=0;i<(int)resid->size();i++){
    //Tools::print("Cell number: %d\n",i);
    //continuity
    //Tools::print("continuity res.: %e\n",Resid[i][0]);
    norm[0]+= pow((*resid)[i][0],2);
    //x-mom
    //Tools::print("x-mom. res.: %e\n",Resid[i][1]);
    norm[1]+= pow((*resid)[i][1],2);
    //energy
    //Tools::print("energy res.: %e\n",Resid[i][2]);
    norm[2]+= pow((*resid)[i][2],2);

  }
  norm[0] = sqrt(norm[0]/imax); //normalizing by imax
  norm[1] = sqrt(norm[1]/imax);
  norm[2] = sqrt(norm[2]/imax);

  //Tools::print("Continuity res.: %e\n",norm[0]);
  //Tools::print("X-Momentum res.: %e\n",norm[1]);
  //Tools::print("Energy res.: %e\n",norm[2]);

  return norm;
}

//---------------------------------------------------------
double SpaceVariables1D::ComputeNormAvg(array<double,3> &Norms){

  double norms_avg = (Norms[0]+Norms[1]+Norms[2]) / 3.0;
  return norms_avg;

}

//---------------------------------------------------------
double SpaceVariables1D::ComputeRampValue(array<double,3> CurrentNorms,array<double,3> InitNorms,double FinalVal){
  
  //Note: using a log10 function as the ramping function
  
  double InitVal = ComputeNormAvg(InitNorms);
  double CurrentVal = ComputeNormAvg(CurrentNorms);
  double p = 30.0; //used to accelerate or deaccelerate the ramping fcn.

  double ramp_val = (log10(InitVal) - log10(CurrentVal)) / (log10(CurrentVal) - log10(FinalVal));
  ramp_val = pow(ramp_val,p);

  ramp_val = std::max(0.0,std::min(1.0,ramp_val));

  return ramp_val;

}

//---------------------------------------------------------
void SpaceVariables1D::OutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,const char *filename){


  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
    double M; //temp. variable for Mach number
    myfile<<"Primitive Variable Solutions (Including Ghost Cells)"<<endl;
    myfile<<"Cell#"<<"  "<<"Density(kg/m^3)"<<"  "<<"Velocity(m/s)"<<"  "<<"Pressure(Pa)"<<"  "<<"Mach Number"<<endl;
    for (int i=0;i<(int)field->size();i++){
  
      //Computing Mach Number for checking the initial conditions
      M = euler->GetMachNumber(field,i);

      myfile<<i<<"  "<<(*field)[i][0]<<"  "<<(*field)[i][1]<<"  "<<(*field)[i][2]<<"  "<<M<<endl;

    }

  
  myfile.close(); //closing file writing to it
  //myfile.flush();

return;
}

//---------------------------------------------------------
void SpaceVariables1D::AllOutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,string filename,bool cond,int iter,vector<double> &xcoords){

  cell_number = (int)field->size()-4; //number of interior cells
  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false)
    myfile<<"variables= \"cell index\" \"rho(kg/m^3)\" \"u(m/s)\"  \"Press(N/m^2)\" \"Mach\" \"Xcoords\""<<endl;

//Repeat the following each time you want to write out the solution
/*
write(40,*) 'zone T="',num_iter,'" '
write(40,*) 'I=',imax
write(40,*) 'DATAPACKING=POINT'
write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
& DOUBLE DOUBLE )'
  */

  myfile<<"zone T= "<<"\""<<iter<<"\""<<endl;
  myfile<<"I="<<cell_number<<endl;
  myfile<<"DATAPACKING=POINT"<<endl;
  myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;

  double M;

  for (int n=0;n<cell_number;n++){
    //Computing Mach Number for checking the initial conditions
    M = euler->GetMachNumber(field,n+2);
    myfile<<n+1<<"  ";
    myfile<<(*field)[n][0]<<"  "<<(*field)[n][1]<<"  "<<(*field)[n][2]<<"  "<<M<<"  "<<xcoords[n]<<endl;

  }

  
  
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}

//---------------------------------------------------------
void SpaceVariables1D::OutputLocalResiduals(vector<array<double,3>> &Resid,const char *filename){


  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
  myfile<<"Local Residuals"<<endl;
  myfile<<"Cell"<<"  "<<"Continuity"<<"  "<<"Momentum"<<"  "<<"Energy"<<endl;
  for (int n=0;n<(int)Resid.size();n++){
    myfile<<n<<"  "<<Resid[n][0]<<"  "<<Resid[n][1]<<"  "<<Resid[n][2]<<endl;
  }
  
  myfile.close(); //closing file writing to it
  //myfile.flush();

return;
}


//---------------------------------------------------------
void SpaceVariables1D::OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename){

  ofstream myfile;
  myfile.open(filename);
  vector<string> names{"continuity","x-momentum","energy"};
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

   array<double,3> Source{0.0,S,0.0}; //array vector of sorce terms

    for (int n=0;n<3;n++){
      myfile<<names[0]<<" Residual Fluxes"<<endl;
      myfile<<"F_right: "<<F_right[n]<<"  "<<"F_left: "<<F_left[n]<<"  "<<"S: "<<Source[n]<<"  "<<"D2_right"<<D2_right[n]<<"  "<<"D2_left"<<D2_left[n]<<"  "<<"D4_right"<<D4_right[n]<<"  "<<"D4_left"<<D4_left[n]<<endl;
  
    }
  
  myfile.close(); //closing file writing to it
  //myfile.flush();
   



  return;
}

//---------------------------------------------------------
void SpaceVariables1D::ComputeCellAveragedSol(vector<array<double,3>>* &cell_faces,vector<array<double,3>>* &cell_sols,vector<double> &xcoords){

  //SolFace has one more element than SolCell, due to more faces than cells
  // using trapezoidal to approximate integral
  int rght_face,lft_face; //indexes
  double A_right,A_left;
  for (int n=0;n<(int)cell_sols->size();n++){ //computing sol. at each cell
    lft_face = n;
    rght_face = n + 1; 
    A_right = Tools::AreaVal(xcoords[rght_face]);
    A_left = Tools::AreaVal(xcoords[lft_face]);

    for (int i=0;i<3;i++){ //computing cell average for each primitive variable
      (*cell_sols)[n][i] = A_right* (*cell_faces)[rght_face][i] + A_left* (*cell_faces)[lft_face][i];
      (*cell_sols)[n][i] /= A_right + A_left;
    }

  }

  return;
}

//-----------------------------------------------------------------------
SpaceVariables1D::~SpaceVariables1D(){}
//-----------------------------------------------------------------------

// SPACEVARIABLE2D DEFINITIONS

//-----------------------------------------------------------------------
SpaceVariables2D::SpaceVariables2D()
{}
//-----------------------------------------------------------------------
void SpaceVariables2D::ComputeInteriorCellCenteredCoordinate(vector<double> &xcoords,vector<double> &ycoords,vector<double> &cell_center_xcoords,vector<double> &cell_center_ycoords,int cell_imax){ //computes coord. cell center

  //Cell-Center is just an avg. of the 4 corner nodes of the cell
  double x_avg = 0.0;
  double y_avg = 0.0;
  int id; //cell ID
  int nx_pts = cell_imax + 1; //# of points in i-row
  int node_id = 0;
  
  //Summing up x & y
  for (int cell_id=0;cell_id<(int)cell_center_xcoords.size();cell_id++){

    if ((cell_id % cell_imax == 0) && (cell_id>0)) //goes to next "j" row of cells if cell_imax is reached
      node_id++;
    //btm left node is n index
    id = node_id;
    x_avg += xcoords[id];
    y_avg += ycoords[id];
    //btm right node is n+1 index
    id = node_id+1;
    x_avg += xcoords[id];
    y_avg += ycoords[id];
    //top left node is (n + nx)
    id = node_id + nx_pts;
    x_avg += xcoords[id];
    y_avg += ycoords[id];
    //top right node is (n + nx_pts + 1)
    id = node_id + nx_pts + 1;
    x_avg += xcoords[id];
    y_avg += ycoords[id];

    //avg. nodes
    x_avg /= 4.0;
    y_avg /= 4.0;

    //save to cell_center_coords list
    cell_center_xcoords[cell_id] = x_avg;
    cell_center_ycoords[cell_id] = y_avg;
    
    //reset avg values
    x_avg = 0.0;
    y_avg = 0.0;

    node_id++;
    
 }


  return; 
}
//---------------------------------------------------------
void SpaceVariables2D::ComputeGhostCellCenteredCoordinate(MeshGen2D* &mesh){
  
  double x1=0.0; double y1=0.0;
  double x2=0.0; double y2=0.0;
  array<double,2> avg_xy1,avg_xy2;
  array<array<double,4>,2> cell_coords1,cell_coords2; //1&2 by doing corresponding sides at the same time (e.g. btm&top are done at the same time)

  //NOTE: cells coords data structure cell_coords[0]=all x vals. and vice versa for [1]

  //Btm + Top ghost cells
  for (int j=0;j<2;j++){
    for (int i=0;i<mesh->cell_imax;i++){
      cell_coords1 = mesh->GetGhostCellCoords(i,j,0); //top
      cell_coords2 = mesh->GetGhostCellCoords(i,j,1); //btm
      for (int n=0;n<4;n++){
        x1 += cell_coords1[0][n]; y1 += cell_coords1[1][n];
        x2 += cell_coords2[0][n]; y2 += cell_coords2[1][n];
      }
      avg_xy1[0] = x1/4.0; avg_xy1[1] = y1/4.0;
      avg_xy2[0] = x2/4.0; avg_xy2[1] = y2/4.0;
       
      mesh->top_cellcenter_coords[i+(j*mesh->cell_imax)] = avg_xy1;
      mesh->btm_cellcenter_coords[i+(j*mesh->cell_imax)] = avg_xy2;
    }
  }

  //Left + Right ghost cells
  for (int i=0;i<2;i++){
    for (int j=0;j<mesh->cell_jmax;j++){
      cell_coords1 = mesh->GetGhostCellCoords(i,j,2); //left
      cell_coords2 = mesh->GetGhostCellCoords(i,j,3); //right
      for (int n=0;n<4;n++){
        x1 += cell_coords1[0][n]; y1 += cell_coords1[1][n];
        x2 += cell_coords2[0][n]; y2 += cell_coords2[1][n];

      }
      avg_xy1[0] = x1/4.0; avg_xy1[1] = y1/4.0;
      avg_xy2[0] = x2/4.0; avg_xy2[1] = y2/4.0;
       
      mesh->left_cellcenter_coords[i+(j*mesh->cell_imax)] = avg_xy1;
      mesh->right_cellcenter_coords[i+(j*mesh->cell_imax)] = avg_xy2;
    }
  }
  return;

}
//---------------------------------------------------------
SpaceVariables2D::~SpaceVariables2D(){}
