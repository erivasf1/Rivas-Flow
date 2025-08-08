//User-defined functions
#include "Output.h" 

// OUTPUT DEFINITIONS

Output::Output(){}

//Output::Output(array<double,3>* &f)
  //: field(f) {}

//-----------------------------------------------------------

/*
void Output::PrintResidualNorm(int &cellnum,int &n){

  if (n==1 || n==2 || n==3){ 
    array<double,3> norm = {0.0,0.0,0.0};
    for (int i=0;i<cellnum;i++){

      if (n==1){ //L1 norm case
        norm[0] += abs(field[i][0]); //density
        norm[1] += abs(field[i][1]); //velocity
        norm[2] += abs(field[i][2]); //pressure
      }
      else if (n==2){ //L2 norm case
        norm[0] += pow(field[i][0],2); //density
        norm[1] += pow(field[i][1],2); //velocity
        norm[2] += pow(field[i][2],2); //pressure
      }
      else{ //L inf. case
        if(abs(field[i][0]) > norm[0]) norm[0] = abs(field[i][0]); //density
        if(abs(field[i][1]) > norm[1]) norm[1] = abs(field[i][1]); //velocity
        if(abs(field[i][2]) > norm[2]) norm[2] = abs(field[i][2]); //pressure
      }

    }

    if (n==1){ //L2 norm case continued
      norm[0] = sqrt(norm[0]);
      norm[1] = sqrt(norm[1]);
      norm[2] = sqrt(norm[2]);
    }

    if (n==1 | n==2) //!< Printing out norms
      Tools::print("--L %d Norm Selected\n",n);
    else
      Tools::print("--L infinity Norm Selected\n");

    Tools::print("Density:%f,Velocity:%f,Pressure:%f\n",norm[0],norm[1],norm[2]);

    

  }
  

  else {
    Tools::print("Residual type unknown!\n");
    exit(0); 
  }
}
*/
//-----------------------------------------------------------
void Output::DiscretizationErrorNorms(vector<array<double,3>>* &field,vector<array<double,3>>* &exact_field,vector<array<double,3>>* &errors,SpaceVariables1D* &sols){

  for (int n=0;n<(int)field->size();n++){ //calculating errors
    for (int i=0;i<3;i++)
      (*errors)[n][i] = (*field)[n][i] - (*exact_field)[n][i];
  }
  
  //L2 Norms of Error
  array<double,3> ErrorNorms = sols->ComputeSolutionNorms(errors);

  Tools::print("-------------------------\n");
  Tools::print("Discretization Error Norms\n");
  Tools::print("Density: %e\n",ErrorNorms[0]);
  Tools::print("Velocity: %e\n",ErrorNorms[1]);
  Tools::print("Pressure: %e\n",ErrorNorms[2]);

  return;
}

//-----------------------------------------------------------
void Output::CalculateOrderofAccuracy(const char *filename_read,const char *filename_write){

  ifstream myfileread(filename_read);
  ofstream myfilewrite(filename_write);

  if (!myfileread){ //Error Handling
    cerr<<"Error Opening Discretization Error Norms File!"<<endl;
    return; 
  }

  std::string line;

  vector<double> CellSize;
  vector<double> Density;
  vector<double> Velocity;
  vector<double> Pressure;

  // Reading Discreization Error File(.txt)
  while (std::getline(myfileread,line)){ //reading the line as a string

    std::stringstream ss(line); //reading the line as a string
    std::string label;
    double value;

    if (line.find("Cell Size:") != std::string::npos) { //found Cell Size

      ss >> label >> label >> value; 
      CellSize.push_back(static_cast<int>(value));
    }

    else if (line.find("Density:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Density.push_back(value);
    }

    else if (line.find("Velocity:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Velocity.push_back(value);
    }
    else if (line.find("Pressure:") != std::string::npos) { //found Cell Size
      ss >> label >> value; 
      Pressure.push_back(value);
    }

  }

  // Calculating Observed Order of Accuracy
  //using Section 4 Slide 31 Notes to calc. order of accuracy (p)
  // NOTE: arrangement of PHat lists start from coarsest and go to finest grids

  vector<double> PHat_density((int)CellSize.size()-1,0); //order of accuracy value
  vector<double> PHat_velocity((int)CellSize.size()-1,0); //order of accuracy value
  vector<double> PHat_pressure((int)CellSize.size()-1,0); //order of accuracy value

  double r = 2.0; //mesh refinement factor
  for (int n=0;n<(int)CellSize.size()-1;n++){ 
    PHat_density[n] = (log(Density[n]/Density[n+1])) / log(r);
    PHat_velocity[n] = (log(Velocity[n]/Velocity[n+1])) / log(r);
    PHat_pressure[n] = (log(Pressure[n]/Pressure[n+1])) / log(r);
 }

  vector<double> h; //grid spacing 
  h.push_back(1.0); //1st element is the finest grid
  for (int i=1;i<=(int)CellSize.size()-2;i++) //-2 b/c not evaluating coarsest mesh
    h.push_back(h[i-1]*r); //r times the previous mesh spacing


  reverse(h.begin(),h.end()); //reversing order to match with phat calc.


  // Outputting Observed Order of Accuracy in .dat format
  if (!myfilewrite){ //Error Handling
    cerr<<"Error Opening Output for Observed Order of Accuracy File!"<<endl;
    return; 
  }

  myfilewrite<<"variables= \"grid spacing(h)\" \"Phat(density)\" \"Phat(velocity)\"  \"Phat(Pressure)\""<<endl;

  myfilewrite<<"zone T= "<<"\""<<0<<"\""<<endl;
  myfilewrite<<"I="<<(int)h.size()<<endl;
  myfilewrite<<"DATAPACKING=POINT"<<endl;
  myfilewrite<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;

  for (int n=0;n<(int)h.size();n++)
    myfilewrite<<h[n]<<"  "<<PHat_density[n]<<"  "<<PHat_velocity[n]<<"  "<<PHat_pressure[n]<<"  "<<endl;

  
  myfilewrite.close();

  return;

}

//-----------------------------------------------------------
void Output::OutputPrimitiveVariables(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax){

  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false){ //start of .dat file -- printing initial parameters
    myfile<<"TITLE = \" 2D Field Solutions \""<<endl;
    myfile<<"VARIABLES = \"X\",\"Y\",\"Rho\",\"U\",\"V\",\"P\""<<endl;
    //myfile<<"variables= \"cell index\" \"rho(kg/m^3)\" \"u(m/s)\"  \"Press(N/m^2)\" \"Mach\" \"Xcoords\""<<endl;
  }

//Repeat the following each time you want to write out the solution
/*
write(40,*) 'zone T="',num_iter,'" '
write(40,*) 'I=',imax
write(40,*) 'DATAPACKING=POINT'
write(40,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE
& DOUBLE DOUBLE )'
  */

  myfile<<"ZONE T="<<"\""<<iter<<"\""<<endl; //Now adding zone specific info.
  myfile<<"I="<<imax<<", "<<"J="<<jmax<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;
  //myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  myfile<<"VARLOCATION=([3-6]=CELLCENTERED)"<<endl; //-> tells Tecplot this is cell-centered val (must be size (imax-1)*(jmax-1) size


  // Saving all primitive variables in their own corresponding vector
  vector<double> all_rho,all_u,all_v,all_p;

  for (int i=0;i<cell_number;i++){
    all_rho.push_back((*field)[i][0]);
    all_u.push_back((*field)[i][1]);
    all_v.push_back((*field)[i][2]);
    all_p.push_back((*field)[i][3]);
  }

  int count = 0;
  // Writing Xcoords
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Ycoords
  for (int n=0;n<(int)ycoords.size();n++){
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }


  // Writing Rho
  for (int n=0;n<(int)all_rho.size();n++){
    count++;
    myfile<<std::setw(15)<<all_rho[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  
  // Writing U 
  for (int n=0;n<(int)all_u.size();n++){
    count++;
    myfile<<std::setw(15)<<all_u[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing V 
  for (int n=0;n<(int)all_v.size();n++){
    count++;
    myfile<<std::setw(15)<<all_v[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing P
  for (int n=0;n<(int)all_p.size();n++){
    count++;
    myfile<<std::setw(15)<<all_p[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}

//-----------------------------------------------------------
void Output::OutputManufacturedSourceTerms(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax){

  std::ofstream myfile(filename,(cond==true) ? ios::app : ios::out); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  if (cond==false){ //start of .dat file -- printing initial parameters
    myfile<<"TITLE = \" 2D Field Solutions \""<<endl;
    myfile<<"VARIABLES = \"X\",\"Y\",\"Continuity\",\"X-Momentum\",\"Y-Momentum\",\"Energy\""<<endl;
  }

  myfile<<"ZONE T="<<"\""<<iter<<"\""<<endl; //Now adding zone specific info.
  myfile<<"I="<<imax<<", "<<"J="<<jmax<<endl;
  myfile<<"DATAPACKING=BLOCK"<<endl;
  //myfile<<"DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )"<<endl;
  myfile<<"VARLOCATION=([3-6]=CELLCENTERED)"<<endl; //-> tells Tecplot this is cell-centered val (must be size (imax-1)*(jmax-1) size


  // Saving all primitive variables in their own corresponding vector
  vector<double> cont,xmom,ymom,energy;

  for (int i=0;i<cell_number;i++){
    cont.push_back((*field)[i][0]);
    xmom.push_back((*field)[i][1]);
    ymom.push_back((*field)[i][2]);
    energy.push_back((*field)[i][3]);
  }

  int count = 0;
  // Writing Xcoords
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Ycoords
  for (int n=0;n<(int)ycoords.size();n++){
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }


  // Writing Rho
  for (int n=0;n<(int)cont.size();n++){
    count++;
    myfile<<std::setw(15)<<cont[n];
    if (count % 4 == 0)
      myfile<<endl;
  }
  
  // Writing U 
  for (int n=0;n<(int)xmom.size();n++){
    count++;
    myfile<<std::setw(15)<<xmom[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing V 
  for (int n=0;n<(int)ymom.size();n++){
    count++;
    myfile<<std::setw(15)<<ymom[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing P
  for (int n=0;n<(int)energy.size();n++){
    count++;
    myfile<<std::setw(15)<<energy[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;
}
//-----------------------------------------------------------
void Output::OutputGhostCells(string filename,vector<double> &xcoords,vector<double> &ycoords,int &Nx,int &Ny){

  std::ofstream myfile(filename); //true for append
  //myfile.open(filename);

  if (!myfile){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  //Structured and number of coords. spec.
  myfile<<"1"<<endl;
  myfile<<Nx<<"\t"<<Ny<<endl;

  int count = 0;
  // Writing Xcoords
  for (int n=0;n<(int)xcoords.size();n++){
    count++;
    myfile<<std::setw(15)<<xcoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  // Writing Ycoords
  for (int n=0;n<(int)ycoords.size();n++){
    count++;
    myfile<<std::setw(15)<<ycoords[n];
    if (count % 4 == 0)
      myfile<<endl;
  }

  return;

}
//-----------------------------------------------------------

Output::~Output(){}
