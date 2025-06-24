//User-defined functions
#include "MeshGen.h" 

// MESHGENBASE DEFINITIONS
//-----------------------------------------------------------
MeshGenBASE::MeshGenBASE(){}
//-----------------------------------------------------------
double MeshGenBASE::GetCellVolume(int cell_id){
  return (double)cell_id * 0.0; //return 0.0 by default
}

//-----------------------------------------------------------
void MeshGenBASE::GenerateMesh(){}

//-----------------------------------------------------------
void MeshGenBASE::ReadMeshFile(){}

//-----------------------------------------------------------
void MeshGenBASE::OutputMesh(){}

//-----------------------------------------------------------
void MeshGenBASE::GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id){}
//-----------------------------------------------------------
MeshGenBASE::~MeshGenBASE(){}
//-----------------------------------------------------------

// MESHGEN1D DEFINITIONS

MeshGenNozzle::MeshGenNozzle(double &a,double &b,int &c)
  : xmin(a), xmax(b) {

  cellnumber = c;
}

//-----------------------------------------------------------
double MeshGenNozzle::GetCellVolume(int cell_id){

  //using trapezoidal rule since areas are stored at cell faces
  double area_leftface = Tools::AreaVal(xcoords[cell_id]); 
  double area_rightface = Tools::AreaVal(xcoords[cell_id+1]); 
  double area_characteristic = 0.5*(area_leftface+area_rightface); //avg. of left and right face cell area
  //double DArea = (dx/2.0) * abs(area_rightface - area_leftface);
  //previous:double DArea = (dx/2.0) * abs(area_rightface - area_leftface);

  double vol = area_characteristic * dx;
  return vol;

}

//-----------------------------------------------------------
void MeshGenNozzle::GenerateMesh() {

  int facenum = cellnumber + 1; //number of faces
  xcoords = Tools::RetrievePoints(xmin,xmax,facenum);  

  return;

}
//-----------------------------------------------------------
void MeshGenNozzle::OutputNozzleAreas(vector<double> &coords,const char *filename){

  //Computing Areas
  vector<double> Areas((int)coords.size());
  for (int n=0;n<(int)coords.size();n++){
    Areas[n] = Tools::AreaVal(coords[n]);
  }


  //Printing out in filename
  ofstream myfile;
  myfile.open(filename);
  //ofstream myfile(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }
    myfile<<"Areas"<<endl;
    myfile<<"XCoord"<<"  "<<"Area"<<endl;
    for (int i=0;i<(int)coords.size();i++){
  
      myfile<<coords[i]<<"  "<<Areas[i]<<endl;

    }

  
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;


}

//-----------------------------------------------------------
MeshGenNozzle::~MeshGenNozzle(){}

//-----------------------------------------------------------
MeshGen2D::MeshGen2D(const char* name) 
  : filename(name) {

  ReadMeshFile(); //extract nodal coords. of mesh
}
  
//-----------------------------------------------------------
void MeshGen2D::ReadMeshFile(){

  ifstream myfileread(filename);
  
  if (!myfileread) { //Error Handling
    cerr<<"No Mesh File Provided!"<<endl;
    //cerr<<"Error Opening Mesh File "<<filename<<" !"<<endl; 
    return;
  }

  std::string line;
  std::getline(myfileread,line); //!< skips 1st line

  if (std::getline(myfileread,line)){ //!< calling 2nd line
    std::istringstream iss(line); //converts string into stream
    iss >> Nx >> Ny >> Nz; //setting values of 2nd line to imax,jmax,kmax, respectively
  }
  else{
    cerr<<"Mesh File does not contain a second line!"<<endl; 
    return;
  }

  double val; //value of i,j,&k indices
  int total_pts = Nx*Ny*Nz;

  vector<double> xcoords_orig,ycoords_orig; //needed for duplicated xcoords and ycoords when k>1 (i.e. 3D structured grid)
  
  //int pt_ct = 1;
  
  for (int j=0;j<2;j++){
    for (int i=0;i<total_pts;i++){
      myfileread >> val;
      if (j==0)
        xcoords_orig.push_back(val);  
      //else if (j==1)
        //ycoords_orig.push_back(val);
      else
        ycoords_orig.push_back(val);
    }
  }

  //Case if 3D structured mesh is provided -- to extract the x & y coords in 2D plane
  if (Nz>0){
    for(int n=0;n<(int)Nx*Ny;n++){
      xcoords.push_back(xcoords_orig[n]);
      ycoords.push_back(ycoords_orig[n]);
    }
  }
    
  //skipping first line for now
  //2nd line: 1st int refers to imax and 2nd int refers to jmax

  cellnumber = (Nx-1) * (Ny-1); //1 more faces than each dir.

}
//-----------------------------------------------------------
void MeshGen2D::OutputMesh(){

  ifstream myfileread(filename); 
  ofstream myfilewrite("Mesh.dat",ios::out);

  if (!myfilewrite){
    cerr<<"Error: Could not Open \"Mesh.dat File\"!"<<endl;
    return;
  }

  //Writing header to Mesh Output File
  myfilewrite<<"TITLE = \" 2D Structured Mesh \""<<endl;
  myfilewrite<<"VARIABLES = \"X\",\"Y\",\"Z\""<<endl;
  myfilewrite<<"ZONE I="<<Nx<<", J="<<Ny<<", K="<<Nz<<" DATAPACKING=BLOCK"<<endl;
  
  
  //Reading vals. of input mesh file & writing to output mesh file
  double val;
  string line;
  std::getline(myfileread,line);//!< skipping 1st and 2nd line
  std::getline(myfileread,line);

  //std::istringstream iss(line);
  int count = 0;
  while(myfileread>>val){
    count++;
    myfilewrite<<std::setw(15)<<val;
    if (count % 4 == 0)
      myfilewrite<<endl;
 
  }
    
  //myfilewrite.close(); 

}

//-----------------------------------------------------------
void MeshGen2D::GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id){

  //Left Boundary
  //Left Boundary
  //Bottom Boundary
  //(btm_id == 0) ? ReflectGhostCoords(0) : ExtendGhostCoords();
  ReflectGhostCoords(0); //temp
  //Left Boundary
  
  return;

}

//-----------------------------------------------------------
void MeshGen2D::ReflectGhostCoords(int tag){

  //using the symmetry method: https://doi-org.ezproxy.lib.vt.edu/10.2514/3.11983
  //symmetry imposed by creating ghost cell nodes as same sign in x but different sign in y (for top + btm cases), vice versa in right and left cases
  
  //Note: there are 2 layers of ghost nodes that comprise the ghost cell domains
  double x,x1,x2,x3; //x&y: reflected pt. && x1&y1: 1st interior pt. && x2&y2: 2nd interior pt. && x3&y3: interior pt. used to determine reflected axis
  double y,y1,y2,y3;
  int pt_id1,pt_id2,pt_id3; //pt ids in loop
  double p1x,p1y;
  double p2x,p2y;
  double d1,d2;
  double pdotp; //dot products of point vectors
  double theta;
  int ghost_id; 

  //btm boundary case (tag=0) 
  if (tag == 0){
    for (int j=0;j<2;j++){
      for (int i=0;i<Nx;i++){
        //Step 1: determining pts + retrieving coords 
        pt_id1 = i + (j*Nx);
        pt_id2 = i + ((j+1)*Nx);//ptid2 is j+1 from ptid1
        pt_id3 = (i<Nx-1) ? (i+1) + (j*Nx) : (i-1) + (j*Nx); //use i-1 pt. for last i node
        x1 = xcoords[pt_id1]; y1 = ycoords[pt_id1]; //1st interior pt
        x2 = xcoords[pt_id2]; y2 = ycoords[pt_id2]; //2nd interior pt
        x3 = xcoords[pt_id3]; y3 = ycoords[pt_id3]; //3rd interior pt
        //Step 2: computing dist. vectors + mag.(axis of reflection and interior pt line (1&2)
        p1x = x2-x1; p1y = y2-y1;
        p2x = x3-x1; p2y = y3-y1; 
        d1 = sqrt(pow(p1x,2.0) + pow(p1y,2.0));//mag. of interior pt. line (L2 norm)
        d2 = sqrt(pow(p2x,2.0) + pow(p2y,2.0));//mag. of reflection axis line (L2 norm)
        //Step 3: computing angle
        pdotp = Tools::ComputeDotProduct(p1x,p1y,p2x,p2y);
        theta = acos(pdotp / (d1*d2)); 
        //Step 4: computing coords of reflected pt.
        if (j>0) //for 2nd layer, the last ghost node layer are computed from the previously computed layer
          ghost_id = i + (j-1)*Nx;
        if (i==Nx-1){ //for nodes that are at the imax loc.
          x = (j<1) ? d1*cos(-theta) + x1 : d1*cos(-theta) + btm_xcoords[ghost_id]; 
          y = (j<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + btm_ycoords[ghost_id];
        }
        else {
          x = (j<1) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + btm_xcoords[ghost_id]; 
          y = (j<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + btm_ycoords[ghost_id];
        }
        //x = (j<1) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + btm_xcoords[ghost_id]; 
        //y = (j<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + btm_ycoords[ghost_id];
        btm_xcoords.push_back(x);
        btm_ycoords.push_back(y);
        
      }
    }

  }
  //TODO:right boundary case -- tag=1 case
  //TODO:top boundary case -- tag=2 case
  //TODO:btm boundary case -- tag=3 case

  return;

}

//-----------------------------------------------------------
MeshGen2D::~MeshGen2D(){}
//-----------------------------------------------------------
