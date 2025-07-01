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
void MeshGenBASE::GenerateGhostCells(int ,int ,int ,int ){}
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
array<array<double,4>,2> MeshGen2D::GetCellCoords(int &i,int &j){

  //int cell_imax = Nx-2;
  double x1,x2,x3,x4; 
  double y1,y2,y3,y4; 
  int cell_id = i + (j*Nx);
  int pt_id1,pt_id2,pt_id3,pt_id4;
  
  array<array<double,4>,2> cell_coords;
    
  pt_id1 = (j==0) ? cell_id : cell_id+j; //for j>0, cell nodes and pt nodes are offsetted by the # of j rows above j=0

  pt_id2 = cell_id+1; pt_id3 = cell_id + Nx;
  pt_id4 = cell_id + (Nx+1);

  x1 = xcoords[pt_id1]; y1 = ycoords[pt_id1];
  x2 = xcoords[pt_id2]; y2 = ycoords[pt_id2];
  x3 = xcoords[pt_id3]; y3 = ycoords[pt_id3];
  x4 = xcoords[pt_id4]; y4 = ycoords[pt_id4];
  
  //saving cell's xcoords to list
  cell_coords[0][0] = x1;cell_coords[0][1] = x2;
  cell_coords[0][2] = x3;cell_coords[0][3] = x4;
  //saving cell's ycoords to list
  cell_coords[1][0] = y1;cell_coords[1][1] = y2;
  cell_coords[1][2] = y3;cell_coords[1][3] = y4;

  return cell_coords; 
}

//-----------------------------------------------------------
void MeshGen2D::GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id){

  //Bottom Boundary
  (btm_id == 1) ? ReflectGhostCoords(0) : ExtendGhostCoords(0);
  //Top Boundary
  (top_id == 1) ? ReflectGhostCoords(1) : ExtendGhostCoords(1);
  //Left Boundary
  (left_id == 1) ? ReflectGhostCoords(2) : ExtendGhostCoords(2);
  //Right Boundary
  (right_id == 1) ? ReflectGhostCoords(3) : ExtendGhostCoords(3);
  
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

  //btm boundary case (tag==0)
  if (tag == 0) {

    for (int j=0;j<2;j++){ //j increasing
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
        btm_xcoords.push_back(x);
        btm_ycoords.push_back(y);
        
      }
    }

  }
  //TODO:top boundary case (tag==1)
  else if (tag==1){

    for (int j=Ny-1;j>Ny-3;j--){ //j decreasing 
      for (int i=0;i<Nx;i++){
        //Step 1: determining pts + retrieving coords 
        pt_id1 = i + (j*Nx);
        pt_id2 = i + ((j-1)*Nx);//ptid2 is j-1 from ptid1
        pt_id3 = (i<Nx-1) ? (i+1) + (j*Nx) : (i-1) + (j*Nx); //use i-1 pt. for last i node
        x1 = xcoords[pt_id1]; y1 = ycoords[pt_id1]; //1st interior pt
        x2 = xcoords[pt_id2]; y2 = ycoords[pt_id2]; //2nd interior pt
        x3 = xcoords[pt_id3]; y3 = ycoords[pt_id3]; //3rd interior pt
        //Step 2: computing dist. vectors + mag.(axis of reflection and interior pt line (1&2))
        p1x = x2-x1; p1y = y2-y1;
        p2x = x3-x1; p2y = y3-y1; 
        d1 = sqrt(pow(p1x,2.0) + pow(p1y,2.0));//mag. of interior pt. line (L2 norm)
        d2 = sqrt(pow(p2x,2.0) + pow(p2y,2.0));//mag. of reflection axis line (L2 norm)
        //Step 3: computing angle
        pdotp = Tools::ComputeDotProduct(p1x,p1y,p2x,p2y);
        theta = acos(pdotp / (d1*d2)); 
        //Step 4: computing coords of reflected pt.
        if (j<Ny-1) //for 2nd layer, the last ghost node layer are computed from the previously computed layer
          ghost_id = i;
          //ghost_id = i + (j-1)*Nx;
        if (i==Nx-1){ //for nodes that are at the imax loc.
          x = (j>Ny-2) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + top_xcoords[ghost_id]; 
          y = (j>Ny-2) ? d1*sin(-theta) + y1 : d1*sin(-theta) + top_ycoords[ghost_id];
        }
        else { //Nodes less than imax
          x = (j>Ny-2) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + top_xcoords[ghost_id]; 
          y = (j>Ny-2) ? d1*sin(-theta) + y1 : d1*sin(-theta) + top_ycoords[ghost_id];
        }
        top_xcoords.push_back(x);
        top_ycoords.push_back(y);

        }
      }
  }
  //TODO:left boundary case (tag==2)
  else if (tag == 2) {

    for (int i=0;i<2;i++){ //i increasing 
      for (int j=0;j<Ny;j++){
        //Step 1: determining pts + retrieving coords 
        pt_id1 = i + (j*Nx);
        pt_id2 = (i+1) + ((j)*Nx);//ptid2 is j+1 from ptid1
        pt_id3 = (j<Ny-1) ? i + ((j+1)*Nx) : i + ((j-1)*Nx); //use i-1 pt. for last i node
        x1 = xcoords[pt_id1]; y1 = ycoords[pt_id1]; //1st interior pt
        x2 = xcoords[pt_id2]; y2 = ycoords[pt_id2]; //2nd interior pt
        x3 = xcoords[pt_id3]; y3 = ycoords[pt_id3]; //3rd interior pt
        //Step 2: computing dist. vectors + mag.(axis of reflection and interior pt line (1&2))
        p1x = x2-x1; p1y = y2-y1;
        p2x = x3-x1; p2y = y3-y1; 
        d1 = sqrt(pow(p1x,2.0) + pow(p1y,2.0));//mag. of interior pt. line (L2 norm)
        d2 = sqrt(pow(p2x,2.0) + pow(p2y,2.0));//mag. of reflection axis line (L2 norm)
        //Step 3: computing angle
        pdotp = Tools::ComputeDotProduct(p1x,p1y,p2x,p2y);
        theta = acos(pdotp / (d1*d2)); 
        //Step 4: computing coords of reflected pt.
        if (i>0) //for 2nd layer, the last ghost node layer are computed from the previously computed layer
          ghost_id = (i-1) + (j*Nx);
        if (j==Ny-1){ //for nodes that are at the jmax loc.
          x = (i<1) ? d1*cos(-theta) + x1 : d1*cos(-theta) + left_xcoords[ghost_id]; 
          y = (i<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + left_ycoords[ghost_id];
        }
        else {
          x = (i<1) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + left_xcoords[ghost_id]; 
          y = (i<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + left_ycoords[ghost_id];
        }
        left_xcoords.push_back(x);
        left_ycoords.push_back(y);

        }
      }


  }
  //TODO:right boundary case (tag==3)
  else if (tag==3) {

    for (int i=Nx-1;i>Nx-3;i--){ //i decreasing 
      for (int j=0;j<Ny;j++){
        //Step 1: determining pts + retrieving coords 
        pt_id1 = i + (j*Nx);
        pt_id2 = (i-1) + (j*Nx);//ptid2 is i-1 from ptid1
        pt_id3 = (j<Ny-1) ? i + ((j+1)*Nx) : i + ((j-1)*Nx); //use i-1 pt. for last i node
        x1 = xcoords[pt_id1]; y1 = ycoords[pt_id1]; //1st interior pt
        x2 = xcoords[pt_id2]; y2 = ycoords[pt_id2]; //2nd interior pt
        x3 = xcoords[pt_id3]; y3 = ycoords[pt_id3]; //3rd interior pt
        //Step 2: computing dist. vectors + mag.(axis of reflection and interior pt line (1&2))
        p1x = x2-x1; p1y = y2-y1;
        p2x = x3-x1; p2y = y3-y1; 
        d1 = sqrt(pow(p1x,2.0) + pow(p1y,2.0));//mag. of interior pt. line (L2 norm)
        d2 = sqrt(pow(p2x,2.0) + pow(p2y,2.0));//mag. of reflection axis line (L2 norm)
        //Step 3: computing angle
        pdotp = Tools::ComputeDotProduct(p1x,p1y,p2x,p2y);
        theta = acos(pdotp / (d1*d2)); 
        //Step 4: computing coords of reflected pt.
        if (i<Nx-1) //for 2nd layer, the last ghost node layer are computed from the previously computed layer
          ghost_id = j;
        if (j==Ny-1){ //for nodes that are at the jmax loc.
          x = (i<1) ? d1*cos(-theta) + x1 : d1*cos(-theta) + right_xcoords[ghost_id]; 
          y = (i<1) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + right_ycoords[ghost_id];
        }
        else { //nodes that are from 0<j<jmax-1
          x = (i>Nx-2) ? -d1*cos(-theta) + x1 : -d1*cos(-theta) + right_xcoords[ghost_id]; 
          y = (i>Nx-2) ? -d1*sin(-theta) + y1 : -d1*sin(-theta) + right_ycoords[ghost_id];
        }
        right_xcoords.push_back(x);
        right_ycoords.push_back(y);

        }
      }
      



  }
    //error handling!
    else {
      cerr<<"Unidentified tag number!"<<endl;  
      return;
    }

  return;

}

//-----------------------------------------------------------
void MeshGen2D::ExtendGhostCoords(int tag){

  //algoritm: find x&y diff. between corresponding interior point and boundary and then subtract that diff. to boundary point to find ghost coord
  double x,x1,x2; //x&y: reflected pt. && x1&y1: 1st interior pt. && x2&y2: 2nd interior pt. && x3&y3: interior pt. used to determine reflected axis
  double y,y1,y2;
  int pt_id1,pt_id2; //pt ids in loop
  double dx,dy;

  //btm boundary
  if (tag==0){
    for (int j=0;j<2;j++){
      for (int i=0;i<Nx;i++){
        //extracting pts
        pt_id1 = i + (j*Nx); pt_id2 = i + (j+1)*Nx; 
        x1 = xcoords[pt_id1];x2 = xcoords[pt_id2];
        y1 = ycoords[pt_id1];y2 = ycoords[pt_id2];
        //computing x+y diffs.
        dx = x2-x1;
        dy = y2-y1;
        //computing ghost cell coord.
        x = x1 - dx; y = y1 - dy;
        //saving ghost cell coords to list
        btm_xcoords.push_back(x); 
        btm_ycoords.push_back(y); 
      }
    }
  }

  //top boundary
  if (tag==1){
    for (int j=Ny-1;j<Ny-3;j--){
      for (int i=0;i<Nx;i++){
        //extracting pts
        pt_id1 = i + (j*Nx); pt_id2 = i + (j-1)*Nx; 
        x1 = xcoords[pt_id1];x2 = xcoords[pt_id2];
        y1 = ycoords[pt_id1];y2 = ycoords[pt_id2];
        //computing x+y diffs.
        dx = x2-x1;
        dy = y2-y1;
        //computing ghost cell coord.
        x = x1 - dx; y = y1 - dy;
        //saving ghost cell coords to list
        top_xcoords.push_back(x); 
        top_ycoords.push_back(y); 
      }
    }
  }

  //left boundary
  if (tag==2){
    for (int i=0;i<2;i++){
      for (int j=0;j<Ny;j++){
        //extracting pts
        pt_id1 = i + (j*Nx); pt_id2 = (i+1) + (j*Nx); 
        x1 = xcoords[pt_id1];x2 = xcoords[pt_id2];
        y1 = ycoords[pt_id1];y2 = ycoords[pt_id2];
        //computing x+y diffs.
        dx = x2-x1;
        dy = y2-y1;
        //computing ghost cell coord.
        x = x1 - dx; y = y1 - dy;
        //saving ghost cell coords to list
        left_xcoords.push_back(x); 
        left_ycoords.push_back(y); 
      }
    }
  }

  //right boundary
  if (tag==3){
    for (int i=Nx-1;i<Nx-3;i--){
      for (int j=0;j<Ny;j++){
        //extracting pts
        pt_id1 = i + (j*Nx); pt_id2 = (i-1) + (j*Nx); 
        x1 = xcoords[pt_id1];x2 = xcoords[pt_id2];
        y1 = ycoords[pt_id1];y2 = ycoords[pt_id2];
        //computing x+y diffs.
        dx = x2-x1;
        dy = y2-y1;
        //computing ghost cell coord.
        x = x1 - dx; y = y1 - dy;
        //saving ghost cell coords to list
        right_xcoords.push_back(x); 
        right_ycoords.push_back(y); 
      }
    }
  }


  return;
}

//-----------------------------------------------------------
double MeshGen2D::GetInteriorCellArea(int &i,int &j,int side){
  //side: top = 0, btm = 1, left = 2, right = 3
  double node1_x,node1_y,node2_x,node2_y;
  int node1_id,node2_id;
  double area;
  if (side == 0){ //top side case
    node1_id = i + ((j+1)*Nx);  
    node2_id = (i+1) + ((j+1)*Nx);  
    node1_x = xcoords[node1_id]; node1_y = ycoords[node1_id];
    node2_x = xcoords[node2_id]; node2_y = ycoords[node2_id];
    
    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));
  }

  else if (side == 1){ //btm side case
    node1_id = i + (j*Nx);  
    node2_id = (i+1) + (j*Nx);  
    node1_x = xcoords[node1_id]; node1_y = ycoords[node1_id];
    node2_x = xcoords[node2_id]; node2_y = ycoords[node2_id];
    
    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));
  }

  else if (side == 2){ //left side case
    node1_id = i + (j*Nx);  
    node2_id = i + ((j+1)*Nx);  
    node1_x = xcoords[node1_id]; node1_y = ycoords[node1_id];
    node2_x = xcoords[node2_id]; node2_y = ycoords[node2_id];
    
    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));
  }

  else if (side == 3){ //right side case
    node1_id = (i+1) + (j*Nx);  
    node2_id = (i+1) + ((j+1)*Nx);  
    node1_x = xcoords[node1_id]; node1_y = ycoords[node1_id];
    node2_x = xcoords[node2_id]; node2_y = ycoords[node2_id];
    
    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));
  }

  else { //error handling
    cerr<<"Error:Unkown side # when retrieving area fcn."<<endl;

  }

  return area;
}

//-----------------------------------------------------------
double MeshGen2D::GetGhostCellArea(int &i,int &j,int side){ //retrieves the area of the specified side of the domain

  //side: top = 0, btm = 1, left = 2, right = 3
  double node1_x,node1_y,node2_x,node2_y;
  int node1_id,node2_id;
  double area;

  node1_id = i + (j*Nx); node2_id = (i+1) + (j*Nx);
  if (side == 0){ //top side case
    node1_x = top_xcoords[node1_id]; node1_y = top_xcoords[node1_id];
    node2_x = top_xcoords[node2_id]; node2_y = top_xcoords[node2_id];

    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));

  } 

  else if (side == 1){ //btm side case
    node1_x = btm_xcoords[node1_id]; node1_y = btm_xcoords[node1_id];
    node2_x = btm_xcoords[node2_id]; node2_y = btm_xcoords[node2_id];

    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));

  }

  else if (side == 2){ //left side case
    node1_x = left_xcoords[node1_id]; node1_y = left_xcoords[node1_id];
    node2_x = left_xcoords[node2_id]; node2_y = left_xcoords[node2_id];

    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));
  }

  else if (side == 3){ //right side case
    node1_x = right_xcoords[node1_id]; node1_y = right_xcoords[node1_id];
    node2_x = right_xcoords[node2_id]; node2_y = right_xcoords[node2_id];

    area = sqrt(pow(node2_x-node1_x,2.0) + pow(node2_y-node1_y,2.0));

  }

  else { //error handling
    cerr<<"Error:Unkown side # when retrieving area fcn."<<endl;
  }

  return area;
}

//-----------------------------------------------------------
array<double,2> MeshGen2D::ComputeOutwardUnitVector(int &i,int &j,int side){

  //side refering to the side of the cell
  double nx,ny; //unit normal vector components
  double x1,y1,x2,y2;
  double area;
  array<array<double,4>,2> cell_coords = GetCellCoords(i,j);

  if (side == 0){ //top side case
    x1 = cell_coords[0][2]; x2 = cell_coords[0][3];
    y1 = cell_coords[1][2]; y2 = cell_coords[1][3];
    area = GetInteriorCellArea(i,j,0);
    
    nx = y2-y1 / area; ny = x2-x1 / area;
    
  }

  else if (side == 1){ //btm side case
    x1 = cell_coords[0][0]; x2 = cell_coords[0][1];
    y1 = cell_coords[1][0]; y2 = cell_coords[1][1];
    area = GetInteriorCellArea(i,j,1);
    
    nx = y2-y1 / area; ny = x2-x1 / area;

  }

  else if (side == 2){ //left side case
    x1 = cell_coords[0][0]; x2 = cell_coords[0][2];
    y1 = cell_coords[1][0]; y2 = cell_coords[1][2];
    area = GetInteriorCellArea(i,j,2);
    
    nx = y2-y1 / area; ny = x2-x1 / area;

  }

  else if (side == 3){ //right side case
    x1 = cell_coords[0][1]; x2 = cell_coords[0][3];
    y1 = cell_coords[1][1]; y2 = cell_coords[1][3];
    area = GetInteriorCellArea(i,j,3);
    
    nx = y2-y1 / area; ny = x2-x1 / area;

  }

  else {
    cerr<<"Unknown side # specification!"<<endl;
  }
  
  array<double,2> normals{nx,ny};
  return normals;

}

//-----------------------------------------------------------
MeshGen2D::~MeshGen2D(){}
//-----------------------------------------------------------
