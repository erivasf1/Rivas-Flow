// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"
#include "MeshAccess.hpp"
#include <fstream>
#include <iostream>

class MeshGenBASE {

  public:
  vector<double> xcoords,ycoords; //domain nodes
  int cellnumber;
  int Nx,Ny,Nz;
  int cell_imax,cell_jmax;

  vector<double> right_xcoords,left_xcoords,top_xcoords,btm_xcoords;
  vector<double> right_ycoords,left_ycoords,top_ycoords,btm_ycoords;

  vector<array<double,4>> top_cells,btm_cells,right_cells,left_cells; //ghost cells sub-domains

  vector<array<double,2>> right_cellcenter_coords,left_cellcenter_coords,top_cellcenter_coords,btm_cellcenter_coords;

  MeshGenBASE();

  virtual double GetCellVolume(int cell_id);
  virtual void GenerateMesh(double xmin,double xmax,double ymin,double ymax);
  virtual void ReadMeshFile();
  virtual void OutputMesh();
  virtual void GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id); //-- adds the extra ghost cells to the domain virtual double GetInteriorCellArea(int &i,int &j,int side);
  virtual double GetInteriorCellArea(int i,int j,int side);
  virtual double GetInteriorCellFaceArea(int i,int j,int dir);
  virtual double GetCellVolume(int i,int j);
  virtual array<array<double,4>,2> GetGhostCellCoords(int i,int j, int tag);
  virtual void ComputeGhostCellCenteredCoordinate();
  virtual array<double,2> ComputeOutwardUnitVector(int i,int j,int side);
  virtual array<double,4> GetGhostCellVarVec(int i,int j,int side);
  
  virtual ~MeshGenBASE();

};


class MeshGenNozzle : public MeshGenBASE { //creates a uniform mesh (in x)
  double xmin,xmax;
  double dx;
  

  public:

  MeshGenNozzle(double &a, double &b, int &c);

  double GetCellVolume(int cell_id) override;
  //static double GetCellVolume(int cell_id) override;
 
  void GenMesh(); //TODO: polymorphism here after 2D implementation

  void OutputNozzleAreas(vector<double> &coords,const char *filename);

  ~MeshGenNozzle();


};

class MeshGen2D : public MeshGenBASE { //reads in a non-uniform 2D mesh
 

   //vector<array<double,4>> top_cells[Nx-1],btm_cells[Nx-1],right_cells[Nx-1],left_cells[Nx-1]; //ghost cells sub-domains

  //NOTE: Make indexing of ghost cells consistent with the interior cells for the pertaining indice
  //size of top and btm cells are Nx*2 & size of right and left cells are 2*Ny

  const char* filename;

  public:


  MeshGen2D(const char* name);

  MeshGen2D(double xmin,double xmax,double ymin,double ymax,int Nx,int Ny); //constructor for true Cartesian testing

  void ReadMeshFile() override; //assings values to x and y coords list

  void OutputMesh() override; //generates a true Cartesian mesh. Use only for TESTING!!!

  void GenerateMesh(double xmin,double xmax,double ymin,double ymax) override;

  array<array<double,4>,2> GetCellCoords(int &i,int &j); //fcn. to retrieve coords; indexing: [btm_left,btm_right,top_left,top_right]!
  array<array<double,4>,2> GetGhostCellCoords(int i,int j,int tag) override; 
  //double GetCellVolume(int cell_id) override;
  
  void GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id) override; //main fcn. that generates ghost nodes to each pertaining bounds
  void ReflectGhostCoords(int tag); //used for slip wall treatments
  void ExtendGhostCoords(int tag); //used for inflow and outflow treatments
  double GetInteriorCellArea(int i,int j,int side) override; //retrieves the area of the specified side of an interior cell (really just a length); side =0(top),1(btm),2(left),3(right)
  double GetInteriorCellFaceArea(int i,int j,int dir) override; //retrieves the cell internal area given a specified dir. in computational coordinates (dir:0=i or 1=j)
  double GetGhostCellArea(int i,int j,int side); //side refers to the boundary of the domain (instead of the side of the cell)
  void ComputeGhostCellCenteredCoordinate() override; //computes+stores the ghost cell center coordinates 
  double GetCellVolume(int i,int j) override;
  array<double,2> ComputeOutwardUnitVector(int i,int j,int side) override;
  array<double,4> GetGhostCellVarVec(int i,int j,int side) override; //extracts vector of primitive vars. from the specified ghost cell
  void AssignGhostCellVarVec(int i,int j,int side,array<double,4> &res);

  ~MeshGen2D();

};

#endif
