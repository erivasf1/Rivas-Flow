// Responsible for creating a Mesh (1D for now)
#ifndef _MESHGEN_H_
#define _MESHGEN_H_
#include "ExactNozzle.h"
#include <fstream>
#include <iostream>

class MeshGenBASE {

  public:
  vector<double> xcoords,ycoords; //domain nodes
  //vector top_xcoords,top_ycoords; //top ghost cell nodes
  //vector btm_xcoords,btm_ycoords; //btm ghost cell nodes
  //vector right_xcoords,right_ycoords; //right ghost cell nodes
  //vector left_xcoords,left_ycoords; //left ghost cell nodes
  int cellnumber;
  int Nx,Ny,Nz;

  MeshGenBASE();

  virtual double GetCellVolume(int cell_id);
  virtual void GenerateMesh();
  virtual void ReadMeshFile();
  virtual void OutputMesh();
  virtual void GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id); //-- adds the extra ghost cells to the domain 
  
  virtual ~MeshGenBASE();

};


class MeshGenNozzle : public MeshGenBASE { //creates a uniform mesh (in x)
  double xmin,xmax;
  double dx;
  

  public:

  MeshGenNozzle(double &a, double &b, int &c);

  double GetCellVolume(int cell_id) override;
  //static double GetCellVolume(int cell_id) override;
 
  void GenerateMesh() override;

  void OutputNozzleAreas(vector<double> &coords,const char *filename);

  ~MeshGenNozzle();


};

class MeshGen2D : public MeshGenBASE { //reads in a non-uniform 2D mesh
 
  vector<double> right_xcoords,left_xcoords,top_xcoords,btm_xcoords;
  vector<double> right_ycoords,left_ycoords,top_ycoords,btm_ycoords;

  //vector<array<double,4>> top_cells,btm_cells,right_cells,left_cells;

  //NOTE: Make indexing of ghost cells consistent with the interior cells for the pertaining indice
  //size of top and btm cells are Nx*2 & size of right and left cells are 2*Ny

  const char* filename;

  public:
  MeshGen2D(const char* name);

  void ReadMeshFile() override;

  void OutputMesh();

  //array<double,2> GetCoords(int ghost_domain,int cell_id); //fcn. to retrieve coords
  //double GetCellVolume(int cell_id) override;
  
  void GenerateGhostCells(int left_id,int right_id,int btm_id,int top_id) override; //main fcn. that generates ghost nodes to each pertaining bounds
  void ReflectGhostCoords(int tag); //used for slip wall treatments
  void ExtendGhostCoords(int tag); //used for inflow and outflow treatments
  //double GetArea(int &cell_id,int &side); //retrives the area of the specified side of a cell
  //void Farfield

  ~MeshGen2D();

};

#endif
