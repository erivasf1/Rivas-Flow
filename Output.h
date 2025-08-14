// Responsible for outputting results in a clean format
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <fstream> 
#include <cmath> 

#include "ExactNozzle.h"
#include "MeshGen.h"
#include "DataManager.h"
//#include "TimeIntegrator.h"

using namespace std;

//Forward Declarations
class SpaceVariables1D;
class SpaceVariables2D;

class Output {
//  array<double,3>* &field;

  public:
  //Output(array<double,3>* &field);
  Output();
  
  void PrintResidualNorm(int &cellnum,int &n);

  void DiscretizationErrorNorms(vector<array<double,3>>* &field,vector<array<double,3>>* &exact_field,vector<array<double,3>>* &errors,SpaceVariables1D* &sols);

  void CalculateOrderofAccuracy(const char *filename_read,const char *filename_write); //creates a new file containing the order of accuracy value given the discretization error file.txt
 
  void OutputPrimitiveVariables(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax); //outputs primitive variables in tecplot format

  void OutputManufacturedSourceTerms(vector<array<double,4>>* &field,string filename,bool cond,int iter,vector<double> &xcoords,vector<double> &ycoords,int cell_number,int imax,int jmax); //outputs primitive variables in tecplot format

  void OutputGhostCoords(string filename,vector<double> &xcoords,vector<double> &ycoords,int Nx,int Ny); //for visualizing the ghost cells
  void OutputGhostCells(vector<array<double,4>>* &ghost_cell,string filename,vector<double> &xcoords,vector<double> &ycoords,vector<double> &ghost_xcoords,vector<double> &ghost_ycoords,int Nx,int Ny,int ghost_Nx,int ghost_Ny,int side); //for visualizing the ghost cells

  //void ConvertToDatFile(const char*filename_read,const char *filename_write); //TODO: creates a .dat file of a given .txt file

  ~Output();


};




#endif
