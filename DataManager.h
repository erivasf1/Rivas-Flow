// Responsible for creating a Mesh (1D for now)
#ifndef _DATAMANAGER_H_
#define _DATAMANAGER_H_
//#include "MeshGen.h"
#include "EulerOperator.h"

//Forward Declarations
class Euler1D;
class Euler2D;
class Euler2DMMS;

//SpaceVariable Class stores the values of ALL flow quantities at all points

//SpaceVariable Class for 1D Problems
class SpaceVariables1D {
  int cell_number;

  public:

  SpaceVariables1D();

  array<double,3> ComputeSolutionNorms(vector<array<double,3>>* &resid);

  double ComputeNormAvg(array<double,3> &Norms);

  double ComputeRampValue(array<double,3> CurrentNorms,array<double,3> InitNorms,double FinalVal); //outputs a ramping function val.

  void OutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,const char *filename);

  void AllOutputPrimitiveVariables(vector<array<double,3>>* &field,Euler1D* &euler,string filename,bool cond,int iter,vector<double> &xcoords);

  void OutputLocalResiduals(vector<array<double,3>> &Resid,const char *filename); 

  void OutputResidualTerms(array<double,3> F_right,array<double,3> F_left,double S,array<double,3> D2_right,array<double,3> D2_left,array<double,3> D4_right,array<double,3> D4_left,const char* filename); //Not needed for now

  void ComputeCellAveragedSol(vector<array<double,3>>* &cell_faces,vector<array<double,3>>* &cell_sols,vector<double> &xcoords); //area evaluated for quasi-steady 1D nozzle case

  ~SpaceVariables1D();

};


//SpaceVariable Class for 2D Problems
class SpaceVariables2D { //TODO: Will replace SpaceVariables1D

  public:

  SpaceVariables2D();
  
  void ComputeInteriorCellCenteredCoordinate(vector<double> &xcoords,vector<double> &ycoords,vector<double> &cell_center_xcoords,vector<double> &cell_center_ycoords,int imax); //approximates cell-center avg. by averaging corner nodes of cell
  void ComputeGhostCellCenteredCoordinate(MeshGen2D* &mesh); //TODO:approximates cell-center avg. by averaging corner nodes of cell

  ~SpaceVariables2D();

};

#endif
