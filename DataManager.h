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
  void ComputeGhostCellCenteredCoordinate(MeshGen2D* &mesh); //For Interior Cells Only
  array<double,4> ComputeSolutionNorms(vector<array<double,4>>* &resid);
  double ComputeNormAvg(array<double,4> &norms);

  ~SpaceVariables2D();

};

//Vector class for 2D Problems
class Vector {

  double v[4];

  public:
  Vector() {v[0] = 0.0;v[1] = 0.0;v[2] = 0.0;v[3] = 0.0;}
  Vector(double d1,double d2,double d3,double d4) {v[0] = d1; v[1] = d2; v[2] = d3; v[3] = d4;}

  Vector operator+(const Vector& vec2) {Vector res; for(int i=0;i<4;i++) res.v[i] = v[i] + vec2.v[i]; return res;}

  Vector operator-(const Vector& vec2) {Vector res; for(int i=0;i<4;i++) res.v[i] = v[i] - vec2.v[i]; return res;}

  Vector operator*(const Vector& vec2) {Vector res; for(int i=0;i<4;i++) res.v[i] = v[i] * vec2.v[i]; return res;}

  Vector operator*(double& val)const {Vector res; for(int i=0;i<4;i++) res.v[i] = v[i] * val; return res;}

  Vector operator*=(double val)const {Vector res;for (int i=0;i<4;i++) res.v[i] *= val; return res;}

  Vector operator/(const Vector& vec2) {Vector res; for(int i=0;i<4;i++) res.v[i] = v[i] / vec2.v[i]; return res;}

  double operator[](int i) const {return v[i];}
  double &operator[](int i) {return v[i];}

  ~Vector(){};

};

// scalar * vector (non-member fcn.) -- to allow "scalar * Vector"
inline Vector operator*(double val, const Vector& vec) {
    return vec * val;
}


#endif
