//User-defined functions
#ifndef _EXACTNOZZLE_H_
#define _EXACTNOZZLE_H_
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdarg.h>
#include <array>
#include <fstream>

using namespace std;

class Tools{
  public:

    static void print(const char format[],...);

    static void printWithPrecision(int precision,const char* format,...);

    static vector<double> RetrievePoints(double xmin,double xmax,int pt_num);

    static double AreaVal(double x); // Returns the area val. at a specified x-location
    static double ComputeDotProduct(double &x1,double &y1,double &x2,double &y2);
  
};

class SuperSonicNozzle : public Tools{

  double area,area_star,stag_pressure,stag_temp;
  bool cond;

  constexpr static double gamma = 1.4; // specific heat ratio
  double Ru = 8314.0; // J/(kmol*K) -- universal gas constant
  double MolMass = 28.96; // kg/kmol

  double tol = 1.0e-14; //tolerance for Mach # computation
  double maxiter = 1e5;

  public:

  SuperSonicNozzle(double &a,double &b,double &c,double &d,bool &e);
 
  static double GetPhi(double M); //Used for Mach Number computation
  static double GetF(double Phi,double ABar,double M); //Used for Mach Number computation
 
  static double GetFPrime(double Phi,double ABar,double M); //Used for Mach Number computation

  void OutputAllMachNumbers(const char* filename,int num); //outputs specified"num" value of arbitrary Mach numbers to visualize the non Mach Number sol. profile

  double ComputeMachNumber(); //Computes mach number using area-mach number relationship via Newton's method
  void ComputeExactSol(array<double,3> &sol); //Computes the exact solution for density,pressure, velocity & Mach Number
  void PrintResults(vector<double> &sol); //Used for printing out results


  ~SuperSonicNozzle();

};

#endif
