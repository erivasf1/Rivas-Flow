//User-defined functions
#include "ExactNozzle.h" 

// TOOLS DEFINITIONS
//-----------------------------------------------------------------------
void Tools::print(const char format[],...) {

  va_list Argp;
  va_start(Argp, format);
  vprintf(format, Argp);
  va_end(Argp);
  return;

}

//-----------------------------------------------------------------------
void Tools::printWithPrecision(int precision,const char* format,...) {
  //Code generated from ChatGPT
  std::ostringstream output;
    output << std::fixed << std::setprecision(precision);

    // Initialize the variadic argument list
    va_list args;
    va_start(args, format);

    const char* ptr = format;
    while (*ptr != '\0') {
        if (*ptr == '%' && *(ptr + 1) == 'f') { // Handle floating-point formatting
            double value = va_arg(args, double);
            output << value;
            ptr++; // Skip 'f'
        } else {
            output << *ptr;
        }
        ptr++;
    }

    va_end(args);

    // Print the final formatted string
    std::cout << output.str() << std::endl;

}




//-----------------------------------------------------------------------
vector<double> Tools::RetrievePoints(double xmin,double xmax,int pt_num){

  vector<double> pts;
  //pts.push_back(xmin); -- used in HW1

  double dx = abs(xmax-xmin) / (pt_num-1.0); 
  for (int k=0;k<pt_num;k++)
    pts.push_back(xmin + k * dx); 
    //pts.push_back(xmin+(dx*(k+1))); -- used in HW1
   
  //debug:
  //print("Points:%f,%f,%f,%f,%f\n",pts[0],pts[1],pts[2],pts[3],pts[4]);

  return pts;
}

double Tools::AreaVal(double x){

  //A(x) = 0.2 + 0.4[1+sin(pi(x-0.5)]
  double PI = M_PI;
  double res = 0.2 + 0.4*(1+sin(PI*(x-0.5)));
  return res;
}
//-----------------------------------------------------------------------
double Tools::ComputeDotProduct(double &x1,double &x2,double &y1,double &y2){

  double res = (x1*x2) + (y1*y2);

  return res;

}
//-----------------------------------------------------------------------

// SUPERSONICNOZZLE DEFINITIONS
SuperSonicNozzle::SuperSonicNozzle(double &a,double &b,double &c,double &d,bool &e) //constructor
    : area(a), area_star(b), stag_pressure(c), stag_temp(d), cond(e){}


SuperSonicNozzle::~SuperSonicNozzle(){} //destructor


//---------------------------------------------------------
double SuperSonicNozzle::GetPhi(double M) {

  double phi = 2.0/(gamma+1.0);
  phi *= 1.0 + ((gamma-1.0)/2.0) * pow(M,2.0);
  //phi *= (1.0 + ((gamma-1.0)/2.0)) * pow(M,2.0);
  return phi;
}

//---------------------------------------------------------
double SuperSonicNozzle::GetF(double Phi,double ABar,double M) {

  double f = pow(Phi,(gamma+1.0)/(gamma-1.0));
  f -= pow(ABar*M,2.0);
  return f;

}

//---------------------------------------------------------
double SuperSonicNozzle::GetFPrime(double Phi,double ABar,double M) {

  double fprime = 2.0 * M;
  fprime *= pow(Phi,2.0/(gamma-1.0)) - pow(ABar,2.0);
  return fprime;

}

//---------------------------------------------------------
void SuperSonicNozzle::OutputAllMachNumbers(const char* filename,int num){

  ofstream myfile;
  myfile.open(filename);

  if (!myfile.is_open()){ //checking if file opened successfully
    cerr<<"Error: Could Not Open File "<<filename<<endl;
    return;
  }

  vector<double> MachNumbers = Tools::RetrievePoints(0.0,4.0,num); //arbitrary num. of Mach Numbers


  double phi,F;
  double area_bar = area/area_star;

  myfile<<"Mach Number"<<"  "<<"F(Residual)"<<endl;
  for (int n=0;n<num;n++){
    phi = GetPhi(MachNumbers[n]);  
    F = GetF(phi,area_bar,MachNumbers[n]);
  
    myfile<<MachNumbers[n]<<"  "<<F<<endl;
  
  }
  
  myfile.close(); //closing file writing to it
  //myfile.flush();

  return;

}

//---------------------------------------------------------
double SuperSonicNozzle::ComputeMachNumber(){
//F[m] = [2/gamma+1(1+(gamma-1)/2)M^2]^[(gamma+1)/(gamma-1)]
// F'[M] = 2M[Phi^(2/gamma-1)-\Bar{A}^2] & Phi = [2/gamma+1(1+(gamma-1)/2)M^2]

  double M0,M1; //old (initial guess) and converged sol., respectively
  double ABar = area/area_star; //ratio of present area to throat area
  //print("Abar:%f\n",Abar);
  double F,Fprime; //fcn. and fcn. derivative
  double phi;
  double resid; // residual

  //debug:
  //print("area: %f and tol.: %f\n",area,tol);
  
  M0 = (cond == true) ? 0.1:4.0; // initial guess of 0.1 for subsonic & 4.0 for supersonic region
  //print("M0 initial guess:%f\n",M0); //for debugging
  M1 = M0; //initial guess & resid.
  resid = 1.0e5;  

  //Computing initial residual
  double phi_init = GetPhi(M0);
  double F_init = GetF(phi_init,ABar,M0);
  
  //Employing Hybrid Bracketing (Bisection) Newton's Method to solve nonlinear equation

  double Ma,Mb; //end brackets for M 
  double Fa,Fb; //function vals(Fs) for end brackets
  double phi_a,phi_b; //end bracket phi values

  Ma = (cond == true) ? 1.0e-4:1.0; //1st end bracket (true for subsonic case)
  Mb = (cond == true) ? 1.0:5.0; //2nd end bracket (true for subsonic case)
  for (int i=0;i<maxiter;i++) {
    if (resid <= tol) break;

    // End interval function values
    phi_a = GetPhi(Ma);
    phi_b = GetPhi(Mb);
    Fa = GetF(phi_a,ABar,Ma);
    Fb = GetF(phi_b,ABar,Mb);
    if ((Fa > 0.0 && Fb > 0.0) || (Fa < 0.0 && Fb < 0.0)){
      print("End intervals are same functional sign. No solution exists!\n"); 
      return 0.0;
    }
    

    // Stepping via Newton's Method
    phi = GetPhi(M1);
    F = GetF(phi,ABar,M1);
    Fprime = GetFPrime(phi,ABar,M1);
    M1 -= F/Fprime;

    // Checking if M1 is outside interval
    if (M1>Mb || M1<Ma)
      M1 = (Mb+Ma) / 2.0; //M1 computed via bi-section method
    

    //checking residual
    phi = GetPhi(M1);
    F = GetF(phi,ABar,M1);
    resid = abs(F/F_init);

    // updating end brackets
    if ((Fa > 0.0 && F > 0.0) || (Fa < 0.0 && F < 0.0)) //updating Ma case
      Ma = M1;

    if ((Fb > 0.0 && F > 0.0) || (Fb < 0.0 && F < 0.0)) //updating Mb case
      Mb = M1;
    //resid = abs(F);
    //F_init = F;

    //debug mode
    //print("M1: %f & residual:%f\n",M1,resid);
  }
  return M1;

}

//---------------------------------------------------------
void SuperSonicNozzle::ComputeExactSol(array<double,3> &sol){

  //vector<double> sol;
  double T,P,Rho,V;
  double M,Psi;
  double R = Ru/MolMass; //specific gas constant
  //stag_pressure *= 1000.0;

  //print("Point: %f and Area: %f\n",pts[n],area);
  M = ComputeMachNumber();
  //print("Mach Number: %f\n",M);
  Psi = 1.0 + ((gamma-1.0)*0.5 * pow(M,2.0));

  //Temperature
  // Psi = 1 + (gamma-1/2)M^2
   T = stag_temp/Psi;
   //print("Temperature ratio: %f\n",T/stag_temp);
  
  //Pressure -- convert from kPa
  // P = P0/[Psi^(gamma/gamma-1)] 
   P = stag_pressure / pow(Psi,gamma/(gamma-1.0));
   //print("Pressure ratio: %f\n",P/stag_pressure);

  //Density
  // Rho = P / (RT)
   Rho = P / (R*T);
  // double Rho_t = stag_pressure / (R*stag_temp);
   //print("Density ratio: %f\n",Rho/Rho_t);
  
  //Velocity
  // V = sqrt(gamma*R*T)*M
   V = abs(M * sqrt(gamma*R*T));

  //Appending flow quantities sol. vector of point
  /*sol.push_back(Rho);
  sol.push_back(V);
  sol.push_back(P);
  sol.push_back(M);
  */
  sol[0] = Rho;
  sol[1] = V;
  sol[2] = P;
    //debug:
    //print("Point: %f\n",pts[n]);
    //print("Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",Rho,V,P,M);

  
  
  //PrintResults(sol);
}

void SuperSonicNozzle::PrintResults(vector<double> &sol) {

  /*print("---------\n");
  print("RESULTS\n");
  print("---------\n");*/


  //print("Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",sol[0],sol[1],sol[2],sol[3]);
  printWithPrecision(14,"Rho[kg/m^3] = %f,Velocity[m/s] = %f,Pressure[kPa] = %f,Mach Number = %f\n",sol[0],sol[1],sol[2],sol[3]);
  print("-------\n");

}



