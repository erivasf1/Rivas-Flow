//1D FVM Data Allocation

1) Mesh:
  1.1) in user-provided mesh specs the following are given:
    1.1a) x-range (xmin and xmax)
    1.1b) cellnum (number of cells)
  1.2) the following is outputted from GenerateMesh fcn.:
    1.2a) xcoords: list of xcoord locations for every CELL FACE including xmin and xmax (pretty much a linspace of xmin,xmax,face# which = cellnum+1)

2) Data indexing
  2.1) Before adding ghost cells (SetBoundaryConditions fcn.)
    2.1a) cell_list[0,....,cellnum-1] -- only consisting of interior nodes
    2.1b) identifying corresponding cell faces: cell_i = cell_list[i]; right face = xcoords[i+1] & left face = xcoords[i]
  2.2) After adding 4 ghost cells (2 per side of nozzle)
    2.2a) cell_list[0,....,cellnum+3] -- consisting of interior + ghost cells 
    2.2b) identifying corresponding cell faces: cell_i = cell_list[i]; right face = xcoords[i-1] & left face = xcoords[i-2]

CURRENTLY:
Copies labeled "TEST" are copies of its corresponding file with the main nuance being the absence of pointers. Try getting this to run first (hopefully converge) and then use pointers

3) Under-relaxation Algorithm (enforced after the 1st timestep & is equation specific for now)
  3.1) Forward Euler Advance to obtain new primtive variable solutions at next time step (this will be stored as an intermediate step)
  3.2) Compute the Residual of the intermediate step
  3.3) Check if residual norm of intermediate step is in "bad state" (i.e. if greater than previous residual by a factor of C (typical values are 1.1-1.5))
  3.4a) if residual norm is in "bad state", then recompute Euler Advance with a smaller CFL (i.e. reduce time-step)
  3.4b) if residual norm is in "good state", then assign new time step values to the intermediate values
  3.5) Repeat

THINGS LEARNED:
1) Residuals really are just checks on how well the cell(or domain) satisfies the governing equations (i.e. conservation laws). If the residual is high then solver will want to "greatly" change the values in the cell in hopes of it generating a lesser residual

TODO:
1) Figure out a way to best store Field vectors 

//2D FVM Solver Extension Idea:
1) Mesh class will be only class to be split into a 1D version and 2D version 
  1.1) 1D version - generates a uniform 1D cartesian grid, really is for the quasi-steady nozzle
  1.2) 2D version - reads in a structured Plot3D file and extracts all coords and cell number
2) Everything else will be re-written to extend into 2D.
  2.1) Ramification - for 1D scenarios, there will be a "dummy" variable that won't be solved for, however this is justified by a cleaner written source code
  2.2) There should still be a tag/variable that indicates if the flow is 1D or 2D.
3) There should be an input file to clean up code as well. Consider using something like IoData.cpp/.h

TODO!!!:
1) Modify all pertaining fcns. to handle Field Variables of vector,array size of 4
2) Create BASE classes for EulerOperator Class, to handle all requirements (1D quasi-steady nozzle,2D generic flow, 2D MMS flow)
  2.1) EulerBASE
  2.2) Euler1D
  2.3) Euler2D
  2.4) Euler2DMMS
3) Create BASE classes for TimeIntegrator Class, to handle all requirements (1D quasi-steady nozzle,2D generic flow, 2D MMS flow)
  3.1) ExplicitBASE
  3.2) EulerExplicit 
  3.3) RungeKutta2
  3.4) RungeKutta4
--BOUNDARY CONDITIONS
1) Manufactured Sols. Case
  1.1) Top - Outflow
  1.2) Bottom - Inflow
  1.3) Left - Inflow
  1.4) Right - Outflow
1) 30 deg. inlet:
  1.1) Top - Slip Wall
  1.2) Bottom - Slip Wall
  1.3) Left - Inflow
  1.4) Right - Outflow
1) Airfoil?
  1.1) Top - Slip Wall
  1.2) Bottom - Slip Wall
  1.3) Left - Inflow
  1.4) Right - Outflow
--IDEAS:
1) Set manufactured sol. as initial condition to start off
