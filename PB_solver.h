// PB_SOLVER.h

#ifndef _PB_solver // must be unique name in the project
#define _PB_solver

//Structure containing all of the physical parameters of the system
struct params{
	//all of the physical parameters of the system:
    double z1; //cation valancy 
    double z2; //anion valancy
    double c1b; //bulk concentration of cations (M)
    double c2b; //bulk concentration of anions (M)
    double lb; //bjerrum length in bulk (m)
	double sigma; //plate surface charge (C/m^2) 
	double D; //plate seperation (m)
	double boundary; // boundary condition for poisson equation
	double delta; // spatial incement
	double diam1; //diameter of cation
	double diam2; //diameter of anion
};




//structure for grouping all the electrostatic data
struct profile{
	double *rho1; //cation density
	double *rho2; //anion density
	double *self_energy; //
	double *greens; //spatial greens function
	double *psi;
	double *kappa2;
};



//structure for storing the matrices/vectors for the Ax=B problem
struct axb{
	//A is tridiagonal matrix and will be stored in sparse format
    int dim; //spatial discretization: number of steps between the plates
	int b_cols; // for lapack solver
	int b_rows;// for lapack solver
	double *dl;  //subdiagonal of A
	double *d;   //diagonal of A
	double *du;  //superdiagonal of A
	double *b;   //RHS (total charge density)
	double *x; // unknown (Ax=B)
};


int PB_solver(params* sys, profile* p, int num_points); // prototype declaration of the function in a.cpp

#endif 


