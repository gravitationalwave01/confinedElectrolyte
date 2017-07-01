//This code solves the poisson-boltzmann equation for a system of two plates seperated by an ionic solution
//Physical parameters are hard-coded, so changing the system requires recompilation (for now)
//Matrix inversion is done with LAPACK, when compiling you must link to LAPACK and BLAS
//compile with: g++ main.cpp -llapack -lblas -o testprog

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#define USE_MATH_DEFINES

extern "C" void dgtsv_(int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO); //lapack tridiagonal matrix solver



//this function solves the equation Ax=B, where A is the laplacian matrix, x is the electric potential, B is the charge density
//this equation must be solved iteratively until convergence
//the meanings of the variables are specified in main
//TODO: group the system variables into one struct (so that I don't pass them individually)
double iterate(int dim, double *psi, double *d, double *du, double *dl, double *b, double *rho1, double *rho2, double *kappa2, double boundary, double c1b, double c2b, double z1, double z2, double delta, int *b_cols, int *b_rows)
{
	int info = 0; // error code for dgtsv_
	dgtsv_(&dim, b_cols, dl, d, du, b, b_rows, &info); //call to lapack to solve Ax=b. Solution is stored in b vector.
	if (info != 0)
	{
		std::cerr << "PROBLEM SOLVING Ax=B, ERROR CODE " << info << std::endl;
		exit(1);
	}
		
	//update densities & compute whether we have converged
	double error = 0;
	for (int i=0; i<dim; i++)
	{
		error += pow(psi[i] - b[i],2); //very crude measure of convergence
		psi[i] = b[i];
		rho1[i] = c1b * z1 * exp(-z1*psi[i]);
		rho2[i] = c2b * z2 * exp(z2*psi[i]);
		kappa2[i] = 4 * M_PI * (z1*z1*rho1[i] + z2*z2*rho2[i]);
	}

	//dgtsv_ overwrites the d,dl, and du arrays -- we need to repopulate them.
	for (int i=0;i<dim-1;i++)
	{
		d[i] = -2.-kappa2[i]*delta*delta; 
		dl[i] = 1.;
		du[i] = 1.;
	}
	d[0] = -1. - kappa2[0]*delta*delta;
	d[dim-1] = -1. - kappa2[dim-1]*delta*delta;

	//update the RHS
	for (int i=0; i<dim;i++)
		b[i] = delta*delta*(4.*M_PI*(rho2[i] - rho1[i]) - kappa2[i]*psi[i]);
	b[0] = boundary;
	b[dim-1] =  -boundary; 

	return error;
}





int main()
{
	//Start by initializing some system parameters and constants:
    double z1 = 1.; //cation valancy 
    double z2 = 1.; //anion valancy
    double c1b = .1; //bulk concentration of cations (M)
    double c2b = .1; //bulk concentration of anions (M)
    double kb = 1.3806e-23; //boltzmann constant (m^2kgs^-2K^-1)
    double T = 250.; //temperature (K)
    double beta = 1/(kb*T); // 1/kT @ T=300K (1/J)
    double dielec = 80. * 8.85418e-12; // dielectric constant(F/m) for water
    double e_charge = 1.60218e-19; //electron charge (coulombs)
    double Av = 6.022e23; //avogadro's number
    double lb = beta*pow(e_charge,2)/(4*M_PI*dielec); //bjerrum length in bulk (m)
	double sigma = 0.02; //plate surface charge (e/nm^2) 
	double D = 10.e-9; //plate seperation (m)
    
	//First: make all physical quantities unitless
	// make the concentrations unitless 
    c1b *= 1.e3 * Av * lb*lb*lb;
    c2b *= 1.e3 * Av * lb*lb*lb;

	//make the charge density unitless:
	sigma *= pow(lb*1.e9,2);

	//make seperation dimensionless:
	D /= lb;


    int dim = 100; //spatial discretization: number of steps between the plates
	int b_cols = 1; // for lapack solver
	int b_rows = dim;

    //initialize variables for the Ax=B solver
	//A is tridiagonal matrix and will be stored in sparse format
	double *dl;  //subdiagonal of A
	double *d;   //diagonal of A
	double *du;  //superdiagonal of A
	double *b;   //RHS (total charge density)
	double *psi; // electric potential 
	double *rho1, *rho2; //density of cations/anions
	double *kappa2; //screening length squared
	dl = new double[dim-1]; 
	du = new double[dim-1]; 
	d = new double[dim]; 
	b = new double[dim];
	psi = new double[dim];
	rho1 = new double[dim];
	rho2 = new double[dim];
	kappa2 = new double[dim];

	double delta = D / (dim - 1); // discretization step size between the plates
	double boundary = 4. * M_PI * delta * sigma; // term for the discontinuity at the surface 

	//now we initialize the laplacian and the b vector
	//our first guess for the potential psi is psi increases linearly from one plate to the other 
	for (int i=0;i<dim;i++)
	{
		psi[i] = -sigma*D/2. + i*delta*sigma;
		rho1[i] = c1b * z1 * exp(-z1*psi[i]);
		rho2[i] = c2b * z2 * exp(z2*psi[i]);
		kappa2[i] = 4 * M_PI * (z1*z1*rho1[i] + z2*z2*rho2[i]);
		d[i] = -2. - kappa2[i]*delta*delta; //we subtract kappa2*psi from both sides of the equation to avoid singular matrix
		b[i] = delta*delta*(4.*M_PI*(rho2[i] - rho1[i]) - kappa2[i]*psi[i]); 
	}
	
	for (int i=0;i<dim-1;i++)
	{
		//Now we build the laplacian matrix
		du[i] = 1.;
		dl[i] = 1.; 
	}


	//now we incorporate the boundary conditions:	
	d[0] = -1.-kappa2[0]*delta*delta;  //again, subtracting kappa2*psi from both sides to avoid singularity
	d[dim-1] = -1.-kappa2[dim-1]*delta*delta; 
	
	b[0] = boundary;
	b[dim-1] = -boundary; 

	//below is the main loop of the code
	//basically, just solve Ax=B repeatedly until x stops changing (this change is measured by variable error)
	double error = 1.;
	int numiter = 0;
	while (error > 1.e-12 || numiter < 8) //for sanity I insist on >8 iterations. Technically unneccessary
	{
		//for (int i=0; i<dim;i++)//debugging
		//	std::cout << -D/2+i*delta << " " << rho1[i] << " " << rho2[i] << " " << psi[i] << " " << b[i] << std::endl;

		if (numiter > 1) std::cout <<"error is " << error <<std::endl;

		error = iterate(dim, psi, d, du, dl, b, rho1, rho2, kappa2, boundary, c1b, c2b, z1, z2, delta, &b_cols, &b_rows);
		if (error >10) {std::cerr << "large error  -- system diverges " << std::endl; break;} //TODO: catch and resolve divergences better
		numiter++;
	}

	//output results
	//TODO: format the output better.
	std::ofstream output;
	output.open("output.dat");
	for (int i=0;i<dim;i++)
		output << -D/2+i*delta << " " << rho1[i] << " " << rho2[i] << " " << psi[i] << std::endl;

	delete [] dl;
	delete [] du;
	delete [] d;
	delete [] b;
	delete [] psi;

    return(0);
}

