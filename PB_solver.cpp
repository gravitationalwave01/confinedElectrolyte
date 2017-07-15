//compilation:
// g++ -L/usr/bin/ -llapack -lopenblas -c PB_solver.cpp -o PB_solver.o

//This code solves the poisson-boltzmann equation for a system of two plates seperated by an ionic solution
//At the moment, the plates must have opposite charge. This is easily changed by 
//  setting b[dim-1] = +boundary rather than -boundary
//Physical parameters are passed in via structures 
//Matrix inversion is done with LAPACK, when compiling you must link to LAPACK and BLAS

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include "PB_solver.h"

#define USE_MATH_DEFINES

extern "C" void dgtsv_(int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO); //lapack tridiagonal matrix solver



//this function solves the equation Ax=B, where A is the laplacian matrix, x is the electric potential, B is the charge density
//this equation must be solved iteratively until convergence
//the meanings of the variables are specified in main
double iterate(params *sys, axb *solver, profile *p)
{
	int info = 0; // error code for dgtsv_
	dgtsv_(&solver->dim, &solver->b_cols, solver->dl, solver->d, solver->du, solver->b, &solver->b_rows, &info); //call to lapack to solve Ax=b. Solution is stored in b vector.
	if (info != 0)
	{
		std::cerr << "CRITICAL PROBLEM SOLVING Ax=B, ERROR CODE " << info << std::endl;
		std::cerr << " .... Please see the LAPACK DGTSV solver documentation to understand this error" << std::endl;
		exit(1);
	}
		
	//update densities & compute whether we have converged
	double error = 0;
	for (int i=0; i<solver->dim; i++)
	{
		error += pow(solver->x[i] - solver->b[i],2); //very crude measure of convergence
		solver->x[i] = solver->b[i];
		p->rho1[i] = sys->c1b * sys->z1 * exp(-sys->z1*solver->x[i]);
		p->rho2[i] = sys->c2b * sys->z2 * exp(sys->z2*solver->x[i]);
		p->kappa2[i] = 4 * M_PI * (sys->z1*sys->z1*p->rho1[i] + sys->z2*sys->z2*p->rho2[i]);
	}
	error /= solver->dim;
	
	//dgtsv_ overwrites the d,dl, and du arrays -- we need to repopulate them.
	for (int i=0;i<solver->dim-1;i++)
	{
		solver->d[i] = -2.-p->kappa2[i]*sys->delta*sys->delta; 
		solver->dl[i] = 1.;
		solver->du[i] = 1.;
	}
	solver->d[0] = -1. - p->kappa2[0]*sys->delta*sys->delta;
	solver->d[solver->dim-1] = -1. - p->kappa2[solver->dim-1]*sys->delta*sys->delta;

	//update the RHS
	for (int i=0; i<solver->dim;i++)
		solver->b[i] = sys->delta*sys->delta*(4.*M_PI*(p->rho2[i] - p->rho1[i]) - p->kappa2[i]*solver->x[i]);
	solver->b[0] = sys->boundary;
	solver->b[solver->dim-1] =  -sys->boundary; 

	return error;
}


int PB_solver(params* sys, profile* p, int num_points)
{
	std::cout << "entered PB solver" << std::endl;
	axb solver;
    solver.dim = num_points; //spatial discretization: number of steps between the plates
	solver.b_cols=1; // for lapack solver
	solver.b_rows=solver.dim; //for lapack solver
	solver.dl=new double[solver.dim-1];  //subdiagonal of A
	solver.d=new double[solver.dim];   //diagonal of A
	solver.du=new double[solver.dim-1];  //superdiagonal of A
	solver.b=new double[solver.dim];   //RHS (total charge density)
	solver.x=p->psi; //Ax=b. The x is the electric potential
							//note: we are assigning the pointer rather than copying

	sys->delta = sys->D / (solver.dim - 1); // discretization step size between the plates
	sys->boundary = 4. * M_PI * sys->delta * sys->sigma; // term for the discontinuity at the surface 


	


	bool diverge = false;
	//for (int iter=0;iter<num_iter;iter++)
	{
		diverge = false;

		//now we initialize the laplacian and the b vector
		//our first guess for the potential psi is psi increases linearly from one plate to the other 
		for (int i=0;i<solver.dim;i++)
		{
			//solver.x points to sys.psi (it's the electric potential)
			solver.x[i] = -sys->sigma*sys->D/2. + i*sys->delta*sys->sigma;
			p->rho1[i] = sys->c1b * sys->z1 * exp(-sys->z1*solver.x[i]);
			p->rho2[i] = sys->c2b * sys->z2 * exp(sys->z2*solver.x[i]);
			p->kappa2[i] = 4 * M_PI * (sys->z1*sys->z1*p->rho1[i] + sys->z2*sys->z2*p->rho2[i]);
			solver.d[i] = -2. - p->kappa2[i]*sys->delta*sys->delta; //we subtract kappa2*psi from both sides of the equation to avoid singular matrix
			solver.b[i] = sys->delta*sys->delta*(4.*M_PI*(p->rho2[i] - p->rho1[i]) - p->kappa2[i]*solver.x[i]); 
		}
		
		for (int i=0;i<solver.dim-1;i++)
		{
			//Now we build the laplacian matrix
			solver.du[i] = 1.;
			solver.dl[i] = 1.; 
		}


		//now we incorporate the boundary conditions:	
		solver.d[0] = -1.-p->kappa2[0]*sys->delta*sys->delta;  //again, subtracting kappa2*psi from both sides to avoid singularity
		solver.d[solver.dim-1] = -1.-p->kappa2[solver.dim-1]*sys->delta*sys->delta; 
		
		solver.b[0] = sys->boundary;
		solver.b[solver.dim-1] = -sys->boundary; 

		//below is the main loop of the code
		//basically, just solve Ax=B repeatedly until x stops changing (this change is measured by variable error)
		double error = 1.;
		int numiter = 0;
		while (error > 1.e-12 || numiter < 8) //for sanity I insist on >8 iterations. Technically unneccessary
		{
			//for (int i=0; i<solver.dim;i++)//debugging
			//	std::cout << -D/2+i*delta << " " << p->rho1[i] << " " << p->rho2[i] << " " << psi[i] << " " << b[i] << std::endl;

		//	if (numiter >= 1) std::cout <<"error is " << error <<std::endl;

			error = iterate(sys, &solver, p);
			if (numiter >= 50) {diverge = true; std::cerr << "large error  -- system diverges " << std::endl; break;} //TODO: catch and resolve divergences better
			numiter++;
		}
		if (diverge)
		{
			std::cerr << "CRITICAL ERROR: PB equation solver diverged. Computation aborted." << std::endl;
			exit(1);
		}

/*
		if (!diverge)
		{
			//output results
			//TODO: format the output better.
			std::stringstream ss;
			ss << "./output/density_potential_" << iter << ".dat";
			std::ofstream output(ss.str().c_str());
			for (int i=0;i<solver.dim;i++)
				output << -sys->D*sys->lb/2+i*sys->delta*sys->lb << " " << p->rho1[i] << " " << p->rho2[i] << " " << solver.x[i]*kb*T/e_charge << " " << std::endl;
			cap << sys->sigma*e_charge/pow(sys->lb,2) << " " << sys->sigma*e_charge/pow(sys->lb,2)/((solver.x[solver.dim-1]-solver.x[1])*kb*T/e_charge) << std::endl;
		}
*/

	}


	delete [] solver.dl;
	delete [] solver.du;
	delete [] solver.d;
	delete [] solver.b;

    return(0);
}

