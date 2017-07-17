#include<iostream>
#include<fstream>
#include<math.h>
#include<sstream>
#include "PB_solver.h"

#define _USE_MATH_DEFINES

//The goal here is to implement the WKB method for an attractive plate.
//We first solve for the Greens function assuming that the ionic strength is constant 
//(and equal to the bulk ionic strength) 
//Then we plug the varying ionic strength into the resulting epxression for the greens function


// compile with:
// g++ -std=c++11 -I/usr/lib -llapack -lopenblas full_WKB.cpp PB_solver.cpp -o ./full_WKB
extern "C" void dgtsv_(int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO);


//The plate seperation is D
//We never want to evaluate at the plate wall
//getpos converts the index position i into a position
//between the plates (evenly spaced)
double getpos(int i, int num_points, double D) //units of lb or ion diameter
{
	return (-D/2.0 + (double)(i) / ((double)num_points-1.) * D);
}
	
//returns true of the given index is within a half-diameter of the wall
bool HS(int index, int num_points, double D, double diam)
{
	return (D/2. - fabs(getpos(index,num_points,D)) < diam/2.);
}



double Jkx(double x, double kappa,double D, double f)
{
	double ans = 0; //This represents the sum (eq 7 in Rui's 2013 paper). We include enough terms in the sum so that it converges
	double oldJkx = 1;
	int m=1;
	while (m < 10 || fabs(ans-oldJkx) > 1.0e-10)
	{
		//std::cout << m*D + 2*x << std::endl;
		oldJkx = ans;
		if (m % 2 == 0)
			ans += pow(f,m)*2.0/(m*D) * exp(-kappa * m * D);
		else
			ans += pow(f,m)*(exp(-kappa*(m*D+2.0*x))/(m*D+2*x) + exp(-kappa*(m*D-2.0*x))/(m*D-2*x));
		m+=1;
	}
	//std::cout << "m is " << m << "  J is " << ans << std::endl;
	return ans;
}



//This function computes the charging integral F_fl:
//everything in this function is unitless
double Ffl(double x, double kappa, double D, double f)
{
	double Ffl = 0;
	double numsteps = 400.0;
	double dalpha = 1.0/(numsteps+1.0);
	int counter =0;
	for (double alpha = 0; alpha < 1; alpha += dalpha)
	{
		counter++;
		if (counter%1000==0) std::cout << " ........... in FFL, counter is " << counter << std::endl; 
		Ffl += Jkx(x,alpha*kappa,D,f)*alpha - alpha*kappa;
	}
    std::cout << "FFL term is " << Ffl << std::endl;
	return kappa*kappa * Ffl * dalpha / (4*M_PI);
}	


int main()
{
	//Start by initializing some system parameters and constants:
    double kb = 1.3806e-23; //boltzmann constant (m^2kgs^-2K^-1)
    double T = 300.; //temperature (K)
    double beta = 1/(kb*T); // 1/kT @ T=300K (1/J)
    double dielec = 80. * 8.85418e-12; // dielectric constant(F/m) for water
    double e_charge = 1.60218e-19; //electron charge (coulombs)
    double Av = 6.022e23; //avogadro's number
	
	params sys;
	bool point_charge = false; //set to true if you want point charges
							   //length unit is then the bjerrum length
    sys.z1 = 1.; //cation valancy 
    sys.z2 = 1.; //anion valancy
	sys.diam1 = 5.e-10; //diameter of ions [for repulsion from wall]
	sys.diam2 = 5.e-10;
    sys.c1b = 1.e-3; //bulk concentration of cations (M)
    sys.c2b = 1.e-3; //bulk concentration of anions (M)
    sys.lb = beta*pow(e_charge,2)/(4*M_PI*dielec); //bjerrum length in bulk (m)
	sys.sigma = 0.; //plate surface charge (C/m^2) 
	sys.D = 3.e-9; //plate seperation (m)

	//First: make all physical quantities unitless
	double lu = std::min(sys.diam1,sys.diam2);
	if (point_charge) lu = sys.lb;

	// make the concentrations unitless 
    sys.c1b *= 1.e3 * Av * lu*lu*lu;
    sys.c2b *= 1.e3 * Av * lu*lu*lu;

	//make the charge density unitless:
	sys.sigma *= pow(lu,2) / e_charge;

	//make seperation and ion radii dimensionless:
	sys.D /= lu;
	sys.diam1 /= lu;
	sys.diam2 /= lu;

//	std::cout << "lb is " << sys.lb << std::endl;
//	std::cout << "D is " << sys.D << std::endl;
	std::cout << "diameter1 is " << sys.diam1 << std::endl;
	std::cout << "diameter2 is " << sys.diam2 << std::endl;


	double epsilS = 80.;
	double epsilP = 2.5;
	double kappab2 = 4*M_PI*(sys.z1*sys.z1 * sys.c1b + sys.z2*sys.z2 * sys.c2b); //bulk inverse screening length 


	//For each plate seperation D we will:
	//1) guess a greens function by assuming kappa = kappab everywhere
	//   and incorporating image charges & correlations
	//2) Compute resulting self energy, ion concentrations, potential
	//3) compute a new guess for the greens function
	//4) repeat 2 and 3 until convergence
	//5) compute the resulting free energy per area G
	//6) compute the resulting pressure dG/dD

	FILE* fenergy;
	fenergy = fopen("./output/test/total_energy.dat","w");
	fprintf(fenergy,"seperation    energy \n");
	double sigma_init = sys.sigma;
	double c1b_init = sys.c1b;
	int cursig = 0;
	double f = (epsilS - epsilP)/(epsilS + epsilP);
	//for (sys.D = 1.; sys.D < 30; sys.D += 1.) //D is unitless (multiples of lb)
	//for (sys.sigma = -sigma_init*10; sys.sigma < sigma_init * 10; sys.sigma += sigma_init*0.1)
    for (sys.D = 1.5; sys.D < 20; sys.D += 0.1) //in units of ion diameter (or bjerrum length)
	{
        std::cout << "D is " << sys.D << std::endl;
	
		//open the output file: (TODO: generate the required directories)
		std::stringstream ss;
		ss << "./output/test/ion_profile_D" << cursig << ".dat" ;
		cursig++;
		FILE* output;
		output = fopen(ss.str().c_str(),"w");
		fprintf(output,"pos      rho1       rho2       psi       u      kappa\n");
	

		int num_points=(int)(sys.D*2.*5.); //spatial discretization
		std::cout << "num_points is " << num_points << std::endl;
		std::cout << "Exclusion zone ends at " << -sys.D/2.0+sys.diam1/2.0 << std::endl;

		profile p;
		p.rho1 = new double[num_points];
		p.rho2 = new double[num_points];
		p.greens = new double[num_points];
		p.psi = new double[num_points];
		p.self_energy1 = new double[num_points];
		p.self_energy2 = new double[num_points];
		p.kappa2 = new double[num_points];
		
		//first guess:
		for (int i=0;i<num_points;i++)
		{

			if (HS(i,num_points,sys.D,sys.diam1))
				p.rho1[i]=0;
			else 
				p.rho1[i] = sys.c1b;
			
			if (HS(i,num_points,sys.D,sys.diam2))
				p.rho2[i]=0;
			else 
				p.rho2[i] = sys.c2b;

			p.kappa2[i] = 4*M_PI*(sys.z1*sys.z1 * p.rho1[i] + sys.z2*sys.z2 * p.rho2[i]);

			p.self_energy1[i] = 0.5*sys.z1*sys.z1*(Jkx(getpos(i,num_points,sys.D),sqrt(p.kappa2[i]),sys.D,f) - sqrt(p.kappa2[i]));
			p.self_energy2[i] = 0.5*sys.z2*sys.z2*(Jkx(getpos(i,num_points,sys.D),sqrt(p.kappa2[i]),sys.D,f) - sqrt(p.kappa2[i]));
			
			p.psi[i] = 0;
		}
		
		//Now we iterate:
		double error = 1.;
		int counter =0;
		while (error > 1.e-12 || counter < 5)
		{
		    counter++;
/*
		std::stringstream fuck;
		fuck << "./output/debug" << counter  << ".out";
		FILE* debug;
		debug = fopen(fuck.str().c_str(),"w");
		fprintf(debug,"pos  rho1   rho2   kappa2  u     psi\n");

			fprintf(debug," ..... new iteration .......\n");
*/
			
			//Solve modified PB equation
			PB_solver(&sys,&p,num_points);

			// We now have psi, rho1, rho2, kappa2 all populated with solutions to the PB equation
			// Now we guess new rho1 and rho2 by taking into account self energy and the psi from PB equation
			error = 0;
			for (int i=0;i<num_points;i++)
			{
				//compute new charge densities
				if (HS(i,num_points,sys.D,sys.diam1))
					p.rho1[i]=0;
				else 
					p.rho1[i] = sys.c1b * exp(-0.5*sys.z1*sys.z1*sqrt(kappab2) - sys.z1*p.psi[i]- p.self_energy1[i]);;
			
				if (HS(i,num_points,sys.D,sys.diam2))
					p.rho2[i]=0;
				else 
					p.rho2[i] = sys.c2b * exp(-0.5*sys.z2*sys.z2*sqrt(kappab2) + sys.z2*p.psi[i] - p.self_energy2[i]);

				p.kappa2[i] = 4*M_PI*(sys.z1*sys.z1 * p.rho1[i] + sys.z2*sys.z2 * p.rho2[i]);

				//compute the new self energy (error is computed as the change in self energy from last iteration)
				double tmp;
				if (HS(i,num_points,sys.D,sys.diam1))
					tmp = 0;
				else 
					tmp = 0.5*sys.z1*(Jkx(getpos(i,num_points,sys.D),sqrt(p.kappa2[i]),sys.D,f) - sqrt(p.kappa2[i])); 
	
				error += pow(p.self_energy1[i] - tmp, 2);
				p.self_energy1[i] = tmp;
				if (HS(i,num_points,sys.D,sys.diam1))
					p.self_energy2[i] = 0;
				else 
					p.self_energy2[i] = 0.5*sys.z2*(Jkx(getpos(i,num_points,sys.D),sqrt(p.kappa2[i]),sys.D,f) - sqrt(p.kappa2[i])); 


			}
			error /= num_points;
			std::cout << "WKB ERROR IS " << error << std::endl;
		}	
		
		//Now we have self consistently determined self energy, kappa2, rho1, rho2
		//Next compute the total grand energy (per unit area) using a charging integral (see Rui's 2013 WKB paper).
		double *dpsi = new double[num_points]; //derivative of psi 
		dpsi[0] = (p.psi[1]-p.psi[0])/sys.D/(num_points+1);
		for (int i=1;i<num_points-1;i++)
			dpsi[i] = (p.psi[i+1]-p.psi[i-1])/(2.*sys.D/(num_points+1));
		dpsi[num_points-1] = (p.psi[num_points-1]-p.psi[num_points-2])/sys.D/(num_points+1);

//		std::cout << "..... computed derivative of psi" << std::endl;

		double G = 0;//total grand energy beta*G
		//potential is 1/8pi * dpsi/dz 
		for (int i=0;i<num_points;i++)
		{
 //           std::cout << "-------------------------" << std::endl;
			G += 1./8./M_PI *pow(dpsi[i],2);
//			std::cout << "added psi " << G << std::endl;

			if (!HS(i,num_points,sys.D,sys.diam1))
				G += p.rho1[i]*(log(p.rho1[i]/sys.c1b)+0.5*sys.z1*sys.z1*sqrt(kappab2) - 1);
			if (!HS(i,num_points,sys.D,sys.diam2))
				G += p.rho2[i]*(log(p.rho2[i]/sys.c2b)+0.5*sys.z2*sys.z2*sqrt(kappab2) - 1);

//			std::cout << "added rho " << G << std::endl;
//			std::cout << ".......took log of " << p.rho2[i]/sys.c2b << std::endl;
			if (!HS(i,num_points,sys.D,sys.diam2)) // WARNING: THIS ONLY WORKS WHEN DIAMETERS ARE THE SAME!!!
				G += Ffl(getpos(i,num_points,sys.D),sqrt(p.kappa2[i]),sys.D,f);
//			std::cout << "added ffl " << G << std::endl;
		}
//		std::cout << "..... computed energy G" << std::endl;
		
		fprintf(fenergy,"%f %f %f\n",sys.D,G*sys.D/num_points,-0.5*sys.z2*sys.z2*sqrt(kappab2)-0.5*sys.z1*sys.z1*sqrt(kappab2));
        fflush(fenergy);

		for (int i=0;i<num_points;i++)
			fprintf(output,"%f %f %f %f %f %f\n",getpos(i,num_points,sys.D),p.rho1[i]/sys.c1b,p.rho2[i]/sys.c2b,p.psi[i], p.self_energy1[i]+0.5*sqrt(p.kappa2[i])*sys.lb, p.kappa2[i]);		

		fclose(output);
	delete [] p.rho1;
	delete [] p.rho2;
	delete [] p.greens;
	delete [] p.psi;
	delete [] p.self_energy1;
	delete [] p.self_energy2;
	delete [] p.kappa2;
	}
//	std::cout << "terminated succesfully" << std::endl;
	exit(0);
}

	
		
