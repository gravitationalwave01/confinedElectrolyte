#include<iostream>
#include<fstream>
#include<math.h>
#include<sstream>

#define _USE_MATH_DEFINES

//The goal here is to implement the WKB method for an attractive plate.
//We first solve for the Greens function assuming that the ionic strength is constant 
//(and equal to the bulk ionic strength) 
//Then we plug the varying ionic strength into the resulting epxression for the greens function


extern "C" void dgtsv_(int* N, int* NRHS, double* DL, double* D, double* DU, double* B, int* LDB, int* INFO);

double Jkx(double x, double kappa,double D, double f)
{
	double ans = 0; //This represents the sum (eq 7 in Rui's 2013 paper). We include enough terms in the sum so that it converges
	double oldJkx = 1;
	int m=1;
	while (m < 10 || fabs(ans-oldJkx) > 1.0e-4)
	{
		//std::cout << m*D + 2*x << std::endl;
		oldJkx = ans;
		if (m % 2 == 0)
			ans += pow(f,m)*2.0/(m*D) * exp(-kappa * m * D);
		else
			ans += pow(f,m)*(exp(-kappa*(m*D+2.0*x))/(m*D+2*x) + exp(-kappa*(m*D-2.0*x))/(m*D-2*x));
		m=m+1;
	}
	std::cout << "m is " << m << "  J is " << ans << std::endl;
	return ans;
}



//This function computes the charging integral F_fl:
//everything in this function is unitless
double Ffl(double x, double kappa, double D)
{
	double Ffl = 0;
	double numsteps = 5.0;
	double dalpha = 1.0/(numsteps+1.0);
	for (double alpha = 0; alpha < 1; alpha += dalpha)
	{
		Ffl += Jkx(x,alpha*kappa,D,-1)*alpha - alpha*kappa;
	}
	return kappa*kappa * Ffl * dalpha / (4*M_PI);
}	

//The plate seperation is D
//We never want to evaluate at the plate wall
//getpos converts the index position i into a position
//between the plates (evenly spaced)
double getpos(int i, int num_points, double D) //units of lb
{
	return (-D/2.0 + (double)(i+1) / (num_points+2) * D);
}
	
int main()
{
	//Start by initializing some constants:
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
	double lb = beta*pow(e_charge,2)/(4*M_PI*dielec); //bjerrum length (m)

	// make the concentrations unitless	
	c1b *= 1.e3 * Av * lb*lb*lb;
	c2b *= 1.e3 * Av * lb*lb*lb;
	std::cout << "the bulk concentration is " << c1b << std::endl;

	double* cation; //concentration of cations (unitless)
	double* anion; //concentration of anions (unitless)
	int num_points; //spatial discretization
	double* psi; //electrostatic potential
	double* greens; //spatial greens function (unitless)
	double* self_energy; //spatial self energy (unitless)
	double* kappa; //spatial inverse screening length (unitless)


	//let's test linking to lapack:
	double *d,*du,*dl,*b;
	int dim = 10;
	int b_row = dim;
	int b_col = 1;
	int info = 0;
	for (int i=0;i<dim-1;i++)
	{
		d[i] = 2.;
		du[i] = -1.;
		dl[i] = -1.;
		b[i] = 1;
	}
	d[dim-1]=2.;
	b[dim-1]=1.;

	dgtsv_(&dim, &b_row, dl, d, du, b, &b_col, &info);
	std::cout << "SUCCESS!" << std::endl;
	exit(0);	

	//For each plate seperation D we will:
	//1) guess a greens function by assuming kappa = kappab everywhere
	//   and incorporating image charges & correlations
	//2) Compute resulting self energy, ion concentrations, potential
	//3) compute a new guess for the greens function
	//4) repeat 2 and 3 until convergence
	//5) compute the resulting free energy per area G
	//6) compute the resulting pressure dG/dD

	for (double D = 2.; D < 40; D += 5.) //D is unitless (multiples of lb)
	{
		std::stringstream ss;
		ss << "ion_profile_D" << D;
		std::ofstream output;
		output.open(ss.str().c_str());
		num_points = 1000;
		cation = new double[num_points];
		anion = new double[num_points];
		greens = new double[num_points];
		psi = new double[num_points];
		self_energy = new double[num_points];
		kappa = new double[num_points];
		double epsilS = 80;
		double epsilP = 2.5;
		double f = (epsilS - epsilP)/(epsilS + epsilP);
		double kappab = sqrt(4*M_PI*(z1*z1 * c1b + z2*z2 * c2b)); //bulk inverse screening length 
		
		//first guess:
		for (int i=0;i<num_points;i++)
		{
			kappa[i] = kappab; 
			self_energy[i] = 0.5*(Jkx(getpos(i,num_points,D),kappab,D,f) - kappab); 
			double tmp = c1b * exp(-0.5*z1*z1*lb*kappab-self_energy[i]);
			cation[i] = tmp;
			anion[i] = c2b * exp(-0.5*z2*z2*lb*kappab-self_energy[i]);
			psi[i] = 0;
		}
		
		//Now we iterate:
		double error = 1.;
		while (error > 1.e-8)
		{
			std::cout << error << std::endl;
			error = 0;
			for (int i=0;i<num_points;i++)
			{
			
				kappa[i] = sqrt(4*M_PI*(z1*z1*cation[i]+z2*z2*anion[i]));
				self_energy[i] = 0.5*(Jkx(getpos(i,num_points,D),kappa[i],D,f) - kappa[i]); 
				double tmp = c1b * exp(-0.5*z1*z1*lb*kappa[i]-self_energy[i]);
				std::cout << self_energy[i] << std::endl; 
				error += pow(cation[i] - tmp, 2);
				cation[i] = tmp;
				anion[i] = c2b * exp(-0.5*z2*z2*lb*kappa[i]-self_energy[i]);
				psi[i] = 0;
			}
		}	


		for (int i=0;i<num_points;i++)
			output << getpos(i,num_points,D) << " " << cation[i] << " " << anion[i] << std::endl;


	}

}

	
		
