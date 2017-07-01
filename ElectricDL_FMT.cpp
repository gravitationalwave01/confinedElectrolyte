//This code computes the ion profile, electric potential, and total free energy for a ionic solution confined between
//two plates. Level of theory is Renormalized Gaussian Fluctuation Theory (described in ZG Wang's 2010 PRE paper)
//Free energy is computed using charging method (see R. Wang's 2015 JCP paper)
//Physical parameters of the space are hardcoded and can be changed under "initialization" section (just search for it)
//

# include <iostream>
# include <iomanip>
# include <fstream>
# include <string>
# include <math.h>
# include <cstring>
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
//# include <cstdio>
//# include <cstdlib>

using namespace std;

//define some useful structures
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// IonType is for setting up the basic features of the cations and anions in the system
struct IonType {
	int    Z[2];          //Z[0] means the valence of cations;
	double Dia[2];        //Dia[0] means the diameter of cations; (the unit is Angstrom)
};

// SysType is for setting up the system
struct SysType {
    double cs;            //the concentration of the salt; the the unit is M
	double sigma;         //the surface charge density; the unit is C/m^2
	double epsilonS;      //the relative dielectric constant for solvent
	double epsilonP;      //the relative dielectric constant for the planar surface
	double T;             //the temperature
};

// GridType is for setting up parameter for numerical method
struct GridType {
	int Nki;               //number of grids for the integration in k space
	int Nz1i;              //number of grids for the integration for z
    int Nz2i;              //number of grids for the integration for z'
};

// IteraType is for setting up the details for iteration
struct IteraType {
	int maxItera;         //maximum iteration steps
    double mixCoe;        //the mixture coefficient of the picard iteration
	double relTolPsi;     //maximum relative error for the potential
	double relTolRho;     //maximum relative error for the ion density
};

// ConsType is used to show the universal constant
struct ConsType {
	double epsilon0;      //the dielectric constant for the vacuum
	double e0;            //the elementary charge
	double kB;            //Bolzman constant
	double Pi;
	double NA;            //Avogodro number
};


//chemecal potential from hard sphere interaction
struct MuHards {
	double pos;
	double neg;
};


// external potential between hard wall and hard sphere
struct ExterHW {
	double pos;
	double neg;
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// define the global identifiers
const ConsType cons ={ 8.854187817e-12, 1.60217653e-19, 1.3806505e-23, 3.141592653589793, 6.02214129e23};

double *z = NULL;                  // z axis
double rho0[2];             //rho[0] means the bulk density of cations (scaled ones)
double kappa0, kappa02;     //inverse bulk Debye length (scaled one)
double *kappa = NULL;     //local screening length (inverse local Debye length) (scaled)
double *kappa2=NULL;
double *kappaI=NULL;
double *rho[2]={NULL,NULL};             //local density for cations and anions(scaled)
double *psi=NULL;
double *psi0=NULL;         //electrostatic potential (scaled)

double epsilS, epsilP;      //scaled permitivity
double D[2], D2[2];         //scaled diameters of cations and anions
double Radii[2], Radii2[2];
double sigma;               //scaled surface charge density
double lUnit;               //the lenght unit is the diameter of the smaller ions
double maxk;                //the upper bound for the integration in the k space (scaled by kappa0)
double maxz1;               //the upper bound for the integration in the k space (scaled by 1/kappa0)
double maxz2;               //the upper bound for the integration in the k space (scaled by 1/kappa0)
double dz1, dz2, dk, deta;        //the step lengths
double dz22;                //dz22=dz2*dz2
double boundary, coef;      //boundary is for solving poisson equation
double Pi1, Pi2, Pi4;
double const1, const2, const3, const4, const5;

int    Z1, Z2, Neta;              //the valences of cations and anions
int    rMinz1[2], rMinz2[2];//index of the left contact position
int    rMaxz1[2], rMaxz2[2];//index of the right contact position //Added by AV
int    dMinz1[2], maxNz1e, maxNz1ee;
int    Nk, Nz1, Nz2;

//end define global identifier
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//global function
double **WeightD(int iz, double **wd);    //function for calculate the weighted density with FMT method
double **DerivatPhi(double **Deri);  //function for calculate the derivative of Phi with respect to weighted density
MuHards  MuHardsphere(int iz, double **DPh); //function for calculate the local chemical potential from hard sphere
MuHards  MuHardbulk();                       //function for calculate the bulk chemical potential from hard sphere
ExterHW  ExternalHW(int iz);                 //External potential between hard wall and hard sphere
double   GreenB(double k);                   //function for solving the Green's function in bulk in Fourier space
double   Green(double k, int zn);            //function for solving the Green's function G(k,z,z) in Fourier space
void     Poisson(double dz);                 //function for solving the Poisson equation
double   GreenBEta(double k, double eta); // computing bulk greens function for charging method
double   GreenEta(double k, int zn, double eta); //computing full greens function for charging method
double   GreenFree(double k, int zn); //computing bulk greens function for 0 ions 
double   *DeltaPhi(double *deltaPhi); // increment 
double   **GreenFast(double **GreenKZn);
double   **GreenBEtaFast(double **GreenBEtaK);
double   ***GreenEtaFast(double ***GreenEtaKZ);
double   **GreenFreeFast(double **GreenFreeKZ);
double   *GreenBFast(double *GreenBK);

//end define global function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////




int main()
{
  int     zn1, zn2, kn, zn2i, Nz11e, Nz12e, Ncs, Ncs1, TPY, Nr1, Nr2;//index of the z grid, z' grid and k grid
  int     i1, i2, i3, i4, iCount, iterN, etan,rMinz2mi,rMinz1ma, Nz2t;
  double  z1, z2, inz, eta;                    //coordinate for z, z'
  double  Z12, Z22;                            //square of valency
  double  k1, k2, k3;
  double  err, errRho, errPsi, maxErrRho, maxErrPsi,iterMixCoe;
  double  errU, maxErrU;
  double  sumIntb, sumIntz;
  double  beta,bulkcs,bulkcs1;
  double  sUnit, vUnit;                       //the unit of area and volume, respectively.
  double  in1,in2,in3,in4;
  double  greenb, greenz,adsorp1,adsorp2,adsorp3,adsorp4;
  double  u10, u20, u1z, u2z, u1, u2, maxU1, maxU2;
  double  mu1Bulk, mu2Bulk, mu1, mu2,u1zNew,u2zNew;
  double  rho1New, rho2New, maxRho1, maxRho2,minRadii,maxRadii;
  double  **dPhi=NULL;
  double  **GreenKZn=NULL;
  double  **GreenBEtaK=NULL;
  double  ***GreenEtaKZ=NULL;
  double  **GreenFreeKZ=NULL;


  double  *rho11=NULL;
  double  *rho22=NULL;
  double  *zzz= NULL;
  double  *psi00=NULL;

  double  *GreenBK=NULL;
  double  *deltaPhi=NULL;
  double  coez, coerho, coee, sigmai, delta11, delta12, delta13 ,delta14;
  double  delta21, delta22, delta23 ,delta24, delta31, delta32, delta33 ,delta34;
  double  *deltaGz, *deltaG0, *deltaGzk, *deltaG0k, *G0k;
  double  *UU1, *UU2, *deltaG0eta, *deltaGzeta;
  double  *Ener1, *Ener3, *Ener4, *Ener5;
  double  sumIntz1,sumIntz2, temp1, temp2, temp3, EnergyTotal,sumIntz3,dz00,endTime;
  char    filename[30]={0},num[5]={0};
  clock_t start, finish;


  start=clock();
  //initialize
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  bool hard_sphere = true; //set to true if you want to allow hard sphere interactions
  IonType   ion  = {{1, 1}, {5.0, 5.0}};    //initialize the parameters of ions, the unit is Angstr
  SysType   sys  = {0.0001, 0, 80.0, 10000, 250.0}; //initialize the parameters of the system
  IteraType iter = {10000, 0.01, 1.0e-10, 1.0e-8};    //initialize the parameters of iteration
  GridType  grid = {200, 1000, 20000};        //initialize the integral grid numbers

//  IonType   ion  = {{1, 1}, {2.5, 2.5}};    //initialize the parameters of ions, the unit is Angstrom
//  SysType   sys  = {0.1, 1.60217653e-3, 80.0, 2.5, 298.0}; //initialize the parameters of the system
//  IteraType iter = {2000, 0.01, 1.0e-10, 1.0e-8};    //initialize the parameters of iteration
//  GridType  grid = {200, 1000, 20000};        //initialize the integral grid numbers
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


  MuHards   muHs;
  MuHards   muHs0;
  ExterHW   extHw;

  Nk    = grid.Nki;
  Nz1   = grid.Nz1i;
  Nz2t  = grid.Nz2i;
  Nz2   = Nz2t + 1;

  //typedef double ArrayN[grid.Nz1+1];

  //define the Array
  z        = new double[Nz1+1] ();
  zzz      = new double[Nz1+1] ();
  rho11    = new double[Nz1+1] ();
  rho22    = new double[Nz1+1] ();
  psi00    = new double[Nz1+1] ();
  kappa    = new double[Nz1+1] ();
  kappaI   = new double[Nz2+1] ();
  kappa2   = new double[Nz2+1] ();
  psi      = new double[Nz1+1] ();
  psi0     = new double[Nz1+1] ();
  deltaPhi = new double[Nz1+1] ();
  UU1      = new double[Nz1+1] ();
  UU2      = new double[Nz1+1] ();
  deltaGz  = new double[Nz1+1] ();
  deltaG0  = new double[Nz1+1] ();
  Ener1    = new double[Nz1+1] ();
  Ener3    = new double[Nz1+1] ();
  Ener4    = new double[Nz1+1] ();
  Ener5    = new double[Nz1+1] ();
  G0k      = new double[Nk+1] ();
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


  if(ion.Dia[0]<=ion.Dia[1])
  {
    lUnit = ion.Dia[0]*1.0e-10;            //here the unit is meter
    coez  = lUnit/1.0e-9;
    D[0]  = ion.Dia[0]/ion.Dia[0];         //here the unit is ion.Dia[0]
    D[1]  = ion.Dia[1]/ion.Dia[0];
    minRadii = D[0]*0.5;
    maxRadii = D[1]*0.5;
  }
  else
  {
    lUnit   = ion.Dia[1]*1.0e-10;         //here the unit is meter
    coez    = lUnit/1.0e-9;
    D[0]    = ion.Dia[0]/ion.Dia[1];      //here the unit is ion.Dia[1]
    D[1]    = ion.Dia[1]/ion.Dia[1];
	minRadii = D[1]*0.5;
	maxRadii = D[0]*0.5;
  }
  Radii[0]   = 0.5*D[0];
  Radii[1]   = 0.5*D[1];
  Radii2[0]  = Radii[0]*Radii[0];
  Radii2[1]  = Radii[1]*Radii[1];
  D2[0]      = D[0]*D[0];
  D2[1]      = D[1]*D[1];

  sUnit      = lUnit*lUnit;
  vUnit      = lUnit*lUnit*lUnit;
  beta       = 1.0/(sys.T*cons.kB);
  epsilS     = cons.epsilon0*sys.epsilonS*lUnit/(beta*cons.e0*cons.e0);
  epsilP     = cons.epsilon0*sys.epsilonP*lUnit/(beta*cons.e0*cons.e0);
  sigma      = sys.sigma*sUnit/cons.e0;
  coee       = cons.e0/sUnit;
  Z1         = ion.Z[0];
  Z2         = ion.Z[1];
  Z12        = Z1*Z1;
  Z22        = Z2*Z2;
  Pi1        = cons.Pi;
  Pi2        = Pi1 + Pi1;
  Pi4        = Pi2 + Pi2;
  coerho     = 1.0/(cons.NA*vUnit*1.0e3);
  deta       = 0.01;
  Neta       = int(1.0/deta + 1.0e-8);
  deltaG0eta = new double[Neta+1] ();
  deltaGzeta = new double[Neta+1] ();

  maxk       = 100.0;
  maxz1      = 20.0; //AV length of the domain in A
  maxz2      = 20.0;
  //maxz1      = 20.0*maxRadii;
  //maxz2      = maxz1*2.0;

  ofstream outFile4;
  outFile4<<setiosflags(ios::left)<<setiosflags(ios::fixed);
  outFile4.open("Contact_Density.txt");

  ofstream outFile5;
  outFile5<<setiosflags(ios::left)<<setiosflags(ios::fixed);
  outFile5.open("Adsorption_Potential.txt");



  strcpy(filename,"Total_Energy");
  //sprintf(num,"%d",Ncs);
  //strcat(filename,num);
  strcat(filename,".txt");
  ofstream outFile2;
  outFile2<<setiosflags(ios::left)<<setiosflags(ios::fixed);
  outFile2.open(filename);

  bulkcs     = sys.cs;
  //bulkcs1    = bulkcs - 0.01;
  TPY        = 0;
// for loop edited by AV on oct 27 2016 -- changed to iterate only once
//  for(Ncs=0; Ncs<=100; Ncs=Ncs+1)
for (Ncs = 0; Ncs<=3; Ncs++)
  {

  rho0[0]    = bulkcs*Z2*cons.NA*vUnit*1.0e3;
  rho0[1]    = bulkcs*Z1*cons.NA*vUnit*1.0e3;
  maxRho1    = rho0[0]*100.0;
  maxRho2    = rho0[1]*100.0;
  maxU1      = maxRho1/rho0[0];
  maxU2      = maxRho2/rho0[1];
  kappa0     = sqrt((Z12*rho0[0]+Z22*rho0[1])/(epsilS));
  kappa02    = kappa0*kappa0;
  dk         = maxk*kappa0/(Nk*1.0);

//this if loop checks that the coefficient C is greater than 0
//C is a coefficient of the discretization of the Greens Function differential equation
//JJ says that its doesn't matter what C is exactly, but the calculation is 
//sped up when C is larger than 1.
  if((1.0-maxk*kappa0*const1)<0.0)
  {
  	cerr<<"C[0] < 0"<<endl;
	cerr << "maxk: " << maxk << endl;
	cerr << "kappa0: " << kappa0 << endl;
	cerr << "const1: " << const1 << endl;  	
 	exit(0);	
  }
  const5 = 1.0e-12/(maxk*kappa0);

  //dz1        = maxz1/(Nz1*1.0);
  //dz2        = maxz2/(Nz2*1.0);
//  dz1        = maxz1/kappa0/(Nz1*1.0); //AV commented nad replaced with whats belowq
  dz1        = maxz1/(double)Nz1;
//  dz1        = minRadii/(double(int(minRadii/dz1+1.0e-8))); //commented out by AV
  dz2        = (Nz1*dz1)/(Nz2t*1.0);
  //dz2        = maxz2/kappa0/(Nz2t*1.0);
  //dz2        = minRadii/(double(int(minRadii/dz2+1.0e-8)));
  //if(dz1<1.0e-4) dz1=1.0e-4;
  //if(dz2<1.0e-5) dz2=1.0e-5;

  dz22       = dz2*dz2;

  coef       = 1.0*dz1*dz1/epsilS;
  const1     = 0.5*dz2*epsilP/epsilS;
  const2     = -dz2/epsilS;
  const3     = dz1/dz2;
  const4     = 0.5*dz2;


  rMinz1[0]  = int(Radii[0]/dz1 + 1.0e-8);
  rMaxz1[0]  = Nz1 - int(Radii[0]/dz1 + 1.0e-8); //Added by AV
  Nr1        = rMinz1[0];
  rMinz1[1]  = int(Radii[1]/dz1 + 1.0e-8);
  rMaxz1[1]  = Nz1 - int(Radii[1]/dz1 + 1.0e-8); //Added by AV
  Nr2        = rMinz1[1];
  rMinz2[0]  = int(Radii[0]/dz2 + 1.0e-8);
  rMinz2[1]  = int(Radii[1]/dz2 + 1.0e-8);
  dMinz1[0]  = int(D[0]/dz1 + 1.0e-8);
  dMinz1[1]  = int(D[1]/dz1 + 1.0e-8);
  Nz11e      = Nz1+dMinz1[0]+4;
  Nz12e      = Nz1+dMinz1[1]+4;

  if(rMinz2[0]>rMinz2[1])
  {
  	rMinz2mi=rMinz2[1];
  	rMinz1ma=rMinz1[0];
  }
  else
  {
  	rMinz2mi=rMinz2[0];
  	rMinz1ma=rMinz1[1];
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(Nz11e>=Nz12e)
  {
  	maxNz1e = Nz11e;
  }
  else
  {
  	maxNz1e = Nz12e;
  }

  if(rMinz1[0]>rMinz1[1])
  {
  	maxNz1ee = Nz1 + rMinz1[0] + 1;
  }
  else
  {
  	maxNz1ee = Nz1 + rMinz1[1] + 1;
  }

  rho[0] = new double[maxNz1e+1] ();
  rho[1] = new double[maxNz1e+1] ();



  dPhi = new double*[6] ();

  for(i1=0; i1<=5; i1++)
  {
  	dPhi[i1] = new double[maxNz1ee+1] ();
  }

  GreenKZn = new double*[Nk+1] ();
  for(i1=0; i1<=Nk; i1++)
  {
  	GreenKZn[i1] = new double[Nz1+1] ();
  }

  GreenBEtaK = new double*[Neta+1] ();
  for(i1=0; i1<=Neta; i1++)
  {
  	GreenBEtaK[i1] = new double[Nk+1] ();
  }

  GreenFreeKZ = new double*[Nk+1] ();
  for(i1=0; i1<=Nk; i1++)
  {
  	GreenFreeKZ[i1] = new double[Nz1+1] ();
  }

  GreenEtaKZ = new double**[Neta+1] ();
  for(i1=0; i1<=Neta; i1++)
  {
  	GreenEtaKZ[i1] = new double*[Nk+1] ();
  	for(i2=0; i2<=Nk; i2++)
  	{
  		GreenEtaKZ[i1][i2] = new double[Nz1+1] ();
	}
  }

  GreenBK = new double[Nk+1] ();
  /*calculate the Green function in bulk in real space
    the integration in the k space is calculated by the Simpon method
  */


  GreenBK    = GreenBFast(GreenBK);
  sumIntb    = 0.0;
  for(i1=1; i1<Nk-1; i1=i1+2)
  {
  	k1      = i1*dk;
  	k2      = k1 + dk;
  	k3      = k2 + dk;
  	sumIntb = sumIntb + GreenBK[i1]*k1/Pi2;
  	sumIntb = sumIntb + GreenBK[i1+1]*k2*2.0/Pi1;
  	sumIntb = sumIntb + GreenBK[i1+2]*k3/Pi2;
  }
  if((Nk-1)%2 != 0)
  {
  	k2      = Nk*dk;
  	k1      = k2 - dk;
  	sumIntb = sumIntb/3.0 + GreenBK[Nk-1]*k1/Pi4 + GreenBK[Nk]*k2/Pi4;
  	greenb  = sumIntb*dk;
  }
  else
  {
  	greenb = sumIntb*dk/3.0;
  }
  u10      = greenb*Z12*0.5;
  u20      = greenb*Z22*0.5;


  if (hard_sphere)
	muHs0   = MuHardbulk();
  else
  { 
	muHs0.pos = 0;	
	muHs0.neg = 0;	
  }

  //sigmai     = sigma;
  sigmai = sigma;
  //for(Ncs=1; Ncs<=3; Ncs=Ncs+1)
  boundary   = -sigmai*dz1/epsilS;

 //initialize the density and potential
  if(TPY==0)
  {

    if(Ncs == 0)
    {
        for(i1=0; i1<=maxNz1e; i1++)
        {
     	    if(i1<rMinz1[0])
            {
    	        rho[0][i1] = 0.0;
        	}
         	else if (i1>rMaxz1[0])
				rho[0][i1] = 0.0;
			else
        	{
        		rho[0][i1] = rho0[0];
        	}

         	if(i1<rMinz1[1])
         	{
        	    rho[1][i1] = 0.0;
        	}
         	else if (i1>rMaxz1[1])
				rho[1][i1] = 0.0;
        	else
        	{
        		rho[1][i1] = rho0[1];
        	}

        	if(i1<=Nz1)
        	{
        		kappa[i1] = sqrt((Z12*rho[0][i1]+Z22*rho[1][i1])/(epsilS));
        		z[i1]    = i1*dz1;
        	    psi[i1]  = 0.0;
        	}
        }

    }
    else
    {
    	Ncs1 = Ncs - 1;
    	strcpy(filename,"Density_Potential_");
        sprintf(num,"%d",Ncs1);
        strcat(filename,num);
        strcat(filename,".txt");
        ifstream inFile1;
        inFile1.open(filename);
        if(inFile1)
        {

            for(i1=0; i1<=Nz1; i1++)
            {
    			inFile1>>in1>>in2>>in3>>in4;
                zzz[i1]      = in1/coez;
                rho11[i1]    = in2/coerho;
                rho22[i1]    = in3/coerho;
                psi00[i1]    = in4;
         	}
    	//////////////////////////////////////////////////////////
        	dz00 = zzz[1]-zzz[0];
            for(i1=0; i1<=maxNz1e; i1++)
            {
            	if(i1<=Nz1)
            	{
            		z[i1]    = i1*dz1;
            		zn1      = int(z[i1]/dz00 + 1.0e-10);

            		if(zn1<Nz1)
            		{

						inz = z[i1]/dz00-zn1;
						if(i1<rMinz1[0])
                        {
                 	        rho[0][i1] = 0.0;
                        }
         				else if (i1>rMaxz1[0]) //AV
							rho[0][i1] = 0.0;
                     	else
                    	{
                    		rho[0][i1] = (rho11[zn1] + (rho11[zn1+1]-rho11[zn1])*inz)*(bulkcs/bulkcs1);
                    	}

                     	if(i1<rMinz1[1])
                    	{
                     	    rho[1][i1] = 0.0;
                     	}
         				else if (i1>rMaxz1[1]) //AV
							rho[1][i1] = 0.0;
                     	else
                    	{
                     		rho[1][i1] = (rho22[zn1] + (rho22[zn1+1]-rho22[zn1])*inz)*(bulkcs/bulkcs1);
                     	}
    					psi[i1]    = psi00[zn1] + (psi00[zn1+1]-psi00[zn1])*inz;
    					kappa[i1]  = sqrt((Z12*rho[0][i1]+Z22*rho[1][i1])/(epsilS));
    				}
    				else
    				{
    					rho[0][i1] = rho0[0];
    					rho[1][i1] = rho0[1];
    					psi[i1]    = 0.0;
    					kappa[i1]  = kappa0;
    				}
    	    	}
        		else
         		{
        			rho[0][i1]  = rho0[0];
                	rho[1][i1]  = rho0[1];
        		}
        	}
       /////////////////////////////////////////////////////////////////

        }
        else
        {
        	cerr<<"fail to open the initial density file"<<endl;
         	exit(0);
    	}
        inFile1.close();
        filename[0] = '\0';
    }
  }
  else
  {
    	strcpy(filename,"Density_Potential_tpy");
        strcat(filename,".txt");
        ifstream inFile0;
        inFile0.open(filename);
        if(inFile0)
        {

            for(i1=0; i1<=Nz1; i1++)
            {
    			inFile0>>in1>>in2>>in3>>in4;
                z[i1]      = in1;
                rho[0][i1] = in2;
                rho[1][i1] = in3;
                psi[i1]    = in4;
				kappa[i1]  = sqrt((Z12*rho[0][i1]+Z22*rho[1][i1])/(epsilS));
         	}

            for(i1=Nz1+1; i1<=maxNz1e; i1++)
            {
        		rho[0][i1]  = rho0[0];
                rho[1][i1]  = rho0[1];
        	}
        }
        else
        {
        	cerr<<"fail to open the initial density file"<<endl;
         	exit(0);
    	}
        inFile0.close();
        filename[0] = '\0';
		TPY = 0;
  }
 //end initialization


  iterN = 0;
  iterMixCoe  = iter.mixCoe;

  strcpy(filename,"Error_Detection_");
  sprintf(num,"%d",Ncs);
  strcat(filename,num);
  strcat(filename,".txt");
  ofstream outFile3;
  outFile3<<setiosflags(ios::left)<<setiosflags(ios::fixed);
  outFile3.open(filename);

  do
  {
	dPhi = DerivatPhi(dPhi);
	maxErrRho = 0.0;
 	maxErrPsi = 0.0;
 	maxErrU   = 0.0;
  //linear interpolation
    for(i1=1; i1<=(Nz2-1); i1++)
    {
	    if(i1<rMinz2mi)
		{
			kappaI[i1] = 0;
    	}
    	else
    	{
    		z2  = dz2*i1;
        	zn1 = int(z2/dz1 + 1.0e-8);
        	z1  = zn1*dz1;
        	inz = (z2 - z1)/dz1;
        	if(zn1<Nz1)
        	{
        		if(zn1==(rMinz1ma-1))
        		{
        			kappaI[i1] = kappa[zn1] + (kappa[zn1]-kappa[zn1-1])*inz;
				}
        		else
        		{
        			kappaI[i1] = kappa[zn1] + (kappa[zn1+1]-kappa[zn1])*inz;
				}
        	}
        	else
        	{
        		kappaI[i1] = kappa0;
        	}

		}

    	kappa2[i1] = kappaI[i1]*kappaI[i1];
    }

 	GreenKZn = GreenFast(GreenKZn);


 	for(i1=1; i1<=Nz1; i1++)
 	{
 		sumIntz = 0.0;
 		//zn2     = int(i1*dz1/dz2 + 0.5);
 		for(i2=1; i2<Nk-1; i2=i2+2)
 		{
 			k1      = i2*dk;
 			k2      = k1 + dk;
 			k3      = k2 + dk;
 			sumIntz = sumIntz + GreenKZn[i2][i1]*k1/Pi2;
 			sumIntz = sumIntz + GreenKZn[i2+1][i1]*k2*2.0/Pi1;
 			sumIntz = sumIntz + GreenKZn[i2+2][i1]*k3/Pi2;
		}
		if((Nk-1)%2 != 0)
		{
			k2      = Nk*dk;
			k1      = k2 - dk;
			sumIntz = sumIntz/3.0 + GreenKZn[Nk-1][i1]*k1/Pi4 + GreenKZn[Nk][i1]*k2/Pi4;
			greenz  = sumIntz*dk;
		}
		else
		{
			greenz  = sumIntz*dk/3.0;
		}


		u1z  = greenz*Z12*0.5;
		u2z  = greenz*Z22*0.5;


		//calculate the local chemical potential
		if (hard_sphere)
		{
			muHs  = MuHardsphere(i1,dPhi); //AV turn off hard shell potential here
			extHw = ExternalHW(i1); //removing collisions between hard shell and the wall
		}
		else 
		{
			muHs.pos = 0;
			muHs.neg = 0;
			extHw.pos = 0; 	
			extHw.neg = 0;
		}
		mu1   = muHs.pos - muHs0.pos + extHw.pos;
		mu2   = muHs.neg - muHs0.neg + extHw.neg;

		rho1New = rho0[0]*exp(u10-u1z-Z1*psi[i1]-mu1);
		rho2New = rho0[1]*exp(u20-u2z+Z2*psi[i1]-mu2);

		if((rho1New>=maxRho1)||(rho1New<=1.0e-30))
		{
			rho1New = 0.0;
		}
		if((rho2New>=maxRho2)||(rho2New<=1.0e-30))
		{
			rho2New = 0.0;
		}

		errRho = fabs((rho1New-rho[0][i1])/(rho0[0]));
		err    = fabs((rho2New-rho[1][i1])/(rho0[1]));
		if(errRho<err)
		{
			errRho = err;
		}
		if(errRho>maxErrRho)
		{
			maxErrRho = errRho;
		}


		if(maxErrRho<1.0e-6)
		{
			iterMixCoe = 0.5;
		}
		else if(maxErrRho<1.0e-5)
		{
			iterMixCoe = 0.4;
		}
		else if(maxErrRho<1.0e-4)
		{
			iterMixCoe = 0.2;
		}
		else if(maxErrRho<1.0e-3)
		{
			iterMixCoe = 0.1;
		}
		else if(maxErrRho<1.0e-2)
		{
			iterMixCoe = 0.05;
		}
		else
		{
			iterMixCoe = iter.mixCoe;
		}

		if(maxErrRho>100) iterMixCoe = 0.001;



		rho[0][i1] = (1.0 - iterMixCoe)*rho[0][i1] + iterMixCoe*rho1New;
		rho[1][i1] = (1.0 - iterMixCoe)*rho[1][i1] + iterMixCoe*rho2New;

		kappa[i1] = sqrt((Z12*rho[0][i1]+Z22*rho[1][i1])/(epsilS));

/*		//update the interpolation of kappaI
		z1 = i1*dz1;
		for(i3=(int((z1-dz1)/dz2)+1); i3<=int((z1+dz1)/dz2); i3++)
		{
			z2  = dz2*i3;
			zn1 = int(z2/dz1);
			inz = z2 - zn1*dz1;
			if(zn1<Nz1)
			{
				kappaI[i3] = kappa[zn1] + (kappa[zn1+1]-kappa[zn1])*inz;
			}
			else
			{
				kappaI[i3] = kappa0;
			}
			kappa2[i3] = kappaI[i3]*kappaI[i3];
		}

		*/


	}

	Poisson(dz1);
	for(i1=0; i1<=Nz1; i1++)
	{
		//psi[i1] = (1.0 - iter.mixCoe)*psi0[i1] + iter.mixCoe*psi[i1];
		errPsi = fabs(psi[i1]-psi0[i1]);
		if(errPsi>maxErrPsi)
		{
			maxErrPsi = errPsi;
		}
	}


	outFile3<<"iter="<<iterN<<"\t"<<"maxErrRho="<<maxErrRho<<"\t"<<"maxErrPsi=";
	outFile3<<maxErrPsi<<endl;
 	iterN = iterN + 1;


 	finish  = clock();
 	endTime = (double)(finish-start)/CLOCKS_PER_SEC;
 	endTime = endTime/3600.0;
	//std::cout << "total time was " << endTime << std::endl;
 	if((iterN==iter.maxItera)||(endTime>47.5))
 	{

 		strcpy(filename,"Density_Potential_tpy");
        strcat(filename,".txt");
        ofstream outFile0;
        outFile0.precision(12);
        outFile0<<setiosflags(ios::left)<<setiosflags(ios::fixed);
        outFile0.open(filename);
        for(i1=0; i1<=Nz1; i1++)
        {
          outFile0<<z[i1]<<"\t"<<rho[0][i1]<<"\t"<<rho[1][i1]<<"\t"<<psi[i1]<<endl;
        }
        outFile0.close();

        outFile0<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
        filename[0] = '\0';
	}
  }
  while(((maxErrPsi>iter.relTolPsi)||(maxErrRho>iter.relTolRho))&&(iterN<iter.maxItera));

  if((maxErrPsi>iter.relTolPsi)||(maxErrRho>iter.relTolRho))
  {
  	cerr<<"exceed the max iteration"<<endl;
    exit(0);
  }

  outFile3.close();
  outFile3<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
  filename[0] = '\0';

  outFile4<<bulkcs<<"\t"<<rho[0][Nr1]<<"\t"<<rho[1][Nr2]<<endl;

  strcpy(filename,"Density_Potential_");
  sprintf(num,"%d",Ncs);
  strcat(filename,num);
  strcat(filename,".txt");
  ofstream outFile1;
  outFile1.precision(12);
  outFile1<<setiosflags(ios::left)<<setiosflags(ios::fixed);
  outFile1.open(filename);
  for(i1=0; i1<=Nz1; i1++)
  {
    outFile1<<z[i1]*coez<<"\t"<<rho[0][i1]*coerho<<"\t"<<rho[1][i1]*coerho<<"\t"<<psi[i1]<<endl;
  }
  outFile1.close();

  outFile1<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
  filename[0] = '\0';

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Calculate the adsorption quantity rho0[0],rho0[1]
  sumIntz = 0.0;
  sumIntz1= 0.0;
  sumIntz2= 0.0;
  sumIntz3= 0.0;


  for(i1=0; i1<Nz1-1; i1=i1+2)
  {

 	sumIntz  = sumIntz + rho[0][i1]-rho0[0];
 	sumIntz  = sumIntz + 4.0*(rho[0][i1+1]-rho0[0]);
 	sumIntz  = sumIntz + rho[0][i1+2]-rho0[0];

	sumIntz1 = sumIntz1 + rho[1][i1]-rho0[1];
 	sumIntz1 = sumIntz1 + 4.0*(rho[1][i1+1]-rho0[1]);
 	sumIntz1 = sumIntz1 + rho[1][i1+2]-rho0[1];

 	sumIntz2 = sumIntz2 + rho[0][i1]-rho0[0]+rho[1][i1]-rho0[1];
 	sumIntz2 = sumIntz2 + 4.0*(rho[0][i1+1]-rho0[0]+rho[1][i1+1]-rho0[1]);
 	sumIntz2 = sumIntz2 + rho[0][i1+2]-rho0[0]+rho[1][i1+2]-rho0[1];

	sumIntz3 = sumIntz3 + Z1*rho[0][i1]-Z2*rho[1][i1];
 	sumIntz3 = sumIntz3 + 4.0*(Z1*rho[0][i1+1]-Z2*rho[1][i1+1]);
 	sumIntz3 = sumIntz3 + Z1*rho[0][i1+2]-Z2*rho[1][i1+2];

  }
  if(Nz1%2 != 0)
  {
	sumIntz  = sumIntz/3.0 + (rho[0][Nz1-1]+rho[0][Nz1])*0.5-rho0[0];
	adsorp1  = sumIntz*dz1;

	sumIntz1 = sumIntz1/3.0 + (rho[1][Nz1-1]+rho[1][Nz1])*0.5-rho0[1];
	adsorp2  = sumIntz1*dz1;

	sumIntz2 = sumIntz2/3.0 + (rho[0][Nz1-1]+rho[0][Nz1]+rho[1][Nz1-1]+rho[1][Nz1])*0.5-rho0[0]-rho0[1];
	adsorp3  = sumIntz2*dz1;

	sumIntz3 = sumIntz3/3.0 + (Z1*(rho[0][Nz1-1]+rho[0][Nz1])-Z2*(rho[1][Nz1-1]+rho[1][Nz1]))*0.5;
	adsorp4  = sumIntz3*dz1;

  }
  else
  {
	adsorp1  = sumIntz*dz1/3.0;
	adsorp2  = sumIntz1*dz1/3.0;
	adsorp3  = sumIntz2*dz1/3.0;
	adsorp4  = sumIntz3*dz1/3.0;
  }

  outFile5<<bulkcs<<"\t"<<adsorp1<<"\t"<<adsorp2<<"\t"<<adsorp3<<"\t"<<adsorp4<<"\t"<<psi[0]<<endl;

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //calculate the total energy per Area

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////


  GreenBEtaK  = GreenBEtaFast(GreenBEtaK);
  GreenEtaKZ  = GreenEtaFast(GreenEtaKZ);
  GreenFreeKZ = GreenFreeFast(GreenFreeKZ);
  //GreenBK     = GreenBFast(GreenBK);

  for(i1=1; i1<=Nz1; i1++)
  {
 	//zn2 = int(i1*dz1/dz2 + 0.5);

	for(i2=0; i2<=Neta; i2=i2+1)
	{
		//eta = i2*deta;
		sumIntz1 = 0.0;
		sumIntz2 = 0.0;
		for(i3=1; i3<Nk-1; i3=i3+2)
		{
			k1        = i3*dk;
 			k2        = k1 + dk;
 			k3        = k2 + dk;
 			if(i2 == 0)
 			{
 			    temp1 = GreenFreeKZ[i3][i1];
			    temp2 =	GreenFreeKZ[i3+1][i1];
			    temp3 =	GreenFreeKZ[i3+2][i1];
		    }
			else
			{
			    temp1 =	GreenEtaKZ[i2][i3][i1];
			    temp2 =	GreenEtaKZ[i2][i3+1][i1];
			    temp3 =	GreenEtaKZ[i2][i3+2][i1];
			}

 			sumIntz1 = sumIntz1 + (temp1-GreenKZn[i3][i1])*k1/Pi2;
    		sumIntz1 = sumIntz1 + (temp2-GreenKZn[i3+1][i1])*k2*2.0/Pi1;
    		sumIntz1 = sumIntz1 + (temp3-GreenKZn[i3+2][i1])*k3/Pi2;

			if(i2 == 0)
 			{
		        temp1 = GreenFreeKZ[i3][i1];
			    temp2 =	GreenFreeKZ[i3+1][i1];
			    temp3 =	GreenFreeKZ[i3+2][i1];
			}
			else
			{
			    temp1 =	GreenBEtaK[i2][i3];
			    temp2 =	GreenBEtaK[i2][i3+1];
			    temp3 =	GreenBEtaK[i2][i3+2];
			}

    		sumIntz2 = sumIntz2 + (temp1-GreenBK[i3])*k1/Pi2;
    		sumIntz2 = sumIntz2 + (temp2-GreenBK[i3+1])*k2*2.0/Pi1;
    		sumIntz2 = sumIntz2 + (temp3-GreenBK[i3+2])*k3/Pi2;
		}

		if((Nk-1)%2 != 0)
	    {
     		k2      = Nk*dk;
    		k1      = k2 - dk;
		    if(i2 == 0)
 			{
 			    temp1 = GreenFreeKZ[Nk-1][i1];
			    temp2 =	GreenFreeKZ[Nk][i1];
			}
			else
			{
			    temp1 =	GreenEtaKZ[i2][Nk-1][i1];
			    temp2 =	GreenEtaKZ[i2][Nk][i1];
			}
    		sumIntz1 = sumIntz1/3.0 + (temp1-GreenKZn[Nk-1][i1])*k1/Pi4 + \
			           (temp2-GreenKZn[Nk][i1])*k2/Pi4;
    		deltaGzeta[i2]  = sumIntz1*dk;

    		if(i2 == 0)
 			{
 			    temp1 = GreenFreeKZ[Nk-1][i1];
			    temp2 =	GreenFreeKZ[Nk][i1];
			}
			else
			{
			    temp1 =	GreenBEtaK[i2][Nk-1];
			    temp2 =	GreenBEtaK[i2][Nk];
			}
    		sumIntz2 = sumIntz2/3.0 + (temp1-GreenBK[Nk-1])*k1/Pi4 + \
			           (temp2-GreenBK[Nk])*k2/Pi4;
    		deltaG0eta[i2]  = sumIntz2*dk;
    	}
    	else
    	{
    		deltaGzeta[i2]  = sumIntz1*dk/3.0;
    		deltaG0eta[i2]  = sumIntz2*dk/3.0;
    	}

	}


	sumIntz1 = 0.0;
 	sumIntz2 = 0.0;
 	for(i2=0; i2<Neta-1; i2=i2+2)
 	{

 		sumIntz1 = sumIntz1 + deltaGzeta[i2];
 		sumIntz1 = sumIntz1 + 4.0*deltaGzeta[i2+1];
 		sumIntz1 = sumIntz1 + deltaGzeta[i2+2];

 		sumIntz2 = sumIntz2 + deltaG0eta[i2];
 		sumIntz2 = sumIntz2 + 4.0*deltaG0eta[i2+1];
 		sumIntz2 = sumIntz2 + deltaG0eta[i2+2];
	}
	if(Neta%2 != 0)
	{
		sumIntz1 = sumIntz1/3.0 + 0.5*(deltaGzeta[Neta-1]+deltaGzeta[Neta]);
		deltaGz[i1] = sumIntz1*deta;
		sumIntz2 = sumIntz2/3.0 + 0.5*(deltaG0eta[Neta-1]+deltaG0eta[Neta]);
		deltaG0[i1] = sumIntz2*deta;
	}
	else
	{
		deltaGz[i1] = sumIntz1*deta/3.0;
		deltaG0[i1] = sumIntz2*deta/3.0;
	}

	sumIntz = 0;

 	for(i2=1; i2<Nk-1; i2=i2+2)
 	{
 		k1      = i2*dk;
 		k2      = k1 + dk;
 		k3      = k2 + dk;
 		sumIntz = sumIntz + GreenKZn[i2][i1]*k1/Pi2;
 		sumIntz = sumIntz + GreenKZn[i2+1][i1]*k2*2.0/Pi1;
 		sumIntz = sumIntz + GreenKZn[i2+2][i1]*k3/Pi2;
	}
	if((Nk-1)%2 != 0)
	{
		k2      = Nk*dk;
		k1      = k2 - dk;
		sumIntz = sumIntz/3.0 + GreenKZn[Nk-1][i1]*k1/Pi4 + GreenKZn[Nk][i1]*k2/Pi4;
		greenz  = sumIntz*dk;
	}
	else
	{
		greenz  = sumIntz*dk/3.0;
	}

	u1z         = greenz*Z12*0.5;
	u2z         = greenz*Z22*0.5;
	UU1[i1]     = u1z;
	UU2[i1]     = u2z;

	//calculate the local chemical potential
  }
  

  mu1Bulk = log(rho0[0]) + muHs0.pos + u10;
  mu2Bulk = log(rho0[1]) + muHs0.neg + u20;


  for(i1=1; i1<=Nz1; i1=i1+1)
  {
	delta11=0.0;
	if(i1<rMinz1[0])
	{
		delta11 = delta11 - rho0[0]*(log(rho0[0])-1.0);
	}
	else
	{
		delta11 = delta11 + rho[0][i1]*(log(rho[0][i1])-1.0) - \
		          rho0[0]*(log(rho0[0])-1.0);
	}

	if(i1<rMinz1[1])
	{
		delta11 = delta11 - rho0[1]*(log(rho0[1])-1.0);
	}
	else
	{
		delta11 = delta11 + rho[1][i1]*(log(rho[1][i1])-1.0) - \
		          rho0[1]*(log(rho0[1])-1.0);
	}
	Ener1[i1] = delta11;
	Ener3[i1] = rho[0][i1]*(mu1Bulk-UU1[i1]) - rho0[0]*(mu1Bulk-u10) +\
	            rho[1][i1]*(mu2Bulk-UU2[i1]) - rho0[1]*(mu2Bulk-u20);
	Ener4[i1] = 0.5*psi[i1]*(Z1*rho[0][i1]-Z2*rho[1][i1]);
	Ener5[i1] = (deltaGz[i1]*kappa[i1]*kappa[i1] -\
	            deltaG0[i1]*kappa02)*0.5*epsilS;

  }
  deltaPhi = DeltaPhi(deltaPhi);

  sumIntz = 0.0;
  for(i1=1; i1<Nz1-1; i1=i1+2)
  {
    sumIntz = sumIntz + Ener1[i1] + deltaPhi[i1] - Ener3[i1] +\
              Ener4[i1] + Ener5[i1];
    sumIntz = sumIntz + (Ener1[i1+1] + deltaPhi[i1+1] - Ener3[i1+1] +\
              Ener4[i1+1] + Ener5[i1+1])*4.0;
    sumIntz = sumIntz + Ener1[i1+2] + deltaPhi[i1+2] - Ener3[i1+2] +\
              Ener4[i1+2] + Ener5[i1+2];
  }
  if((Nz1-1)%2 != 0)
  {
	sumIntz = sumIntz/3.0 + 0.5*(Ener1[Nz1] + deltaPhi[Nz1] - Ener3[Nz1] +\
              Ener4[Nz1] + Ener5[Nz1] +Ener1[Nz1-1] + deltaPhi[Nz1-1] - Ener3[Nz1-1] +\
              Ener4[Nz1-1] + Ener5[Nz1-1]);

	EnergyTotal = sumIntz*dz1;
  }
  else
  {
	EnergyTotal = sumIntz*dz1/3.0;
  }

  outFile2<<sigmai<<'\t'<<bulkcs<<'\t'<<EnergyTotal<<endl;

  //sigmai = sigmai + sigma;
  bulkcs1=bulkcs;
	
  //bulkcs = bulkcs + 0.02; //commented out by AV on Oct 28 2016. replaced with following:
 //epsilP *= 2;
	bulkcs *= 10;
 /* if(Ncs==0)
  {
  	sigmai = 5.0*sigma;
  }
  else if(Ncs==1)
  {
  	sigmai = 10.0*sigma;
  }
  else if(Ncs==2)
  {
  	sigmai = 20.0*sigma;
  }
  else
  {
  	sigmai = 25.0*sigma;
  }*/

}
  outFile2.close();
  outFile2<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
  filename[0] = '\0';
  //delete the memory

  for(i1=0; i1<=5; i1++)
  {
 	delete [] dPhi[i1];
  }
  delete [] dPhi;


  for(i1=0; i1<=Nk; i1++)
  {
  	delete [] GreenKZn[i1];
  }
  delete [] GreenKZn;

  for(i1=0; i1<=Neta; i1++)
  {
  	delete [] GreenBEtaK[i1];
  }
  delete [] GreenBEtaK;


  for(i1=0; i1<=Nk; i1++)
  {
  	delete [] GreenFreeKZ[i1];
  }
  delete [] GreenFreeKZ;


  for(i1=0; i1<=Neta; i1++)
  {
  	for(i2=0; i2<=Nk; i2++)
  	{
  		delete [] GreenEtaKZ[i1][i2];
	}
	delete [] GreenEtaKZ[i1];
  }
  delete [] GreenEtaKZ;
  delete [] GreenBK;
  delete [] rho[0];
  delete [] rho[1];


  outFile4.close();
  outFile4<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);
  outFile5.close();
  outFile5<<resetiosflags(ios::left)<<resetiosflags(ios::fixed);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////



  delete [] z;
  delete [] kappa;
  delete [] kappaI;
  delete [] kappa2;
  delete [] psi;
  delete [] psi0;
  delete [] deltaPhi;
  delete [] UU1;
  delete [] UU2;
  delete [] deltaGz;
  delete [] deltaG0;
  delete [] G0k;
  delete [] deltaG0eta;
  delete [] deltaGzeta;
  delete [] Ener1;
  delete [] Ener3;
  delete [] Ener4;
  delete [] Ener5;
  delete [] zzz;
  delete [] rho11;
  delete [] rho22;
  delete [] psi00;


  return 0;
}



//using finite difference method and chasing method
double GreenB(double k)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green;
	double greenk, k2, bi;
	int    i1, zn2;

	A     = new double[Nz2+1];
	B     = new double[Nz2+1];
	C     = new double[Nz2+1];
	E     = new double[Nz2+1];
	L     = new double[Nz2+1];
	R     = new double[Nz2];
	Y     = new double[Nz2+1];
	green = new double[Nz2+1];

	A[0]   = 0.0;
	B[0]   = -1.0;
	C[0]   = 1.0;
	E[0]   = 0.0;
	A[Nz2] = 1.0;
	B[Nz2] = -1.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	k2     = k*k;
	bi     = 2.0 + (kappa02+k2)*dz22;
	zn2    = Nz2/2;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		B[i1] = -bi;
		C[i1] = 1.0;
		E[i1] = 0.0;
	}

	E[zn2]    = -dz2/epsilS;

	//chasing method
	L[0] = B[0];
	Y[0] = E[0]/L[0];
	for(i1=1; i1<=Nz2; i1++)
	{
		R[i1-1] = C[i1-1]/L[i1-1];
		L[i1]   = B[i1] - A[i1]*R[i1-1];
		Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
	}

	green[Nz2] = Y[Nz2];
	for(i1=Nz2-1; i1>=zn2; i1--)
	{
		green[i1] = Y[i1] - R[i1]*green[i1+1];
	}

	greenk = green[zn2];

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] green;

	return greenk;
}




//using finite difference method and chasing method
double Green(double k, int zn)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green;
	double greenkz, k2, bi;
	int    i1;

	A     = new double[Nz2+1];
	B     = new double[Nz2+1];
	C     = new double[Nz2+1];
	E     = new double[Nz2+1];
	L     = new double[Nz2+1];
	R     = new double[Nz2];
	Y     = new double[Nz2+1];
	green = new double[Nz2+1];

	A[0]   = 0.0;
	B[0]   = -(1.0 + 0.5*k*(epsilP/epsilS)*dz2);
	C[0]   = 1.0 - 0.5*k*(epsilP/epsilS)*dz2;
	E[0]   = 0.0;
	A[Nz2] = 1.0 - 0.5*k*(epsilP/epsilS)*dz2;
	B[Nz2] = -(1.0 + 0.5*k*(epsilP/epsilS)*dz2);
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	k2     = k*k;
	bi     = 2.0 + k2*dz22;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		C[i1] = 1.0;
		E[i1] = 0.0;
		if(rMinz2[0]>rMinz2[1])
		{
			if(i1<rMinz2[1])
			{
				B[i1] = -bi;
			}
			else
			{
				B[i1] = -(2.0 + (kappa2[i1] + k2)*dz22);
			}
		}
		else
		{
			if(i1<rMinz2[0])
			{
				B[i1] = -bi;
			}
			else
			{
				B[i1] = -(2.0 + (kappa2[i1] + k2)*dz22);
			}
		}

	}

	E[zn]  = -dz2/epsilS;

	//chasing method
	L[0] = B[0];
	Y[0] = E[0]/L[0];
	for(i1=1; i1<=Nz2; i1++)
	{
		R[i1-1] = C[i1-1]/L[i1-1];
		L[i1]   = B[i1] - A[i1]*R[i1-1];
		Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
	}

	green[Nz2] = Y[Nz2];
	for(i1=Nz2-1; i1>=zn; i1--)
	{
		green[i1] = Y[i1] - R[i1]*green[i1+1];
	}

	greenkz = green[zn];

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] green;

	return greenkz;
}



//using finite difference method and chasing method
void Poisson(double dz)
{
	double *A, *B, *C, *E, *L, *R, *Y;
	double dz12, bi;
	int    i1;


	A     = new double[Nz1+1] ();
	B     = new double[Nz1+1] ();
	C     = new double[Nz1+1] ();
	E     = new double[Nz1+1] ();
	L     = new double[Nz1+1] ();
	R     = new double[Nz1] ();
	Y     = new double[Nz1+1] ();

	for(i1=0; i1<=Nz1; i1++)
	{
		psi0[i1] = psi[i1];
	}

	A[0]   = 0.0;
	B[0]   = -1.0;
	C[0]   = 1.0;
	E[0]   = boundary;
	A[Nz1] = 1.0;//AV changed from 0
	B[Nz1] = -1.0;//AV changed from 1.0
	C[Nz1] = 0;
	E[Nz1] = boundary; //AV changed from 0

	dz12   = dz1*dz1;
	bi     = dz12*kappa02;
	for(i1=1; i1<=(Nz1-1); i1++)
	{
		A[i1] = 1.0;
		B[i1] = -2.0-bi;
		C[i1] = 1.0;
		E[i1] = (Z2*rho[1][i1]-Z1*rho[0][i1])*coef - bi*psi0[i1]; //only thing that changes when solving for a delta function
	}


	//chasing method
	L[0] = B[0];
	Y[0] = E[0]/L[0];
	for(i1=1; i1<=Nz1; i1++)
	{
		R[i1-1] = C[i1-1]/L[i1-1];
		L[i1]   = B[i1] - A[i1]*R[i1-1];
		Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
	}

	psi[Nz1] = Y[Nz1];
	for(i1=Nz1-1; i1>=0; i1--)
	{
		psi[i1] = Y[i1] - R[i1]*psi[i1+1];
	}

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;

}



ExterHW ExternalHW(int iz)
{
	//AV turning off hard shell interactions
	ExterHW extp;
	if(iz<rMinz1[0] || iz > rMaxz1[0]) //AV added OR statement
	{
		extp.pos = 1.0e40;
	}
	else
	{
		extp.pos = 0.0;
	}

	if(iz<rMinz1[1] || iz > rMaxz1[1]) //AV added OR statement
	{
		extp.neg = 1.0e40;
	}
	else
	{
		extp.neg = 0.0;
	}

    return extp;
}


double **WeightD(int iz, double **NI)
{
	double   N0I[2], N1I[2], N2I[2], N3I[2], NV1I[2], NV2I[2];
	double   zz, minz, maxz, zz1, zz2, zz3;
	double   coe21, coe22, coe23, coe31, coe32, coe33;
	int      minzn, maxzn, i1, i2, i3;

	zz   = iz*dz1;
	for(i1=0; i1<=1; i1++)
	{
		minz  = zz - Radii[i1];
		maxz  = zz + Radii[i1];

		minzn = int(minz/dz1 + 1.0e-8);
		maxzn = int(maxz/dz1 + 1.0e-8);
		if(minzn<rMinz1[i1])
		{
			minzn = rMinz1[i1];
		}
		if(maxzn> (maxNz1e-1))
		{
			maxzn = maxNz1e-1;
		}
		N2I[i1]  = 0.0;
		N3I[i1]  = 0.0;
		NV2I[i1] = 0.0;
		for(i2=minzn; i2<maxzn-1; i2=i2+2)
		{
  	        zz1      = i2*dz1;
  	        zz2      = zz1 + dz1;
          	zz3      = zz2 + dz1;
          	coe21    = Radii2[i1] - (zz-zz1)*(zz-zz1);
          	coe22    = Radii2[i1] - (zz-zz2)*(zz-zz2);
          	coe23    = Radii2[i1] - (zz-zz3)*(zz-zz3);
          	coe31    = zz - zz1;
          	coe32    = zz - zz2;
          	coe33    = zz - zz3;

  	        N2I[i1]  = N2I[i1] + rho[i1][i2];
  	        N2I[i1]  = N2I[i1] + 4.0*rho[i1][i2+1];
  	        N2I[i1]  = N2I[i1] + rho[i1][i2+2];

  	        N3I[i1]  = N3I[i1] + rho[i1][i2]*coe21;
  	        N3I[i1]  = N3I[i1] + 4.0*rho[i1][i2+1]*coe22;
  	        N3I[i1]  = N3I[i1] + rho[i1][i2+2]*coe23;

  	        NV2I[i1] = NV2I[i1] + rho[i1][i2]*coe31;
  	        NV2I[i1] = NV2I[i1] + 4.0*rho[i1][i2+1]*coe32;
  	        NV2I[i1] = NV2I[i1] + rho[i1][i2+2]*coe33;
        }
        if(((maxzn-minzn)%2 != 0)&&(maxzn != minzn))
        {
  	        zz2      = maxzn*dz1;
  	        zz1      = zz2 - dz1;
  	        coe22    = Radii2[i1] - (zz-zz2)*(zz-zz2);
  	        coe21    = Radii2[i1] - (zz-zz1)*(zz-zz1);
  	        coe32    = zz - zz2;
  	        coe31    = zz - zz1;

  	        N2I[i1]  = N2I[i1]/3.0 + rho[i1][maxzn-1]/2.0 + rho[i1][maxzn]/2.0;

  	        N3I[i1]  = N3I[i1]/3.0 + rho[i1][maxzn-1]*coe21/2.0 + rho[i1][maxzn]*coe22/2.0;

  	        NV2I[i1] = NV2I[i1]/3.0 + rho[i1][maxzn-1]*coe31/2.0 + rho[i1][maxzn]*coe32/2.0;

  	        N2I[i1]  = N2I[i1]*Pi1*D[i1]*dz1;
            N3I[i1]  = N3I[i1]*Pi1*dz1;
            NV2I[i1] = NV2I[i1]*Pi2*dz1;
        }
        else
        {
            N2I[i1]  = N2I[i1]*Pi1*D[i1]*dz1/3.0;
            N3I[i1]  = N3I[i1]*Pi1*dz1/3.0;
            NV2I[i1] = NV2I[i1]*Pi2*dz1/3.0;
		}

        N0I[i1]  = N2I[i1]/(Pi1*D2[i1]);
		N1I[i1]  = N2I[i1]/(Pi2*D[i1]);
		NV1I[i1] = NV2I[i1]/(Pi2*D[i1]);

		NI[i1][0] = N0I[i1];
     	NI[i1][1] = N1I[i1];
    	NI[i1][2] = N2I[i1];
    	NI[i1][3] = N3I[i1];
	    NI[i1][4] = NV1I[i1];
     	NI[i1][5] = NV2I[i1];

	}

    return NI;
}




double **DerivatPhi(double **dPhi)
{
	double **NI;
	double N0, N1, N2, N3, NV1, NV2, Deri[6];
	double N331, N332, N333, N22,N32, NV22;
	int    i1, i2, i3;

	NI = new double*[2];
	for(i1=0; i1<=1; i1++)
	{
		NI[i1] = new double[6];
	}


	for(i1=0; i1<=maxNz1ee; i1++)
	{

		NI = WeightD(i1, NI);

		N0  = NI[0][0]  + NI[1][0];
		N1  = NI[0][1]  + NI[1][1];
		N2  = NI[0][2]  + NI[1][2];
		N3  = NI[0][3]  + NI[1][3];
		NV1 = NI[0][4]  + NI[1][4];
		NV2 = NI[0][5]  + NI[1][5];

		if(N3>=1.0) N3=0.75;

		if(N3<1.0e-10)
		{
			for(i2=0; i2<=5; i2++)
			{
			   dPhi[i2][i1] = 0.0;
			}
		}
		else
		{
			N331    = 1.0-N3;
			N332    = N331*N331;
			N333    = N332*N331;
			N22     = N2*N2;
			N32     = N3*N3;
			NV22    = NV2*NV2;
			Deri[0] = -log(N331);
			Deri[1] = N2/N331;
			Deri[2] = N1/N331+(log(N331)/N3+1.0/N332)*(N22-NV22)/(12.0*Pi1*N3);
			Deri[3] = N0/N331+(N1*N2-NV1*NV2)/N332-(log(N331)/(18.0*Pi1*N32*N3)\
			          +1.0/(36.0*Pi1*N32*N331)+(1.0-3.0*N3)/(36.0*Pi1*N32*N333))\
			          *(N22*N2-3.0*N2*NV22);
			Deri[4] = -NV2/N331;
			Deri[5] = -NV1/N331-(log(N331)/N3+1.0/N332)*N2*NV2/(6.0*Pi1*N3);

			dPhi[0][i1] = Deri[0];
			dPhi[1][i1] = Deri[1];
			dPhi[2][i1] = Deri[2];
			dPhi[3][i1] = Deri[3];
			dPhi[4][i1] = Deri[4];
			dPhi[5][i1] = Deri[5];
		}
	}


	for(i1=0; i1<=1; i1++)
	{
		delete [] NI[i1];
	}

	delete [] NI;


  	return dPhi;
}





MuHards MuHardsphere(int iz, double **DPh)
{
	MuHards MuHs;
	double  muHsz[2], muNI[6];
	double  zz, minz, maxz, zz1, zz2, zz3;
	double  coe21, coe22, coe23, coe31, coe32, coe33;
	int     minzn, maxzn, i1, i2, i3;

	zz   = iz*dz1;
    for(i1=0; i1<=1; i1++)
	{
		minz  = zz - Radii[i1];
		maxz  = zz + Radii[i1];

		minzn = int(minz/dz1 + 1.0e-8);
		maxzn = int(maxz/dz1 + 1.0e-8);
		if(minzn<rMinz1[i1])
		{
			minzn = rMinz1[i1];
		}
		if(maxzn>(Nz1+rMinz1[i1]))
		{
			maxzn = Nz1+rMinz1[i1];
		}

	    for(i3=0; i3<=5; i3++)
	    {
	    	muNI[i3] = 0.0;
		}

		for(i2=minzn; i2<maxzn-1; i2=i2+2)
		{
  	        zz1      = i2*dz1;
  	        zz2      = zz1 + dz1;
          	zz3      = zz2 + dz1;
          	coe21    = Radii2[i1] - (zz-zz1)*(zz-zz1);
          	coe22    = Radii2[i1] - (zz-zz2)*(zz-zz2);
          	coe23    = Radii2[i1] - (zz-zz3)*(zz-zz3);
          	coe31    = zz - zz1;
          	coe32    = zz - zz2;
          	coe33    = zz - zz3;

  	        muNI[0]  = muNI[0] + DPh[0][i2];
  	        muNI[0]  = muNI[0] + 4.0*DPh[0][i2+1];
  	        muNI[0]  = muNI[0] + DPh[0][i2+2];

  	        muNI[1]  = muNI[1] + DPh[1][i2];
  	        muNI[1]  = muNI[1] + 4.0*DPh[1][i2+1];
  	        muNI[1]  = muNI[1] + DPh[1][i2+2];

  	        muNI[2]  = muNI[2] + DPh[2][i2];
  	        muNI[2]  = muNI[2] + 4.0*DPh[2][i2+1];
  	        muNI[2]  = muNI[2] + DPh[2][i2+2];

  	        muNI[3]  = muNI[3] + DPh[3][i2]*coe21;
  	        muNI[3]  = muNI[3] + 4.0*DPh[3][i2+1]*coe22;
  	        muNI[3]  = muNI[3] + DPh[3][i2+2]*coe23;

  	        muNI[4]  = muNI[4] + DPh[4][i2]*coe31;
  	        muNI[4]  = muNI[4] + 4.0*DPh[4][i2+1]*coe32;
  	        muNI[4]  = muNI[4] + DPh[4][i2+2]*coe33;

  	        muNI[5]  = muNI[5] + DPh[5][i2]*coe31;
  	        muNI[5]  = muNI[5] + 4.0*DPh[5][i2+1]*coe32;
  	        muNI[5]  = muNI[5] + DPh[5][i2+2]*coe33;

        }
        if(((maxzn-minzn)%2 != 0)&&(maxzn != minzn))
        {
  	        zz2      = maxzn*dz1;
  	        zz1      = zz2 - dz1;
  	        coe22    = Radii2[i1] - (zz-zz2)*(zz-zz2);
  	        coe21    = Radii2[i1] - (zz-zz1)*(zz-zz1);
  	        coe32    = zz - zz2;
  	        coe31    = zz - zz1;


			muNI[0]  = muNI[0]/3.0 + 0.5*(DPh[0][maxzn-1]+DPh[0][maxzn]);
			muNI[1]  = muNI[1]/3.0 + 0.5*(DPh[1][maxzn-1]+DPh[1][maxzn]);
			muNI[2]  = muNI[2]/3.0 + 0.5*(DPh[2][maxzn-1]+DPh[2][maxzn]);
			muNI[3]  = muNI[3]/3.0 + 0.5*(DPh[3][maxzn-1]*coe21+DPh[3][maxzn]*coe22);
			muNI[4]  = muNI[4]/3.0 + 0.5*(DPh[4][maxzn-1]*coe31+DPh[4][maxzn]*coe32);
			muNI[5]  = muNI[5]/3.0 + 0.5*(DPh[5][maxzn-1]*coe31+DPh[5][maxzn]*coe32);


			muNI[0]  = muNI[0]*dz1/D[i1];
	        muNI[1]  = muNI[1]*dz1/2.0;
		    muNI[2]  = muNI[2]*dz1*Pi1*D[i1];
		    muNI[3]  = muNI[3]*dz1*Pi1;
		    muNI[4]  = muNI[4]*dz1/D[i1];
		    muNI[5]  = muNI[5]*dz1*Pi2;
        }
        else
        {

			muNI[0]  = muNI[0]*dz1/D[i1]/3.0;
	        muNI[1]  = muNI[1]*dz1/6.0;
		    muNI[2]  = muNI[2]*dz1*Pi1*D[i1]/3.0;
		    muNI[3]  = muNI[3]*dz1*Pi1/3.0;
		    muNI[4]  = muNI[4]*dz1/D[i1]/3.0;
		    muNI[5]  = muNI[5]*dz1*Pi2/3.0;
		}


		muHsz[i1] = 0.0;
		for(i3=0; i3<=5; i3++)
		{
			muHsz[i1] = muHsz[i1] + muNI[i3];
	    }

    }

    MuHs.pos = muHsz[0];
    MuHs.neg = muHsz[1];

    return MuHs;
}





MuHards MuHardbulk()
{
	MuHards  mu0;
	double   N0I[2], N1I[2], N2I[2], N3I[2], mu1[2];
	double   N0, N1, N2, N3, dPhi0, dPhi1, dPhi2, dPhi3;
	double   N331, N332, N333, N22, N32, PiD2;
	int      i1;

	for(i1=0; i1<=1; i1++)
	{
		N0I[i1] = rho0[i1];
		N1I[i1] = rho0[i1]*D[i1]/2.0;
		N2I[i1] = rho0[i1]*D2[i1]*Pi1;
		N3I[i1] = rho0[i1]*D2[i1]*D[i1]*Pi1/6.0;
	}

	N0    = N0I[0] + N0I[1];
	N1    = N1I[0] + N1I[1];
	N2    = N2I[0] + N2I[1];
	N3    = N3I[0] + N3I[1];

	N22   = N2*N2;
	N32   = N3*N3;
	N331  = 1.0-N3;
	N332  = N331*N331;
	N333  = N332*N331;

	dPhi0 = -log(N331);
	dPhi1 = N2/N331;
	dPhi2 = N1/N331+N22/(12.0*Pi1)*(log(N331)/N32+1.0/(N3*N332));
	dPhi3 = N0/N331+N1*N2/N332-N22*N2/(36.0*Pi1)*(2.0*log(N331)/(N32*N3)\
	        +(2.0-5.0*N3+N32)/(N32*N333));

	for(i1=0; i1<=1; i1++)
	{
		PiD2   = Pi1*D2[i1];
		mu1[i1] = dPhi0 + 0.5*D[i1]*dPhi1 + PiD2*dPhi2 + PiD2*D[i1]*dPhi3/6.0;
	}
	mu0.pos = mu1[0];
	mu0.neg = mu1[1];
    return mu0;
}


double *DeltaPhi(double *deltaPhi)
{
	double **NI;
	double N0, N1, N2, N3, NV1, NV2, Phi,Phi0;
	double N0I[2], N1I[2], N2I[2], N3I[2];
	double N331, N332, N23,N32, NV22;
	int    i1;

	NI = new double*[2];
	for(i1=0; i1<=1; i1++)
	{
		NI[i1] = new double[6];
	}

//for bulk ///////////////////////////////////////////////
	for(i1=0; i1<=1; i1++)
	{
		N0I[i1] = rho0[i1];
		N1I[i1] = rho0[i1]*D[i1]/2.0;
		N2I[i1] = rho0[i1]*D2[i1]*Pi1;
		N3I[i1] = rho0[i1]*D2[i1]*D[i1]*Pi1/6.0;
	}

	N0    = N0I[0] + N0I[1];
	N1    = N1I[0] + N1I[1];
	N2    = N2I[0] + N2I[1];
	N3    = N3I[0] + N3I[1];

	N331  = 1.0-N3;
	N332  = N331*N331;
	N23   = N2*N2*N2;
	N32   = N3*N3;

	Phi0  = -N0*log(N331) + N1*N2/N331 + (log(N331)/N32+1.0/(N332*N3))\
	        *N23/(36.0*Pi1);

///////////////////////////////////////////////////////////
	for(i1=0; i1<=Nz1; i1++)
	{

		NI = WeightD(i1, NI);

		N0  = NI[0][0]  + NI[1][0];
		N1  = NI[0][1]  + NI[1][1];
		N2  = NI[0][2]  + NI[1][2];
		N3  = NI[0][3]  + NI[1][3];
		NV1 = NI[0][4]  + NI[1][4];
		NV2 = NI[0][5]  + NI[1][5];

		if(N3>=1.0) N3=0.75;

		if(N3<1.0e-10)
		{
			Phi = 0.0;
		}
		else
		{
			N331    = 1.0-N3;
			N332    = N331*N331;
			N23     = N2*N2*N2;
			N32     = N3*N3;
			NV22    = NV2*NV2;

			Phi     = -N0*log(N331) + (N1*N2-NV1*NV2)/N331 + (log(N331)/N32+1.0/(N332*N3))\
	        *(N23-3.0*N2*NV22)/(36.0*Pi1);
		}
		deltaPhi[i1] = Phi - Phi0;
	}


	for(i1=0; i1<=1; i1++)
	{
		delete [] NI[i1];
	}

	delete [] NI;


  	return deltaPhi;

}



//using finite difference method and chasing method
double GreenBEta(double k, double eta)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green;
	double greenk, k2, bi;
	int    i1, zn2;

	A     = new double[Nz2+1];
	B     = new double[Nz2+1];
	C     = new double[Nz2+1];
	E     = new double[Nz2+1];
	L     = new double[Nz2+1];
	R     = new double[Nz2];
	Y     = new double[Nz2+1];
	green = new double[Nz2+1];


	A[0]   = 0.0;
	B[0]   = -1.0;
	C[0]   = 1.0;
	E[0]   = 0.0;
	A[Nz2] = 1.0;
	B[Nz2] = -1.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	k2     = k*k;
	bi     = 2.0 + (kappa02*eta+k2)*dz22;
	zn2    = Nz2/2;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		B[i1] = -bi;
		C[i1] = 1.0;
		E[i1] = 0.0;
	}

	E[zn2]    = -dz2/epsilS;

	//chasing method
	L[0] = B[0];
	Y[0] = E[0]/L[0];
	for(i1=1; i1<=Nz2; i1++)
	{
		R[i1-1] = C[i1-1]/L[i1-1];
		L[i1]   = B[i1] - A[i1]*R[i1-1];
		Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
	}

	green[Nz2] = Y[Nz2];
	for(i1=Nz2-1; i1>=zn2; i1--)
	{
		green[i1] = Y[i1] - R[i1]*green[i1+1];
	}

	greenk = green[zn2];

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] green;

	return greenk;
}

//using finite difference method and chasing method
double GreenEta(double k, int zn, double eta)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green;
	double greenkz, k2, bi;
	int    i1;

	A     = new double[Nz2+1];
	B     = new double[Nz2+1];
	C     = new double[Nz2+1];
	E     = new double[Nz2+1];
	L     = new double[Nz2+1];
	R     = new double[Nz2];
	Y     = new double[Nz2+1];
	green = new double[Nz2+1];




	A[0]   = 0.0;
	B[0]   = -(1.0 + 0.5*k*(epsilP/epsilS)*dz2);
	C[0]   = 1.0 - 0.5*k*(epsilP/epsilS)*dz2;
	E[0]   = 0.0;
	A[Nz2] = 1.0;
	B[Nz2] = -1.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	k2     = k*k;
	bi     = 2.0 + k2*dz22;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
	    C[i1] = 1.0;
    	E[i1] = 0.0;
    	if(rMinz2[0]>rMinz2[1])
    	{
	    	if(i1<rMinz2[1])
    		{
    			B[i1] = -bi;
    		}
    		else
     		{
    			B[i1] = -(2.0 + (kappa2[i1]*eta + k2)*dz22);
     		}
    	}
     	else
     	{
    		if(i1<rMinz2[0])
    		{
    			B[i1] = -bi;
    		}
    		else
    		{
    			B[i1] = -(2.0 + (kappa2[i1]*eta + k2)*dz22);
    		}
    	}

    }

    E[zn]  = -dz2/epsilS;

    //chasing method
    L[0] = B[0];
    Y[0] = E[0]/L[0];
    for(i1=1; i1<=Nz2; i1++)
    {
    	R[i1-1] = C[i1-1]/L[i1-1];
    	L[i1]   = B[i1] - A[i1]*R[i1-1];
    	Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
    }

    green[Nz2] = Y[Nz2];
    for(i1=Nz2-1; i1>=zn; i1--)
    {
    	green[i1] = Y[i1] - R[i1]*green[i1+1];
    }

    greenkz = green[zn];



	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] green;

	return greenkz;
}



double GreenFree(double k, int zn)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green;
	double greenkz, k2, bi;
	int    i1;

	A     = new double[Nz2+1];
	B     = new double[Nz2+1];
	C     = new double[Nz2+1];
	E     = new double[Nz2+1];
	L     = new double[Nz2+1];
	R     = new double[Nz2];
	Y     = new double[Nz2+1];
	green = new double[Nz2+1];


	A[0]   = 0.0;
	B[0]   = -(1.0 + 0.5*k*dz2);
	C[0]   = 1.0 - 0.5*k*dz2;
	E[0]   = 0.0;
	A[Nz2] = -(1.0 + 0.5*k*dz2);
	B[Nz2] = 1.0 - 0.5*k*dz2;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	k2     = k*k;
	bi     = 2.0 + k2*dz22;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
	    C[i1] = 1.0;
    	E[i1] = 0.0;
    	B[i1] = -bi;

    }

    E[zn]  = -dz2/epsilS;

    //chasing method
    L[0] = B[0];
    Y[0] = E[0]/L[0];
    for(i1=1; i1<=Nz2; i1++)
    {
    	R[i1-1] = C[i1-1]/L[i1-1];
    	L[i1]   = B[i1] - A[i1]*R[i1-1];
    	Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
    }

    green[Nz2] = Y[Nz2];
    for(i1=Nz2-1; i1>=zn; i1--)
    {
    	green[i1] = Y[i1] - R[i1]*green[i1+1];
    }

    greenkz = green[zn];



	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] green;

	return greenkz;
}


///////////////////////////////////////////fast-method///////
//using finite difference method and chasing method


double **GreenBEtaFast(double **GreenBEtaK)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green, *AL;
	double greenk, k2, bi, eta,k,const6,kk;
	int    i1, i2, i3, zn2,zn22;

	A     = new double[Nz2+1] ();
	B     = new double[Nz2+1] ();
	C     = new double[Nz2+1] ();
	E     = new double[Nz2+1] ();
	L     = new double[Nz2+1] ();
	R     = new double[Nz2] ();
	Y     = new double[Nz2+1] ();
	AL    = new double[Nz2+1] ();
	green = new double[Nz2+1] ();


	A[0]   = 0.0;
	E[0]   = 0.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	//B[0]   = -1.0;
	//C[0]   = 1.0;
	//A[Nz2] = 1.0;
	//B[Nz2] = -1.0;

	zn2    = Nz2/2;

	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		C[i1] = 1.0;
		E[i1] = 0.0;
	}

	E[zn2] = const2;

	for(i3=1; i3<=Neta; i3++)
	{
		eta = i3*deta;

		for(i2=1; i2<=Nk; i2++)
		{
			k      = i2*dk;
			k2     = k*k;
			kk     = sqrt(kappa02*eta+k2);

        	B[0]   = -(1.0 + kk*const4);
        	C[0]   = 1.0 - kk*const4;
    	 	A[Nz2] = 1.0 - kk*const4;
        	B[Nz2] = -(1.0 + kk*const4);

	        bi     = 2.0 + (kappa02*eta+k2)*dz22;

    	    for(i1=1; i1<=(Nz2-1); i1++)
            {
        		B[i1] = -bi;
        	}

        	//chasing method
        	L[0]   = B[0];
        	Y[0]   = E[0]/L[0];
        	for(i1=1; i1<=Nz2-1; i1++)
        	{
	        	R[i1-1] = C[i1-1]/L[i1-1];
        		L[i1]   = B[i1] - R[i1-1];
        	//	Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
        	    AL[i1]  = -A[i1]/L[i1];
        	}

			R[Nz2-1] = C[Nz2-1]/L[Nz2-1];
        	L[Nz2]   = B[Nz2] - A[Nz2]*R[Nz2-1];
        	AL[Nz2]  = -A[Nz2]/L[Nz2];

            if(AL[Nz2]<0.0)
    		{
    			cerr<<"WRONG !!! AL < 0"<<endl;
    			exit(0);
	    	}

        	//chasing method
            Y[zn2] = E[zn2]/L[zn2];
            zn22 = 0;
         	for(i1=zn2+1; i1<=Nz2; i1++)
        	{
    	     	const6 = const5/AL[Nz2];
				if(Y[i1-1]<const6)
    	     	{
	         		Y[i1]   = 0.0;
    	     		zn22    = i1;
    	     		i1      = Nz2 + 2;
    			}
	    		else
	    		{
	    			Y[i1]   = AL[i1]*Y[i1-1];
	    		}
        	}

		    if(zn22!=0)
    	    {
    	    	green[zn22] = Y[zn22];
            	for(i1=zn22-1; i1>=zn2; i1--)
             	{
	            	green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

             	GreenBEtaK[i3][i2] = green[zn2];
	    	}
    		else
     		{
	    		green[Nz2] = Y[Nz2];
            	for(i1=Nz2-1; i1>=zn2; i1--)
            	{
	            	green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

             	GreenBEtaK[i3][i2] = green[zn2];
     		}

		}
	}

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] AL;
	delete [] green;

	return GreenBEtaK;
}

//using finite difference method and chasing method
double ***GreenEtaFast(double ***GreenEtaKZ)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green, *AL;
	double greenkz, k2, bi,eta, k, const6,kk;
	int    i1,i2,i3,i4, zn2, zn22;

	A     = new double[Nz2+1] ();
	B     = new double[Nz2+1] ();
	C     = new double[Nz2+1] ();
	E     = new double[Nz2+1] ();
	L     = new double[Nz2+1] ();
	R     = new double[Nz2] ();
	Y     = new double[Nz2+1] ();
	AL    = new double[Nz2+1] ();
	green = new double[Nz2+1] ();


	A[0]   = 0.0;
	E[0]   = 0.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;
	//A[Nz2] = 1.0;
	//B[Nz2] = -1.0;

	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
	    C[i1] = 1.0;
    	E[i1] = 0.0;
    }

	for(i4=1; i4<=Neta; i4++)
	{
		eta = i4*deta;
		for(i3=1; i3<=Nk; i3++)
		{
			k      = i3*dk;
    		k2     = k*k;
        	kk     = sqrt(kappa02*eta+k2);

			B[0]   = -(1.0 + k*const1);
        	C[0]   = 1.0 - k*const1;
    	 	A[Nz2] = 1.0 - k*const1;
        	B[Nz2] = -(1.0 + k*const1);

        	bi     = 2.0 + k2*dz22;
        	for(i1=1; i1<=(Nz2-1); i1++)
        	{

            	if(rMinz2[0]>rMinz2[1])
            	{
        	    	if(i1<rMinz2[1])
            		{
            			B[i1] = -bi;
            		}
            		else
            		{
            			B[i1] = -(2.0 + (kappa2[i1]*eta + k2)*dz22);
            		}
             	}
            	else
            	{
            		if(i1<rMinz2[0])
            		{
            			B[i1] = -bi;
             		}
            		else
            		{
            			B[i1] = -(2.0 + (kappa2[i1]*eta + k2)*dz22);
            		}
            	}

            }

            L[0] = B[0];
            Y[0] = E[0]/L[0];
            for(i1=1; i1<=Nz2-1; i1++)
            {
            	R[i1-1] = C[i1-1]/L[i1-1];
            	L[i1]   = B[i1] - R[i1-1];
            //	Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
                AL[i1]  = -A[i1]/L[i1];
            }

            R[Nz2-1] = C[Nz2-1]/L[Nz2-1];
            L[Nz2]   = B[Nz2] - A[Nz2]*R[Nz2-1];
            AL[Nz2]  = -A[Nz2]/L[Nz2];


            if(AL[Nz2]<0.0)
    		{
    			cerr<<"WRONG !!! AL < 0"<<endl;
    			exit(0);
	    	}

			for(i2=1; i2<=Nz1; i2++)
			{
                zn2     = int(i2*const3 + 0.5);
        	   	E[zn2]  = const2;
        	   	Y[zn2]  = E[zn2]/L[zn2];
	            zn22    = 0;
                //chasing method
            	for(i1=zn2 + 1; i1<=Nz2; i1++)
             	{
	         	//R[i1-1] = C[i1-1]/L[i1-1];
	        	//L[i1]   = B[i1] - A[i1]*R[i1-1];
	        	    const6 = const5/AL[Nz2];
    	        	if(Y[i1-1]<const6)
	            	{
	             		Y[i1]   = 0.0;
	             		zn22    = i1;
	            		i1      = Nz2 + 2;
	    			}
	    			else
	     			{
	    				Y[i1]   = AL[i1]*Y[i1-1];
	    			}

             	}
    	        if(zn22!=0)
	            {
	            	green[zn22] = Y[zn22];
            	    for(i1=zn22-1; i1>=zn2; i1--)
                	{
                		green[i1] = Y[i1] - R[i1]*green[i1+1];
                	}

                 	GreenEtaKZ[i4][i3][i2] = green[zn2];
                 	E[zn2] = 0.0;
    			}
	    		else
	    		{
	    			green[Nz2] = Y[Nz2];
            	    for(i1=Nz2-1; i1>=zn2; i1--)
                	{
                		green[i1] = Y[i1] - R[i1]*green[i1+1];
                	}

                	GreenEtaKZ[i4][i3][i2] = green[zn2];
                	E[zn2] = 0.0;
	    		}


			}

		}
	}


	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] AL;
	delete [] green;

	return GreenEtaKZ;
}



double **GreenFreeFast(double **GreenFreeKZ)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green, *AL;
	double greenkz, k2, bi, k, ALN;
	int    i1, i2, i3, kn, zn2,zn22;

	A     = new double[Nz2+1] ();
	B     = new double[Nz2+1] ();
	C     = new double[Nz2+1] ();
	E     = new double[Nz2+1] ();
	L     = new double[Nz2+1] ();
	R     = new double[Nz2] ();
	Y     = new double[Nz2+1] ();
	AL    = new double[Nz2+1] ();
	green = new double[Nz2+1] ();


	A[0]   = 0.0;
	E[0]   = 0.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;
	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
	    C[i1] = 1.0;
    	E[i1] = 0.0;
    }


	for(i3=1; i3<=Nk; i3++)
	{
		k      = i3*dk;
		B[0]   = -(1.0 + k*const4);
    	C[0]   = 1.0 - k*const4;
    	A[Nz2] = 1.0 - k*const4;
    	B[Nz2] = -(1.0 + k*const4);
    	k2     = k*k;
    	bi     = 2.0 + k2*dz22;
		for(i1=1; i1<=(Nz2-1); i1++)
    	{
        	B[i1] = -bi;
        }

		//chasing method
        L[0] = B[0];
        Y[0] = E[0]/L[0];
        ALN  = 0.0;
        for(i1=1; i1<=Nz2-1; i1++)
        {
        	R[i1-1] = C[i1-1]/L[i1-1];
        	L[i1]   = B[i1] - R[i1-1];
        	//Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
        	AL[i1]  = -A[i1]/L[i1];
        }

        R[Nz2-1] = C[Nz2-1]/L[Nz2-1];
        L[Nz2]   = B[Nz2] - A[Nz2]*R[Nz2-1];
        AL[Nz2]  = -A[Nz2]/L[Nz2];

        ALN = AL[Nz2];

        if(ALN<0.0)
		{
			cerr<<"WRONG !!! AL < 0"<<endl;
			exit(0);
		}
		else if(ALN<1.0)
		{
			ALN = 1.0;
		}

        for(i2=1; i2<=Nz1; i2++)
        {
        	zn2    = int(i2*const3 + 0.5);
        	E[zn2] = const2;
        	Y[zn2]  = E[zn2]/L[zn2];
        	zn22 = 0;
        	//chasing method
            for(i1=zn2 + 1; i1<=Nz2; i1++)
            {
 	         	if((Y[i1-1]*ALN)<const5)
	        	{
    	     		Y[i1]   = 0.0;
	         		zn22    = i1;
    	     		i1      = Nz2 + 2;
	    		}
    			else
	    		{
	    			Y[i1]   = AL[i1]*Y[i1-1];
	    		}
            }

	        if(zn22!=0)
	        {
	        	green[zn22] = Y[zn22];
        	    for(i1=zn22-1; i1>=zn2; i1--)
            	{
            		green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

            	GreenFreeKZ[i3][i2] = green[zn2];
            	E[zn2] = 0.0;
			}
			else
			{
				green[Nz2] = Y[Nz2];
        	    for(i1=Nz2-1; i1>=zn2; i1--)
            	{
            		green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

            	GreenFreeKZ[i3][i2] = green[zn2];
            	E[zn2] = 0.0;
			}

            //GreenFreeKZ[i3][i2] = green[zn2];
		}

	}

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] AL;
	delete [] green;

	return GreenFreeKZ;
}


//using finite difference method and chasing method
double *GreenBFast(double *GreenBK)
{
	double *A, *B, *C, *E, *L, *R, *Y, *green, *AL;
	double greenk, k2, bi, k,const6,kk;
	int    i1, i2, zn2, zn22;

	A     = new double[Nz2+1] ();
	B     = new double[Nz2+1] ();
	C     = new double[Nz2+1] ();
	E     = new double[Nz2+1] ();
	L     = new double[Nz2+1] ();
	R     = new double[Nz2] ();
	Y     = new double[Nz2+1] ();
	AL    = new double[Nz2+1] ();
	green = new double[Nz2+1] ();

	A[0]   = 0.0;
	E[0]   = 0.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;
	//B[0]   = -1.0;
	//C[0]   = 1.0;
	//A[Nz2] = 1.0;
	//B[Nz2] = -1.0;
	zn2    = Nz2/2;

	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		C[i1] = 1.0;
		E[i1] = 0.0;
	}
	E[zn2] = const2;

    for(i2=1; i2<=Nk; i2++)
	{
		k    = i2*dk;
		k2   = k*k;
    	kk   = sqrt(kappa02+k2);

    	B[0]   = -(1.0 + kk*const4);
    	C[0]   = 1.0 - kk*const4;
	 	A[Nz2] = 1.0 - kk*const4;
    	B[Nz2] = -(1.0 + kk*const4);

	    bi   = 2.0 + (kappa02+k2)*dz22;
	    for(i1=1; i1<=(Nz2-1); i1++)
    	{
    		B[i1] = -bi;
    	}


    	L[0] = B[0];
    	Y[0] = E[0]/L[0];
     	for(i1=1; i1<=Nz2-1; i1++)
    	{
    		R[i1-1] = C[i1-1]/L[i1-1];
    		L[i1]   = B[i1] - R[i1-1];
    	//	Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
    	    AL[i1]  = -A[i1]/L[i1];
    	}

   		R[Nz2-1] = C[Nz2-1]/L[Nz2-1];
    	L[Nz2]   = B[Nz2] - A[Nz2]*R[Nz2-1];
    	AL[Nz2]  = -A[Nz2]/L[Nz2];

        if(AL[Nz2]<0.0)
		{
			cerr<<"WRONG !!! AL < 0"<<endl;
			exit(0);
		}

		//chasing method
		Y[zn2] = E[zn2]/L[zn2];

        zn22 = 0;
     	for(i1=zn2 + 1; i1<=Nz2; i1++)
     	{
	     	//R[i1-1] = C[i1-1]/L[i1-1];
	     	//L[i1]   = B[i1] - A[i1]*R[i1-1];
	     	const6 = const5/AL[Nz2];
	     	if(Y[i1-1]<const6)
	     	{
	     		Y[i1]   = 0.0;
	     		zn22    = i1;
	     		i1      = Nz2 + 2;
			}
			else
			{
				Y[i1]   = AL[i1]*Y[i1-1];
			}


     	}

	    if(zn22!=0)
	    {
	    	green[zn22] = Y[zn22];
        	for(i1=zn22-1; i1>=zn2; i1--)
        	{
	        	green[i1] = Y[i1] - R[i1]*green[i1+1];
        	}

         	GreenBK[i2] = green[zn2];
		}
		else
		{
			green[Nz2] = Y[Nz2];
        	for(i1=Nz2-1; i1>=zn2; i1--)
        	{
	        	green[i1] = Y[i1] - R[i1]*green[i1+1];
        	}

         	GreenBK[i2] = green[zn2];
		}

	}


	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] AL;
	delete [] green;

	return GreenBK;
}



//using finite difference method and chasing method
double **GreenFast(double **GreenKZn)
{
	double *A, *B, *C, *E, *L, *R, *Y, *AL, *green;
	double greenkz, k2, bi,k,const6,kk;
	int    i1, i2, i3, kn, zn2, zn22;

	A     = new double[Nz2+1] ();
	B     = new double[Nz2+1] ();
	C     = new double[Nz2+1] ();
	E     = new double[Nz2+1] ();
	L     = new double[Nz2+1] ();
	R     = new double[Nz2] ();
	Y     = new double[Nz2+1] ();
	AL    = new double[Nz2+1] ();
	green = new double[Nz2+1] ();

	A[0]   = 0.0;
	E[0]   = 0.0;
	C[Nz2] = 0.0;
	E[Nz2] = 0.0;

	//A[Nz2] = 1.0;
	//B[Nz2] = -1.0;

	for(i1=1; i1<=(Nz2-1); i1++)
	{
		A[i1] = 1.0;
		C[i1] = 1.0;
		E[i1] = 0.0;
	}


	for(i3=1; i3<=Nk; i3++)
	{
	    k      = i3*dk;
    	k2     = k*k;
    	kk     = sqrt(kappa02+k2);

    	B[0]   = -(1.0 + k*const1);
    	C[0]   = 1.0 - k*const1;
	 	B[Nz2] = -(1.0 + k*const1);
    	A[Nz2] = 1.0 - k*const1;

     	bi     = 2.0 + k2*dz22;
    	for(i1=1; i1<=(Nz2-1); i1++)
    	{
    		if(rMinz2[0]>rMinz2[1])
    		{
	    		if(i1<rMinz2[1])
	    		{
	    			B[i1] = -bi;
	    		}
	    		else
	    		{
	    			B[i1] = -(2.0 + (kappa2[i1] + k2)*dz22);
	    		}
	    	}
	    	else
	    	{
	    		if(i1<rMinz2[0])
	    		{
	    			B[i1] = -bi;
	    		}
	    		else
	    		{
	    			B[i1] = -(2.0 + (kappa2[i1] + k2)*dz22);
		    	}
	    	}

     	}


		L[0] = B[0];
        Y[0] = E[0]/L[0];
        for(i1=1; i1<=Nz2-1; i1++)
        {
	        R[i1-1] = C[i1-1]/L[i1-1];
	        L[i1]   = B[i1] - R[i1-1];   //= B[i1] - A[i1]*R[i1-1];
	        //Y[i1]   = (E[i1]-A[i1]*Y[i1-1])/L[i1];
	        AL[i1]  = -A[i1]/L[i1];
        }

	    R[Nz2-1] = C[Nz2-1]/L[Nz2-1];
	    L[Nz2]   = B[Nz2] - A[Nz2]*R[Nz2-1];   //= B[i1] - A[i1]*R[i1-1];
	    AL[Nz2]  = -A[Nz2]/L[Nz2];

    	for(i2=1; i2<=Nz1; i2++)
    	{

	        zn2     = int(i2*const3 + 0.5);
        	E[zn2]  = const2;
        	Y[zn2]  = E[zn2]/L[zn2];
        	//chasing method
        	zn22 = 0;
        	for(i1=zn2 + 1; i1<=Nz2; i1++)
        	{
	         	//R[i1-1] = C[i1-1]/L[i1-1];
	        	//L[i1]   = B[i1] - A[i1]*R[i1-1];
	        	const6 = const5/AL[Nz2];
	        	if(Y[i1-1]<const6)
	        	{
	        		Y[i1]   = 0.0;
	        		zn22    = i1;
	        		i1      = Nz2 + 2;
				}
				else
				{
					Y[i1]   = AL[i1]*Y[i1-1];
				}

         	}
	        if(zn22!=0)
	        {
	        	green[zn22] = Y[zn22];
        	    for(i1=zn22-1; i1>=zn2; i1--)
            	{
            		green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

            	GreenKZn[i3][i2] = green[zn2];
            	E[zn2] = 0.0;
			}
			else
			{
				green[Nz2] = Y[Nz2];
        	    for(i1=Nz2-1; i1>=zn2; i1--)
            	{
            		green[i1] = Y[i1] - R[i1]*green[i1+1];
            	}

            	GreenKZn[i3][i2] = green[zn2];
            	E[zn2] = 0.0;
			}

        }

    }
	///////////////////////////////////////////////////////////////////////////////////////////////

	delete [] A;
	delete [] B;
	delete [] C;
	delete [] E;
	delete [] L;
	delete [] R;
	delete [] Y;
	delete [] AL;
	delete [] green;

	return GreenKZn;
}
