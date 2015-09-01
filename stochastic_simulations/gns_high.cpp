/*  This c++ code generates time series data for a Gillespie simulation for the stochastic model based off of the
    deterministic equations with parameters corresponding to high molecule numbers.
*/

#include <iostream>
#include <string>		/*sring */
#include <sstream>		/*stringstream */
#include <string.h>
#include <iomanip>		/* setprecision */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h> 		/* time */
#include <fstream>		// read and write to files
#include <math.h>       /* log, pow */
using namespace std;

int main()
{
	clock_t program_runtime;
	program_runtime = clock();
	
	/* deterministic bounderies for standard given parameters:
	
							Native							SynEx 
							ak								ak
		Simple, low		0.0007, 0.0033					0.0327, 0.7565	
		
							ak								ak
		Simple, high	0.0007, 0.0033					0.0327, 0.7565
		
							k1								k1
		full, high		0.000419, 0.00219				0.00040, 0.00975
		
	*/
	
	// dimensionless parameters
	float a_k = 0.0000;
	float stepsize_ak = .001;
	
	int moleculesK = 350;
	int moleculesS = 7000;
	int Gamma_k = 25000;
	int Gamma_s = 20;
	float delta = .001;
	float h = 2;
	float p = 5;
	float alpha_k = a_k*Gamma_k*delta;
	float alpha_s = .0004;
	float beta_k = 7.5;
	float beta_s = .06;
	float k_k = 5000;
	float k_s = 833;
	float lambda_k = .0001;
	float lambda_s = .0001;
	
	
	// time parameters
	double t = 0;
	double tau = 0;
	float time_hours = 10000;
	float tmax = time_hours*3600;
	
	// sets random seed for random numbers
  	srand((unsigned)time(NULL));
  	
  	// declaring files to which to write
	ofstream stocfileTKS;
  	
  	// declares array for storing values
	int storeTKSmax = 100000000;
	double storeT[storeTKSmax];
	int storemoleculesK[storeTKSmax];
	int storemoleculesS[storeTKSmax];
	
	// declares counting variables
	int i = 0;
	int j = 0;
	int n = 0;	// number of times copied to vector
	int m = 0;	// number of reactions
	
	// declares other variables
	int samplerate = 50;

	// length of loop for doing multiple values of a_k in one run
  	while (j<1)
  	{
		while (t<=tmax)
		{  	
		
			if ( m%samplerate == 0)
			{
				storeT[n] = t;
				storemoleculesK[n] = moleculesK;
				storemoleculesS[n] = moleculesS;
				n += 1;
			}
			m += 1;
			
			// declaring random variables
			redorand:
			double r1 = ((double) rand() / (RAND_MAX));
			double r2 = ((double) rand() / (RAND_MAX));
			{
				if ((r1 == 0) or (r2 == 0) or (r1 == 1) or (r2 == 1)) {goto redorand;}
			}
			
			// propensities
			double a1 = alpha_k + (beta_k*pow(moleculesK,h))/(pow(k_k,h)+pow(moleculesK,h));
			double a2 = (delta/(1+moleculesK/Gamma_k+moleculesS/Gamma_s)+lambda_k)*moleculesK;
			double a3 = alpha_s + (beta_s/(1+pow(moleculesK/k_s,p)));
			double a4 = (delta/(1+moleculesK/Gamma_k+moleculesS/Gamma_s)+lambda_s)*moleculesS;
	
			double a0 = a1 + a2 + a3 + a4;
			
			// determine time of next event
			tau = ((1)/(a0))*(log(1/r1));
	
			// update time
			t = t + tau;
			
			//determines the outcome of the event
			{
				if (0 <= r2 and r2 < a1/a0) 
				{
					moleculesK += 1;
				}
				else if (a1/a0 <= r2 and r2 < (a1+a2)/a0) 
				{
					moleculesK -= 1;
				}
				else if ((a1+a2)/a0 <= r2 and r2 < (a1+a2+a3)/a0) 
				{
					moleculesS += 1;
				}
				else if ((a1+a2+a3)/a0 <= r2 and r2 < 1)  {
					moleculesS -= 1;
				}
				else
				{
					cout << "ERROR" << "\n";
				}
			}
			
			//cout << t << ",     " << n << ",     " << moleculesK << ",     " << moleculesS << ",     " << a_k << "\n";
		}
		
		// Converts a_k to string
		stringstream ss;
		ss << fixed << setprecision(5) << a_k;
		string string_ak(ss.str());
		string_ak.erase(0,2);
		
		string filename_ak = "gns_TKS_high_" + string_ak + ".txt";
		
		
		
		//writing the t, K, and S values to file
		stocfileTKS.open(filename_ak.c_str(), ios::trunc);
		stocfileTKS.close();
		
		stocfileTKS.open (filename_ak.c_str(), ios::app);
		
		while (i < n)
		{
			if (storeT[i] >= 0)
			{
				stocfileTKS << fixed << setprecision(15) << storeT[i] << ", " << storemoleculesK[i] << ", " << storemoleculesS[i] << ",\n";
			}
			i += 1;
		}
		stocfileTKS.close();
		
		
		a_k = a_k + stepsize_ak;
		j += 1;
		
		
		moleculesK = 350;
		moleculesS = 7000;
		Gamma_k = 25000;
		Gamma_s = 20;
		delta = .001;
		h = 2;
		p = 5;
		alpha_k = a_k*Gamma_k*delta;
		alpha_s = .0004;
		beta_k = 7.5;
		beta_s = .06;
		k_k = 5000;
		k_s = 833;
		lambda_k = .0001;
		lambda_s = .0001;
			
		t = 0;
		tau = 0;
		
	  	srand((unsigned)time(NULL));
	
		memset(storeT, 0, storeTKSmax);
		memset(storemoleculesK, 0, storeTKSmax);
		memset(storemoleculesS, 0, storeTKSmax);
		
		i = 0;
		n = 0;	// number of times copied to vector
		m = 0;	// number of reactions
		
  	}
  	
  	program_runtime = clock() - program_runtime;
  	cout <<"The program took " << ((float)program_runtime)/CLOCKS_PER_SEC << " seconds to run." << "\n";
}
