/*  This c++ code generates time series data for a Gillespie simulation for the stochastic model based off of the
    deterministic equations of the native strain with parameters corresponding to low molecule numbers.
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
		
		
		
		To test: .2 x akd1, .9 x akd1, sqrt(akd1 x akd2), 1.1 x akd2, 5 x akd2
	*/
	
	// initial dimensionless values for K and S
	int K = .0127;
	int S = 100;
	
	// dimensionless parameters
	float a_k = 0.00;
	float stepsize_ak = .001;
	
	float a_s = 0;
	float b_k = 0.3;
	float b_s = 3;
	float k_0 = 0.2;
	double k_1 = 1;
	k_1 = k_1 / 30;
	float Delta_k = 0.1;
	float Delta_s = 0.1;
	float h = 2;
	float p = 5;

	
	// physical parameters
	int Gamma_k = 100;
	int Gamma_s = 1;
	
	// time parameters
	double t = 0;
	double tau = 0;
	float delta = .001;
	float time_hours = 10000;
	float tmax = time_hours*3600;
	
	// relating physical quantities to diensionless quantities, 
	int moleculesK = K*Gamma_k;
	int moleculesS = S*Gamma_s;
	float alpha_k = a_k*Gamma_k*delta;
	float alpha_s = a_s*Gamma_s*delta;
	float beta_k = b_k*Gamma_k*delta;
	float beta_s = b_s*Gamma_s*delta;
	float k_k = k_0*Gamma_k;
	float k_s = k_1*Gamma_k;
	float lambda_k = Delta_k*delta;
	float lambda_s = Delta_s*delta;
	
	// sets random seed for random numbers
  	srand((unsigned)time(NULL));
  	
  	// declaring files to which to write
	ofstream stocfileTKS;
  	
  	// declares array for storing values
	int storeTKSmax = 1000000000;
	double storeT[storeTKSmax];
	int storemoleculesK[storeTKSmax];
	int storemoleculesS[storeTKSmax];
	
	// declares counting variables
	int i = 0;
	int j = 0;
	int n = 0;	// number of times copied to vector
	int m = 0;	// number of reactions
	
	// declares other variables
	int samplerate = 75;

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
		
		string filename_ak = "gns_low_TKS_" + string_ak + ".txt";
		
		
		
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
		
		t = 0;
		tau = 0;
		
		moleculesK = K*Gamma_k;
		moleculesS = S*Gamma_s;
		alpha_k = a_k*Gamma_k*delta;
		alpha_s = a_s*Gamma_s*delta;
		beta_k = b_k*Gamma_k*delta;
		beta_s = b_s*Gamma_s*delta;
		k_k = k_0*Gamma_k;
		k_s = k_1*Gamma_k;
		lambda_k = Delta_k*delta;
		lambda_s = Delta_s*delta;
		
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
