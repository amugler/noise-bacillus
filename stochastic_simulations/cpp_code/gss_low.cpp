/*  This c++ code generates time series data for a Gillespie simulation for the stochastic model based off of the
    deterministic equations for the SynEx strain with parameters corresponding to low molecule numbers.
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
	int M = 0.01;
	
	// dimensionless parameters
	float a_k = 0.0000;
	float stepsize_ak = .001;
	
	float a_m = 0.3;
	float b_k = 15;
	float b_m = 10;
	float rho = 0.5;
	float mu = 1;
	float h = 2;
	float p = 2;

	
	// physical parameters
	int k_k = 10;
	int k_m = k_k*rho;
	
	// time parameters
	double t = 0;
	double tau = 0;
	float lambda = .0001;
	float time_hours = 10000;
	float tmax = time_hours*3600;
	
	// relating physical quantities to diensionless quantities, 
	int moleculesK = K*k_k;
	int moleculesM = M*k_m;
	float alpha_k = a_k*k_k*lambda;
	float alpha_m = a_m*k_m*lambda;
	float beta_k = b_k*k_k*lambda;
	float beta_m = b_m*k_m*lambda;
	float delta = mu*lambda/k_m;
	
	// sets random seed for random numbers
  	srand((unsigned)time(NULL));
  	
  	// declaring files to which to write
	ofstream stocfileTKM;
  	
  	// declares array for storing values
	int storeTKMmax = 10000000;
	double storeT[storeTKMmax];
	int storemoleculesK[storeTKMmax];
	int storemoleculesM[storeTKMmax];
	
	// declares counting variables
	int i = 0;
	int j = 0;
	int n = 0;	// number of times copied to vector
	int m = 0;	// number of reactions
	
	// declares other variables
	int samplerate = 100;

  	while (j<1)
  	{
		while (t<=tmax)
		{  	
		
			if (n <= storeTKMmax and m%samplerate == 0)
			{
				storeT[n] = t;
				storemoleculesK[n] = moleculesK;
				storemoleculesM[n] = moleculesM;
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
			double a1 = alpha_k + beta_k*pow(moleculesK,h)/(pow(k_k,h)+pow(moleculesK,h));
			double a2 = (delta*moleculesM + lambda)*moleculesK;
			double a3 = alpha_m + beta_m*pow(moleculesK,p)/(pow(k_m,p)+pow(moleculesK,p));
			double a4 = lambda*moleculesM;
	
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
					moleculesM += 1;
				}
				else if ((a1+a2+a3)/a0 <= r2 and r2 < 1)  {
					moleculesM -= 1;
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
		ss << fixed << setprecision(4) << a_k;
		string string_ak(ss.str());
		string_ak.erase(1,1);
		
		string filename_ak = "gss_low_TKM_" + string_ak + ".txt";
		
		
		
		//writing the t, K, and M values to file
		stocfileTKM.open(filename_ak.c_str(), ios::trunc);
		stocfileTKM.close();
		
		stocfileTKM.open (filename_ak.c_str(), ios::app);
		
		while (i < n)
		{
			if (storeT[i] >= 0)
			{
				stocfileTKM << fixed << setprecision(15) << storeT[i] << ", " << storemoleculesK[i] << ", " << storemoleculesM[i] << ",\n";
			}
			i += 1;
		}
		stocfileTKM.close();
		
		cout << a_k << "\n";
		a_k = a_k + stepsize_ak;
		j += 1;
		
		t = 0;
		tau = 0;
		
		moleculesK = K*k_k;
		moleculesM = M*k_m;
		alpha_k = a_k*k_k*lambda;
		alpha_m = a_m*k_m*lambda;
		beta_k = b_k*k_k*lambda;
		beta_m = b_m*k_m*lambda;
		delta = mu*lambda/k_m;
		
	  	srand((unsigned)time(NULL));
	
		memset(storeT, 0, storeTKMmax);
		memset(storemoleculesK, 0, storeTKMmax);
		memset(storemoleculesM, 0, storeTKMmax);
		
		i = 0;
		n = 0;	// number of times copied to vector
		m = 0;	// number of reactions
		
  	}
  	program_runtime = clock() - program_runtime;
  	cout <<"The program took " << ((float)program_runtime)/CLOCKS_PER_SEC << " seconds to run." << "\n";
}
