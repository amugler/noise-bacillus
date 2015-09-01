/*  This c++ code generates time series data for a Gillespie simulation for the stochastic model based off of the
    deterministic equations for the SynExSlow strain with parameters corresponding to high molecule numbers.
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
		
	// dimensionless parameters
	float a_k = 0.;
	
	
	int moleculesK = 1;
	int moleculesM = 1;
	int moleculesS = 1;
	
	double alpha_k = a_k*0.04;
	double alpha_m = 0.075;
	double alpha_s = 0.5;
	double beta_k = 7.5;
	double beta_m = 2.5;
	double beta_s = 0.5;
	int k_k = 5000;
	int k_m = 2500;
	int k_s = 500;
	long double delta = .000002;
	long double lambda = .0001;
	int gamma_k = 25000;
	int gamma_s = 20;
	int h = 2;
	int p = 2;
	
	
	
	
	// time parameters
	double t = 0;
	double tau = 0;
	float time_hours = 10000;
	float tmax = time_hours*3600;
	
	// sets random seed for random numbers
  	srand((unsigned)time(NULL));
  	
  	// declaring files to which to write
	ofstream stocfileTKMS;
  	
  	// declares array for storing values
	int storeTKMSmax = 1000000000;
	double storeT[storeTKMSmax];
	int storemoleculesK[storeTKMSmax];
	int storemoleculesM[storeTKMSmax];
	int storemoleculesS[storeTKMSmax];
	
	// declares counting variables
	int i = 0;
	int n = 0;	// number of times copied to vector
	int m = 0;	// number of reactions
	
	// declares other variables
	int samplerate = 100;


		while (t<=tmax)
		{  	
		
			if (n <= storeTKMSmax and m%samplerate == 0)
			{
				storeT[n] = t;
				storemoleculesK[n] = moleculesK;
				storemoleculesM[n] = moleculesM;
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
			double a1 = alpha_k + beta_k*pow(moleculesK,h)/(pow(k_k,h)+pow(moleculesK,h));
			double a2 = (delta*moleculesM*moleculesK)/(1+moleculesK/gamma_k + moleculesS/gamma_s) + lambda*moleculesK;
			double a3 = alpha_m + beta_m*pow(moleculesK,p)/(pow(k_m,p)+pow(moleculesK,p));
			double a4 = lambda*moleculesM;
			double a5 = alpha_s + beta_s*pow(moleculesK,h)/(pow(k_s,h)+pow(moleculesK,h));
			double a6 = (delta*moleculesM*moleculesS)/(1+moleculesK/gamma_k + moleculesS/gamma_s) + lambda*moleculesS;
	
			double a0 = a1 + a2 + a3 + a4 + a5 + a6;
			
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
				else if ((a1+a2+a3)/a0 <= r2 and r2 < (a1+a2+a3+a4)/a0) 
				{
					moleculesM -= 1;
				}
				else if ((a1+a2+a3+a4)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5)/a0) 
				{
					moleculesS += 1;
				}
				else if ((a1+a2+a3+a4+a5)/a0 <= r2 and r2 < 1)  
				{
					moleculesS -= 1;
				}
				else
				{
					cout << "ERROR" << "\n";
				}
			}
			
			// cout << t << ",     " << n << ",     " << moleculesK << ",     " << moleculesM << ",     " << moleculesS << "\n";
		}
		
		// Converts a_k to string
		stringstream ss;
		ss << fixed << setprecision(4) << a_k;
		string string_ak(ss.str());
		string_ak.erase(1,1);
		
		string filename_ak = "gssslow_high_TKMS_" + string_ak + ".txt";
		
		
		
		//writing the t, K, and M values to file
		stocfileTKMS.open(filename_ak.c_str(), ios::trunc);
		stocfileTKMS.close();
		
		stocfileTKMS.open (filename_ak.c_str(), ios::app);
		
		while (i < n)
		{
			if (storeT[i] >= 0)
			{
				stocfileTKMS << fixed << setprecision(15) << storeT[i] << ", " << storemoleculesK[i] << ", " << storemoleculesM[i] << ", " << storemoleculesS[i] << ",\n";
			}
			i += 1;
		}
		stocfileTKMS.close();		
  	program_runtime = clock() - program_runtime;
  	cout <<"The program took " << ((float)program_runtime)/CLOCKS_PER_SEC << " seconds to run." << "\n";
}
