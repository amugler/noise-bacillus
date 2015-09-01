/*  This c++ code generates time series data for a Gillespie simulation for the stochastic full model of the SynEx
    strain without adiabatic elimination of mRNA dynamics with parameters corresponding to high molecule numbers.
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
	
	
	// rates and parameters
	float k1 = 0.00008;
	float stepsize_k1 = .001;
	
	float k2 = .19;
	float k3 = .0019;
	float k4 = .0625;
	float k5 = .2;
	float k6 = .2;
	long double k7 = .00000004;
	float k8 = .005;
	float k9 = .005;
	float k10 = .0001;
	float k11 = .0001;
	float k_k = 5000;
	double k_m = 2500;
	int h = 2;
	int p = 2;
	float omega = 1;
	int moleculesK = 500;
	int moleculesM = 150;
	int mRNA_k = 0;
	int mRNA_m = 0;
	
	
	// time parameters
	double t = 0;
	double tau = 0;
	float time_hours = 10000;
	float tmax = time_hours*3600;
	
	
	// sets random seed for random numbers
  	srand((unsigned)time(NULL));
  	
  	// declaring files to which to write
	ofstream stocfileTKM;
  	
  	// declares array for storing values
	int storeTKMmax = 1000000000;
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


	// length of loop for doing multiple values of k1 in one run
  	while (j<1)
  	{
		while (t<=tmax)
		{  	
		
			if ( m%samplerate == 0)
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
			
			double a1 = k1 + k2*pow(moleculesK,h)/(pow(k_k,h)+pow(moleculesK,h));
			double a2 = k3 + k4*pow(moleculesK,p)/(pow(k_m,p)+pow(moleculesK,p));
			double a3 = k5*mRNA_k;
			double a4 = k6*mRNA_m;
			double a5 = k7/omega*moleculesM*moleculesK;
			double a6 = k8*mRNA_k;
			double a7 = k9*mRNA_m;
			double a8 = k10*moleculesM;
			double a9 = k11*moleculesK;
			
			
			// update sum of propensities
			double a0 = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9;
			
			// determine time of next event
			tau = ((1)/(a0))*(log(1/r1));
	
			// update time
			t = t + tau;
			//determines the outcome of the event
			{
				if (0 <= r2 and r2 < a1/a0) 
				{
					mRNA_k += 1;
				}
				else if (a1/a0 <= r2 and r2 < (a1+a2)/a0) 
				{
					mRNA_m += 1;
				}
				else if ((a1+a2)/a0 <= r2 and r2 < (a1+a2+a3)/a0) 
				{
					moleculesK += 1;
				}
				else if ((a1+a2+a3)/a0 <= r2 and r2 < (a1+a2+a3+a4)/a0)
				{
					moleculesM += 1;
				}
				else if ((a1+a2+a3+a4)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5)/a0) 
				{
					moleculesK -= 1;
				}
				else if ((a1+a2+a3+a4+a5)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5+a6)/a0) 
				{
					mRNA_k -= 1;
				}
				else if ((a1+a2+a3+a4+a5+a6)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5+a6+a7)/a0)
				{
					mRNA_m -= 1;
				}
				else if ((a1+a2+a3+a4+a5+a6+a7)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5+a6+a7+a8)/a0)
				{
					moleculesM -= 1;
				}
				else if ((a1+a2+a3+a4+a5+a6+a7+a8)/a0 <= r2 and r2 < (a1+a2+a3+a4+a5+a6+a7+a8+a9)/a0)
				{
					moleculesK -=1;
				}
				else
				{
					cout << "ERROR" << "\n";
				}
			}
		//	cout << a1 << ", " << a2 << ", " << a3 << ", " << a4 << ", " << a5 << ", " << a6 << ", " << a7 << ", " << a8 << ", " << a9 <<"\n";
		//	cout << t << ",     " << n << ",     " << moleculesK << ",     " << moleculesM << ",     " << test << "\n";
		}
		
		// Converts k1 to string
		stringstream ss;
		ss << fixed << setprecision(5) << k1;
		string string_k1(ss.str());
		string_k1.erase(0,2);
		
		string filename_k1 = "gsf_TKM_" + string_k1 + ".txt";
		
		
		
		//writing the t, K, and S values to file
		stocfileTKM.open(filename_k1.c_str(), ios::trunc);
		stocfileTKM.close();
		
		stocfileTKM.open (filename_k1.c_str(), ios::app);
		
		while (i < n)
		{
			if (storeT[i] >= 0)
			{
				stocfileTKM << fixed << setprecision(15) << storeT[i] << ", " << storemoleculesK[i] << ", " << storemoleculesM[i] << ",\n";
			}
			i += 1;
		}
		stocfileTKM.close();
		
		
		k1 = k1 + stepsize_k1;
		j += 1;
		
		t = 0;
		tau = 0;
		
		moleculesK = 0;
		moleculesM = 0;
		mRNA_k = 0;
		mRNA_m = 0;
		
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
