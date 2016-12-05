/*
|------------------------------------------------------------------|
| Final Project for ORIE 5600, 2016 Fall, by Prof. Minca           |
| Finish Date: 2016-11-29                                          |
| Copyright (C) 2016  Yuheng Zhou                                  |
|                                                                  |
| This program is free software; you can redistribute it and/or    |
| modify it under the terms of the GNU General Public License      |
| as published by the Free Software Foundation; either version 2   |
| of the License, or (at your option) any later version.           |
|                                                                  |
| This program is distributed in the hope that it will be useful,  |
| but WITHOUT ANY WARRANTY; without even the implied warranty of   |
| MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    |
| GNU General Public License for more details.                     |
|                                                                  |
| You should have received a copy of the GNU General Public License|
| along with this program; if not, write to the Free Software      |
| Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA    |
| 02110-1301, USA.                                                 |
|------------------------------------------------------------------|
*/

#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>

#define Pi 3.141592653589793238462643

using namespace std;

class Option
{
public:
	double selft,T,alpha;
	double t,x,K,r,sigma,d1,d2;

	/* ------------------------------------------------------------
		Constructor of Option class.
		Use Initializing list to put value into member variables.
	--------------------------------------------------------------*/
	Option(double t, double x, double T, double K, double r, double sigma, 
		double alpha):t(t),x(x),T(T),K(K),r(r),sigma(sigma),alpha(alpha){};

	/* ------------------------------------------------------------
	 	Price a call option using information in member variables
		- return value: price of call
	--------------------------------------------------------------*/
	double Call(){
		d1 = (log(x/K)+(r+0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t));
		d2 = d1-sigma*sqrt(T-t);
		return x*cdf(d1)-cdf(d2)*K*exp(-r*(T-t));
	}

	/*--------------------------------------------------------------
	 	Calculate delta for any point in the path.
	 	- input t: current time in the path
	 	- input St: current price of stock
	 	- return value: delta
	--------------------------------------------------------------*/
	double Delta(double t, double St){
		double d = (log(St/K)+(r+0.5*sigma*sigma)*(T-t))/(sigma*sqrt(T-t));
		return cdf(d);
	}

	/*--------------------------------------------------------------
	 	Generate 1000 paths to conduct monte carlo simulation
	 		with drift term alpha
	 	- return value: expected value of payoff at time T
	--------------------------------------------------------------*/
	double ExpectedDiscountedPayoffP(){
		int N=1000, NofSample=100;
		default_random_engine generator;
		normal_distribution<double> distribution(alpha*(T-t)/N, sigma*sqrt(1.0*
			(T-t)/N));
		double S, dW, res=0.0;
		for (int i=0; i<NofSample; i++){
			S = x;
			for (int j=0; j<N; j++){
				dW = distribution(generator);
				S += dW*S;
			}
			res += max(S-K,0.0)*exp(-r*(T-t));
		}
		return res/NofSample;
	}

	/*--------------------------------------------------------------
	 	Generate 1000 paths to conduct monte carlo simulation
	 		with drift term r
	 	- return value: expected value of payoff at time T
	--------------------------------------------------------------*/
	double ExpectedDiscountedPayoffPTilde(){
		int N=1000, NofSample=100;
		default_random_engine generator;
		normal_distribution<double> distribution(r*(T-t)/N, sigma*sqrt(1.0*
			(T-t)/N));
		double S, dW, res=0.0;
		for (int i=0; i<NofSample; i++){
			S = x;
			for (int j=0; j<N; j++){
				dW = distribution(generator);
				S += dW*S;
			}
			res += max(S-K,0.0)*exp(-r*(T-t));
		}
		return res/NofSample;
	}

	/*--------------------------------------------------------------
	 	Compare the final value of options and delta hedge
	 	- input n: number of steps
	 	- return value: the hedging error
	--------------------------------------------------------------*/
	vector<double> hedge(int n){
		double X,now,S,delta,dW;
		vector<double> vX, result;
		double aveX = 0.0, sdX = 0.0;
		int NofSample=2000;

		default_random_engine generator;
		normal_distribution<double> distribution(alpha*(T-t)/n, sigma*sqrt(1.0*
			(T-t)/n));
		for (int k=1; k<NofSample; k++){
			delta = 0.0;
			S = x;
			X = Call();
			for (int i=0; i<n; i++){
				now = 1.0*T*i/n;
				dW = distribution(generator);
				X = delta*(1+dW)*S+(X-delta*S)*exp(r*T/n);
				S += dW*S;
				delta = Delta(now,S);
			}
			vX.push_back(X-(max(S-K,0.0)));
			//cout << X << " " << max(S-K,0) << endl;
		}
		aveX = accumulate(begin(vX), end(vX), 0.0);
		aveX /= vX.size();
		vector<double>::iterator it;
		for (it=vX.begin(); it!=vX.end(); it++)
			sdX += pow((*it)-aveX, 2);
		sdX /= vX.size();
		result.push_back(aveX);
		result.push_back(sdX);
		return result;
	}

	/*--------------------------------------------------------------
	 	Compute Phi(), the cdf of normal random variables
	 	- input x: parameter of cdf
	 	- return value: the value of normal cdf
	 	sourse: http://stackoverflow.com/questions/2328258/cumulative
	 	    -normal-distribution-function-in-c-c
	 	Make it incline to improve efficiency.
	--------------------------------------------------------------*/
	inline double cdf(double x){
  		double L,K,w;
  		double const a1=0.31938153, a2=-0.356563782, a3=1.781477937;
  		double const a4=-1.821255978, a5=1.330274429;
  		L = fabs(x);
  		K = 1.0/(1.0+0.2316419*L);
  		w = 1.0-1.0/sqrt(2*Pi)*exp(-L*L/2)*(a1*K+a2*K*K+a3*pow(K,3)
  			+a4*pow(K,4)+a5*pow(K,5));
  		if (x<0) w = 1.0-w;
  		return w;
	}

};

int main(){
	ofstream myfile;
	myfile.open("output.csv");
	Option op(0,100,1,100,0.03,0.2,0.08);
	vector<double> v;
	cout << "The call price is: " << op.Call() << endl;
	//cout << "Result of Simulation: " << op.ExpectedDiscountedPayoffPTilde() << endl;
	/*for (int n=1; n<50; n++){
		v = op.hedge(n);
		myfile << -v[0] << "," << v[1] << endl;
	}*/
	myfile.close();
	return 0;
}