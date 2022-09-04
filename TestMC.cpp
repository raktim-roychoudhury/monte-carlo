// HardCoded.cpp
//
// C++ code to price an option, essential algorithms.
//
// We take CEV model with a choice of the elaticity parameter
// and the Euler method. We give option price and number of times
// S hits the origin.
//
// (C) Datasim Education BC 2008-2011
//

#include "OptionData.hpp" 
#include "NormalGenerator.hpp"
#include "Range.cpp"
#include <cmath>
#include <iostream>

//Function to check accuracy of MC simulations
vector<double> MCAccuracy (const vector<double>& error,double r ,double t){
    double Nsim = error.size();
    vector<double> accuracy(2,0);
    double callsum = 0, callsumsquared = 0;
    for (int i = 0; i < error.size(); i++){
        callsum += error[i];
        callsumsquared += pow(error[i],2);
    }
    double stddev;
    stddev = sqrt((callsumsquared - (pow(callsum,2)/Nsim)) / (Nsim - 1)) * exp(-r*t);
    accuracy[0] = stddev;
    accuracy[1] = stddev/sqrt(Nsim);
    return accuracy;
}


template <typename T>
void print (const std::vector<T>& myList){
    std::cout << std::endl << "Size of vector is " << myList.size() << "\n[";
    typename std::vector<T>::const_iterator i;
    for (i = myList.begin(); i != myList.end(); ++i){
        std::cout << *i <<",";
    }
    std::cout << "]\n";
    
}

namespace SDEDefinition
{ // Defines drift + diffusion + data

	OptionData* data;				// The data for the option MC

	double drift(double t, double X)
	{ // Drift term
	
		return (data->r)*X; // r - D
	}

	
	double diffusion(double t, double X)
	{ // Diffusion term
	
		double betaCEV = 1.0;
		return data->sig * pow(X, betaCEV);
		
	}

	double diffusionDerivative(double t, double X)
	{ // Diffusion term, needed for the Milstein method
	
		double betaCEV = 1.0;
		return 0.5 * (data->sig) * (betaCEV) * pow(X, 2.0 * betaCEV - 1.0);
	}
} // End of namespace


int main()
{
    std:vector<double> Nvec {100,500,1000}; //vector of time step values
    std::vector<double> Nsimvec {30000}; //vector of Nsim values.
    
	std::cout <<  "1 factor MC with explicit Euler\n";
//    Batch 1
    OptionData myOption;
    myOption.K = 65.0;
    myOption.T = 0.25;
    myOption.r = 0.08;
    myOption.sig = 0.3;
    myOption.type = -1;    // Put -1, Call +1
    double S_0 = 60;
    string batch = "Batch 1";
    double call_price = 2.13337, put_price = 5.84628;
    
//    Batch 2
//    OptionData myOption;
//    myOption.K = 100.0;
//    myOption.T = 1;
//    myOption.r = 0.0;
//    myOption.sig = 0.2;
//    myOption.type = 1;    // Put -1, Call +1
//    double S_0 = 100;
//    string batch = "Batch 2";
//    double call_price = 7.96557, put_price = 7.96557;
    
//    Batch 3
//    OptionData myOption;
//    myOption.K = 10.0;
//    myOption.T = 1;
//    myOption.r = 0.12;
//    myOption.sig = 0.5;
//    myOption.type = 1;    // Put -1, Call +1
//    double S_0 = 5;
//    string batch = "Batch 3";
//    double call_price = 0.204058, put_price = 4.07326;
    
//    Batch 4
//    OptionData myOption;
//    myOption.K = 100.0;
//    myOption.T = 30;
//    myOption.r = 0.08;
//    myOption.sig = 0.3;
//    myOption.type = 1;    // Put -1, Call +1
//    double S_0 = 100;
//    string batch = "Batch 4";
//    double call_price = 92.17570, put_price = 1.24750;

    std::cout << batch << endl;

    for (int i = 0; i < Nsimvec.size(); i++){
        for(int j = 0; j < Nvec.size(); j++){
    
//    long N = 500;
	long N = Nvec[j];
//	std::cout << "Number of subintervals in time: ";
//	std::cin >> N;

	// Create the basic SDE (Context class)
	Range<double> range (0.0, myOption.T);
	double VOld = S_0;
    double VNew = 0;

	std::vector<double> x = range.mesh(N);
	

	// V2 mediator stuff
//	long NSim = 50000;
    long NSim = Nsimvec[i];
//	std::cout << "Number of simulations: ";
//	std::cin >> NSim;

	double k = myOption.T / double (N);
	double sqrk = sqrt(k);

	// Normal random number
	double dW;
	double price = 0.0;	// Option price

	// NormalGenerator is a base class
	NormalGenerator* myNormal = new BoostNormal();

	using namespace SDEDefinition;
	SDEDefinition::data = &myOption;

	std::vector<double> res;
	int coun = 0; // Number of times S hits origin

    std::vector<double> error; //Store price values to check accuracy
    
	// A.
	for (long i = 1; i <= NSim; ++i)
	{ // Calculate a path at each iteration
			
//		if ((i/10000) * 10000 == i)
		{// Give status after each 1000th iteration

//				std::cout << i << std::endl;
		}

		VOld = S_0;
		for (unsigned long index = 1; index < x.size(); ++index)
		{

			// Create a random number
			dW = myNormal->getNormal();
				
			// The FDM (in this case explicit Euler)
			VNew = VOld  + (k * drift(x[index-1], VOld))
						+ (sqrk * diffusion(x[index-1], VOld) * dW);

			VOld = VNew;

			// Spurious values
			if (VNew <= 0.0) coun++;
		}
			
		double tmp = myOption.myPayOffFunction(VNew);
        error.push_back(tmp); //adding price for each iteration in vector error.
        
		price += (tmp)/double(NSim);
	}
	
	// D. Finally, discounting the average price
	price *= exp(-myOption.r * myOption.T);
    
    vector<double> accuracy = MCAccuracy(error, myOption.r, myOption.T); //error function

	// Cleanup; V2 use scoped pointer
	delete myNormal;

    std::cout << "Time Steps = " << N << " Nsim = " << NSim << endl;
	std::cout << "Price, after discounting: " << price << ", " << std::endl;
//    std::cout << "Actual Price = " << call_price << std::endl;
//    if (myOption.type == 1) std::cout << "Error = " << (abs(price - call_price)/call_price)*100 << "%" << endl;
//    else std::cout << "Error = " << (abs(price - put_price)/put_price)*100 << "%" << endl;
    std::cout << "Standard Deviation = " << accuracy[0] << " , Standard Error = " << accuracy[1] << endl;
//	std::cout << "Number of times origin is hit: " << coun << endl;
    }
    }
	return 0;
}

