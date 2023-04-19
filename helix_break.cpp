#include <iostream>
#include <cstdlib>  
#include <cmath>
#include "cpgplot.h"
#include "radiation.h"
#include "tumour_growth.h"

int main (){
	//Constants for DSB (const >= 0)
	//	delta is DSB per cell
	//	tau is the repair time constant
	//	gamma is binary reaction rate
	float delta = 2.0;
	float tau = 1.1;
	float gamma = 1.3;

	//Constants Cell Growth (const >= 0)
	//	alpha is the number of unrepairable cells (alpha <= gamma/2)
	//	2*kappa is the rate of lethal misrepair per DSB pair
	float alpha = 0.1;
	float kappa = 0.3;

	float rad = 2.0;

	//n is cell count
	//tg is tumour growth
	float n=50000., tg_val=23.;
	float n1, n2, n3, n4;
	float tg1, tg2, tg3, tg4;

	//CASES
	int iCases;
	std::cout << "Insert numbers for test cases. Options are 1, 2, and 3." << "\n";
	std::cin >> iCases;
	
	switch(iCases) 
	{
		case 1:
			delta = 2.0;
			tau = 1.1;
			gamma = 1.3;
			alpha = 0.1;
			kappa = 0.3;
			rad = 2.0;
			break;
		case 2:
			delta = 0.2;
			tau = 0.1;
			gamma = 1.3;
			alpha = 0.1;
			kappa = 0.3;
			rad = 0.2;
			break;
		case 3:
			delta = 0.0;
			tau = 3.1;
			gamma = 2.3;
			alpha = 0.1;
			kappa = 0.3;
			rad = 0.0;
			break;
		default:
			std::cout << iCases << " is not an allowed choice\n";
	}

	//np represents numbers of cells at given time
	//rp represents DHB average at given time
	//tp represents time
	int steps = 1048576.;
	float np[steps], tg[steps], tp[steps];

	// (Initial)
	np[0] = n;
	tg[0] = tg_val;

	//Time Variables
	double t_start = 0;
	double t_end = .5;
	double t;
	double t_delta = (t_end-t_start)/steps;

	// Graph Initializations
	if (!cpgopen("/XWINDOW")) return 1;
	cpgenv(0.,0.5,0.,12000.,0,1);
	cpglab("time", "cell count & DSB", "ODEs Cell Count over Time with Radiation");

	//Runge Kutta Integration
	for(int i=1; i<steps+1; i++) 
	{
       		tg1 = t_delta * tumour_growth(tg_val, rad, delta, tau, gamma);
		n1 = t_delta * radiation(tg_val, n, rad, alpha, kappa);

       		tg2 = t_delta * tumour_growth(tg_val+(tg1/2), rad, delta, tau, gamma);
		n2 = t_delta * radiation(tg_val+(tg1/2), n, rad, alpha, kappa);

       		tg3 = t_delta * tumour_growth(tg_val+(tg2/2), rad, delta, tau, gamma);
		n3 = t_delta * radiation(tg_val+(tg2/2), n, rad, alpha, kappa);

       		tg4 = t_delta * tumour_growth(tg_val+tg3, rad, delta, tau, gamma);
       		n4 = t_delta * radiation(tg_val+tg3, n, rad, alpha, kappa);


		t += t_delta;		

       		tg_val += (tg1+ 2*tg2 + 2*tg3 + tg4)/6.0;
		n -= (n1+ 2*n2 + 2*n3 + n4)/6.0;
     
		tp[i] = t;
       		tg[i] = tg_val;
       		np[i] = n;
	}

	cpgsci(4);
 	cpgline(steps+1,tp,np);

	cpgsci(6);
	cpgline(steps+1,tp,tg);
	cpgclos();

	return 0;
}
