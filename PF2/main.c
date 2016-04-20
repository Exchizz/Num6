#include <string>
#include <iostream>
#include <fstream>
#include "nr3.h"
#include "cholesky.h"
#include "svd.h"

using namespace std;

#define EPSILON_1 0.80
#define EPSILON_2 0.60

#define T1 1000
#define T2 500

#define SIGMA 1.712e-9

#define d 1
#define w 1

#define ALPHA1 EPSILON_1*SIGMA*pow(T1,4)
#define ALPHA2 (1-EPSILON_1)

#define BETA1 EPSILON_2*SIGMA*pow(T2,4)
#define BETA2 (1-EPSILON_2)

double F(double x, double y, double d_param){
	return 1.0/( 2* pow( d_param*d_param+ pow(x-y,2) ,1.5) );
}



// Number of trapetz
int n = 4;

// numper of points
int Npoints = n + 1;

int main() {
	double a = -0.5;
	double b = 0.5;

	double h =(b-a)/n;

	/* Create Matrix A */
	MatDoub A(2*Npoints,2*Npoints);

	/* Create vector x */
        VecDoub X(2*Npoints);

	/* Create vector */
        VecDoub B(2*Npoints);
	for(int i = 0; i < Npoints; i++) {
		B[i] = ALPHA1;
	}

	for(int i = Npoints; i < 2*Npoints;i++) {
		B[i] = BETA1;
	}


	// Create diagonal matrix
        for(int i = 0; i < Npoints * 2; i++) {
		A[i][i] = 1;
        }


	float x = a;
	float y = a;

	float divisor = 0;

	/* Fill upper right */
	for(int row = 0; row < Npoints; row++){
		for(int col = Npoints; col < 2*Npoints; col++){
			divisor = (y == a || y == b) ? 2 : 1;
			std::cout << "divisor: " << divisor << std::endl;
//			A[row][col] = divisor;
			A[row][col] = -ALPHA2*(h/divisor)*F(x,y,d);
//			cout << "x: " << x << " y:" << y << "\t" ;
			y+=h;
		}
//		std::cout << endl;
		y = a;
		x+=h;
	}

	x = a;
	y = a;

	/* Fill lower left */
	for(int row = Npoints; row < 2*Npoints; row++){
		for(int col = 0; col < Npoints; col++){
			divisor = (x == a || x == b) ? 2 : 1;
//			A[row][col] = divisor;
			A[row][col] = -BETA2*(h/divisor)*F(x,y,d);
//			cout << "x:" << x << " y:" << y << "\t" ;
			x+=h;
		}
		x=a;
//		std:cout << endl;
		y+=h;
	}

	cout << "A Matrix: " << endl;
	A.print();

        SVD svd(A);
        svd.solve(B,X);

	cout << "B vector: " << endl;
	B.print();

	cout << "X vector: " << endl;

	ofstream myfile;

	std::ostringstream ss;
	ss << "results_" << n << ".csv";

	myfile.open ( ss.str() );
	myfile << "x:,";

	for(int i = 0; i < Npoints; i++){
		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
			myfile << a+(i*h) << ",";
	}
	myfile << endl;

	myfile << "u(x),";
	for(int i = 0; i < Npoints; i++){
		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
			myfile << X[i] << ",";
	}
	myfile << endl;
	myfile << "v(y),";

	for(int i = 0; i < Npoints; i++){
		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
			myfile << X[Npoints+i] << ",";
	}
	myfile.close();
}
