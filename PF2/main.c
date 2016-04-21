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




int main(int argc, char** argv) {
	// Number of trapetz
	int n = atoi(argv[1]);
	if(argv[1] == ""){
		std::cout << "default N = 4";
		n=4;
	}


	std::cout << "n: " << n << std::endl;

	// numper of points
	int Npoints = n + 1;

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
			A[row][col] = -ALPHA2*(h/divisor)*F(x,y,d);
			y+=h;
		}
		y = a;
		x+=h;
	}

	x = a;
	y = a;

	/* Fill lower left */
	for(int row = Npoints; row < 2*Npoints; row++){
		for(int col = 0; col < Npoints; col++){
			divisor = (x == a || x == b) ? 2 : 1;
			A[row][col] = -BETA2*(h/divisor)*F(x,y,d);
			x+=h;
		}
		x=a;
		y+=h;
	}

	cout << "A Matrix: " << endl;
	A.print();

        SVD svd(A);
        svd.solve(B,X);

	cout << "B vector: " << endl;
	B.print();

	cout << "X vector: " << endl;

	ofstream u_x_file;
	ofstream v_y_file;

	u_x_file.open ( "results_u_x.csv", ios::app );
	v_y_file.open ( "results_v_y.csv", ios::app );
//	u_x_file << "N  &  -0.5    & -0.25   & 0       & 0.25    & 0.5     &  \\ \cline{1-6}" << std::endl;

//	for(int i = 0; i < Npoints; i++){
//		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
//			u_x_file << ".< & . & . & . & . & . &  \\ \cline{1-6}"";
//	}
//	u_x_file << endl;

	u_x_file << n << " & ";
	for(int i = 0; i < Npoints; i++){
		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
			u_x_file << X[i] << " & ";
	}
	u_x_file << "\\\\ \\cline{1-6}" << std::endl;
	v_y_file << n << " & ";
	for(int i = 0; i < Npoints; i++){
		if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
			v_y_file << X[Npoints+i] << " & ";
	}
	v_y_file << "\\\\ \\cline{1-6}" << std::endl;

	v_y_file.close();
	u_x_file.close();
}
