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

/* Richardson */
#define ALPHA 2
#define K 2

double F(double x, double y, double d_param){
	return 1.0/( 2* pow( d_param*d_param+ pow(x-y,2) ,1.5) );
}

double Q( VecDoub & X, double a, double b , int N, int offset, double k1, double k2){

	double h = (b-a)/N;
	int Npoints = N + 1;

	double sum = 0;

	int divisor = 0;
	for( int elm = 0; elm < Npoints; elm++){
		divisor = (elm == 0 || elm == (Npoints-1)) ? 2 : 1;
		sum += (h/divisor)*(X[offset + elm] - ((k1-X[offset + elm])/-k2));
	}
	return sum;
}



int main(int argc, char** argv) {
	std::vector<double> Q1_values;
	std::vector<double> Q2_values;

	std::vector<double> Q1_values_result;
	std::vector<double> Q2_values_result;

	std::cout.precision(10);

	int current_n = 0;
	for(int richardsson_counter = 2; richardsson_counter > 0; richardsson_counter--){

//		std::cout << "richard counter: " << richardsson_counter << std::endl;
//		std::cout  << "alpha 1: " << ALPHA1 << std::endl;
//		std::cout  << "alpha 2: " << ALPHA2 << std::endl;

		// Number of trapetz
		if(argc == 1){
			std::cout << "HOWTO: <filename> <N>" << std::endl;
			return 1;
		}
		int n = atoi(argv[1]);
		if(argv[1] == ""){
			std::cout << "default N = 4";
			n=4;
		}

		current_n = n;
		n/=richardsson_counter;
//		std::cout << "---------- n: " << n << "------------" << std::endl;

		// numper of points
		int Npoints = n + 1;

		/* Define integral limits */
		double a = -0.5;
		double b = 0.5;

		/* Define stepsize */
		double h = (b-a)/n;

		/* Create Matrix A */
		int zero = 0;
		MatDoub A(2*Npoints,2*Npoints, zero);

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
				A[row][col] = -ALPHA2*(h/divisor)*F(x,y,d);
				y+=h;
			}
			y = a;
			x+=h;
		}

		/* Reset variables */
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

		/* Solve linear system of equations */
//		A.print();
	        SVD svd(A);
	        svd.solve(B,X);

		/* Create file objects */
		ofstream u_x_file;
		ofstream v_y_file;
		ofstream q1_r;
		ofstream q2_r;

		/* Open files for tex files */
		u_x_file.open ( "results_u_x.tex", ios::app );
		v_y_file.open ( "results_v_y.tex", ios::app );

		/* Write contents of X to two seperate files*/
		u_x_file << n << " & ";
		for(int i = 0; i < Npoints; i++){
			if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
				u_x_file << X[i] << " & ";
		}
		v_y_file << n << " & ";
		for(int i = 0; i < Npoints; i++){
			if(a+(i*h) == a || a+(i*h) == -0.25 || a+(i*h) == 0 || a+(i*h) == 0.25 || a+(i*h) == b )
				v_y_file << X[Npoints+i] << " & ";
		}

		/* Calculate Q1 from first Npoints */
		double q1 = Q(X, a, b, n, 0, ALPHA1, ALPHA2);

		/* Calculate Q2 from rest Npoints */
		double q2 = Q(X, a, b, n, Npoints, BETA1, BETA2);

		v_y_file << q1;
		u_x_file << q2;

		u_x_file << "\\\\ \\hline" << std::endl;
		v_y_file << " \\\\ \\hline" << std::endl;
		/* Close files */
		v_y_file.close();
		u_x_file.close();

		Q1_values.push_back(q1);
		Q2_values.push_back(q2);
	}

	//Calculate richardshon based on N/2 and N
	double Q1_R_error = ( Q1_values[1]-Q1_values[0] )/(pow(ALPHA,K)-1);
	double Q2_R_error = ( Q2_values[1]-Q2_values[0] )/(pow(ALPHA,K)-1);

	double Q1_R = Q1_values[1] + Q1_R_error;
	double Q2_R = Q2_values[1] + Q2_R_error;

	// Q1_values[0] = N/2
	// Q1_values[1] = N

	std::cout << current_n <<  " & " << Q1_values[1] << " & " << Q1_R << " & "  << Q1_R_error << " & " << Q2_values[1] << " & " << Q2_R << " & "  << Q2_R_error << "\\\\ \\hline" << std::endl;
//	std::cout << "Q1r: " << Q1_values[1] + (( Q1_values[1]-Q1_values[0] )/(pow(ALPHA,K)-1)) << std::endl;
//	std::cout << "Q2r: " << Q2_values[1] + (( Q2_values[1]-Q2_values[0] )/(pow(ALPHA,K)-1)) << std::endl;
}
