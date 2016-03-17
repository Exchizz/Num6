#include <iostream>
#include <fstream>
#include "nr3.h"
#include "cholesky.h"
#include "svd.h"

using namespace std;
int uvw_print = false;
int main() {
/*
        VecDoub theta1(500);
	VecDoub theta2(500);
	VecDoub x(500);
	VecDoub y(500);
*/
        MatDoub A(1000,4);
        VecDoub z(1000);
        VecDoub q(4);


        ifstream Pont("d2_data");

	double theta1 = 0;
	double theta2 = 0;
	double x = 0;
	double y = 0;


        for(int i = 0; i < 500; i++) {
                Pont >> theta1;
                Pont >> theta2;
                Pont >> x;
                Pont >> y;

		A[i*2][0] = 1;
		A[i*2][1] = 0;
		A[i*2][2] = cos(theta1);
		A[i*2][3] = cos(theta1 + theta2);

		A[i*2+1][0] = 0;
		A[i*2+1][1] = 1;
		A[i*2+1][2] = sin(theta1);
		A[i*2+1][3] = sin(theta1 + theta2);

		z[i*2] = x;
		z[i*2+1] = y;

        }

	SVD svd(A);

	svd.solve(z,q);
	//q.print();
	std::cout << "x0: " << q[0] << " y0: " << q[1] << " a: " << q[2] << " b: " << q[3] << std::endl;
	VecDoub res(1000);
	res = A*q-z;
	std::cout << "Residual:" << res.length() << std::endl;

	if(uvw_print){
		std::cout << "----------U-matrix:------------" << std::endl;
		svd.u.print();
		std::cout << "----------V-matrix:------------" << std::endl;
		svd.v.print();
		std::cout << "----------w-vector:------------" << std::endl;
		svd.w.print();
	}

	std::cout << "rank: " << svd.rank(0.05) << std::endl;
	std::cout << "nullity: " << svd.nullity(0.05) << std::endl;

	for(int j = 0; j < 4; j++){
		double sum = 0;
		for(int i = 0; i < 4; i++){
			sum += pow( svd.v[j][i] / svd.w[i],2 );
		}
		std::cout << "Uncertainty:(" << j << ") standard diviation: " << sqrt(sum) << " variance: " << sum << std::endl;
	}

//	A.print();
//	z.print();

}
