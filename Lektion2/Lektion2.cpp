#include <iostream>
#include <fstream>
#include "nr3.h"
#include "cholesky.h"

using namespace std;

int main() {
	/*
	VecDoub xFilip(82); VecDoub yFilip(82);

	ifstream Filip("src/FilipData.dat");
	for(int i = 0; i < 82; i++) {
		Filip >> yFilip[i];
		Filip >> xFilip[i];
	}
	*/


	VecDoub xPont(40); VecDoub yPont(40);
	MatDoub A(40,3);
	ifstream Pont("PontiusData.dat");
	for(int i = 0; i < 40; i++) {
		Pont >> yPont[i];
		Pont >> xPont[i];
	}


	for(int i = 0; i < 40; i++){
		A[i][0] = 1;
		A[i][1] = xPont[i];
		A[i][2] = xPont[i]*xPont[i];
	}


//	A.print();
        MatDoub AT(40,3);
	MatDoub ATA(40,3);
	VecDoub ATb(40);
	VecDoub x(3);

	AT = A.transpose(); // 3*40

	ATb = AT*yPont; // 3*1

	ATb.print();

	ATA = AT*A;
	Cholesky obj(ATA);
	obj.solve(ATb, x);
	x.print();


return 0;
}
