#include "nr3.h"
#include "ludcmp.h"
#include "qrdcmp.h"
#include "roots_multidim.h"
#include "math.h"

using namespace std;
#define d_value 30
#define n_value 5
#define a_value 40
#define h_value 500
#define k_value 2.5
#define v_value 120
#define w_value 4
#define alpha_value   2*pow(10,-7)

#define x     0
#define p     1
#define theta 2
#define phi   3
#define l0    4
#define l     5
#define a     6
#define h     7

VecDoub F(VecDoub_I X){
	VecDoub ret(8);


	ret[0] = X[a]*(cosh(X[x]/X[a]) -1) - X[p];
	ret[1] = 2*X[a]*sinh(X[x]/X[a]) - X[l];
	ret[2] = 2*X[x]+2*k_value*cos(X[theta]) - d_value;
	ret[3] = X[p]+k_value*sin(X[theta]) - n_value;
	ret[4] = sinh(X[x]/X[a]) - tan(X[phi]);
	ret[5] = (1+v_value/(w_value*X[l0]))*tan(X[phi]) - tan(X[theta]);
	ret[6] = X[l0]*(1+alpha_value*X[h])-X[l];
	ret[7] = ( w_value*X[l0] )/( 2*sin(X[phi]) ) - X[h];


/*	for(int i = 0; i < 7; i++){
		std::cout << ret[i] << "\t";
	}
	std::cout << std::endl;
*/        return ret;
}


/*
	ret[x] = d_value/2 - k_value*cos( X[theta] )-X[x]; //12
	ret[p] = X[a] *(cosh(X[x]/X[a])-1); //10
	ret[theta] = asin((n_value-X[p])/k_value); //13
	ret[phi] = atan( tan( X[theta] / (1+v_value/(w_value*X[l0])) )); //15
	ret[l0] = X[l]/( 1 + alpha_value*h_value ); //16
	ret[l] = 2*X[a]*sinh( X[x]/X[a] ); //11
	ret[a] = X[x]/( asinh( tan(X[phi]) ) ); //14
	ret[h] = (w_value*X[l0])/(2*sin(X[phi])); //17
	return ret;
*/


int main() {
//void newt(VecDoub_IO &x, Bool &check, T &vecfunc) {

	VecDoub X(8);

	// X
	X[0] = 20;

	// P
	X[1] = 2;

	// Theta (deg)
	X[2] = 20*(M_PI/180);

	// phi (deg)
	X[3] = 5*(M_PI/180);

	// L0 [M]
	X[4] = 25;

	// L [M]
	X[5] = 27;

	// a ()
	X[6] = 40;

	// h
	X[7] = 500;



	Bool done_status;
	newt(X, done_status, F);

	for(int i = 0; i < 7; i++){
		std::cout << X[i] << "\t";
	}
	std::cout << std::endl;

}
