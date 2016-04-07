#include "nr3.h"
#include "math.h"
#include "quadrature.h"
//#include "dftintegrate.h"
//#include "fasper.h"
//#include "fourier.h"
#include "interp_1d.h"
#include "romberg.h"
#include <vector>
using namespace std;


//VecDoub F(VecDoub_I x){
Doub F1(Doub x){
	return cos(x*x)*exp(-x);
}

Doub F2(Doub x){
	return sqrt(x)*cos(x*x)*exp(-x);
}

int main() {
	Doub result = 0;

	result = qmidt(F1, 0, 1);
	std::cout << "Result midt A(h): " << result << "\t";
	result = qtrap(F1, 0, 1);
	std::cout << "Result trap A(h): " << result << "\t";
	result = qsimp(F1, 0, 1);
	std::cout << "Result sim A(h): " << result << "\t";
	Midpnt<Doub(Doub)> q1(F1, 0, 1);
	result = qromo(q1);
	std::cout << "Result qpromo: " << result << std::endl;


	result = qmidt(F2, 0, 1);
	std::cout << "Result midt A(h): " << result << "\t";
	result = qtrap(F2, 0, 1);
	std::cout << "Result trap A(h): " << result << "\t";
	result = qsimp(F2, 0, 1);
	std::cout << "Result sim A(h): " << result << "\t";
	Midpnt<Doub(Doub)> q2(F2, 0, 1);
	result = qromo(q2);
	std::cout << "Result qpromo: " << result << std::endl;
}
