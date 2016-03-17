#include "nr3.h"
#include "roots.h"
#include "math.h"
using namespace std;

double f(double x){
	return x-cos(x);
}

class Funcd {
public:
	double operator()(double x){
		return x-cos(x);
	}

	double df(double x){
		return 1+sin(x);
	}
};

int main() {

	std::cout << "Bisection: " << std::endl;
	double result = rtbis(f, 0, M_PI/2, pow(10,-10));
	std::cout << "result: " << result << std::endl;

	std::cout << "Ssecant: " << std::endl;
	result = rtsec(f, 0, M_PI/2, pow(10,-10));
	std::cout << "result: " << result << std::endl;

	std::cout << "False position: " << std::endl;
	result = rtflsp(f, 0, M_PI/2, pow(10,-10));
	std::cout << "result:" << result << std::endl;

	std::cout << "Ridder: " << std::endl;
	result = zriddr(f, 0, M_PI/2, pow(10,-10));
	std::cout << "result:" << result << std::endl;

	Funcd funcd;
	std::cout << "Newton: " << std::endl;
	result = rtnewt(funcd, -M_PI/2, M_PI/2, pow(10,-10));
	std::cout << "Result: " << result << std::endl;

}
