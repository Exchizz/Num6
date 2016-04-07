#include "nr3.h"
#include "math.h"
#include "quadrature.h"
//#include "dftintegrate.h"
//#include "fasper.h"
//#include "fourier.h"
#include "interp_1d.h"
#include "romberg.h"
#include "derule.h"
using namespace std;


template<class T>
Doub qDErule(T &func, const Doub a, const Doub b, const Doub eps=1.0e-10) {
        const Int JMAX=200;
        Doub s,olds=0.0;
        DErule<T> t(func,a,b);
        for (Int j=0;j<JMAX;j++) {
                s=t.next();
		std::cout << "s: " << s << std::endl;
                if (j > 5)
			std::cout << "diff: " << abs(s-olds) <<  " vs. " << eps*abs(olds) << std::endl;
                        if (abs(s-olds) < eps*abs(olds) ||
                                (s == 0.0 && olds == 0.0)) return s;
                olds=s;
        }
        throw("Too many steps in routine qDErule");
}


Doub F1(Doub x){
        return cos(x*x)*exp(-x);
}

Doub F1_derule(Doub x, Doub delta){
        return cos(x*x)*exp(-x);
}

Doub F2(Doub x){
        return sqrt(x)*cos(x*x)*exp(-x);
}

Doub F2_derule(Doub x, Doub delta){
        return sqrt(x)*cos(x*x)*exp(-x);
}


Doub F3(Doub x){
        return 1000*exp(-1/x)*exp(-1/(1-x));
}

Doub F3_derule(Doub x, Doub delta){

	if(x < 0.1){
	        return 1000*exp(-1/(x+delta))*exp(-1/(1-x));
	} else if(x > 0.9){
	        return 1000*exp(-1/x)*exp(-1/(1-(x-delta)));
	} else {
	        return 1000*exp(-1/x)*exp(-1/(1-x));
	}

/*
	if(x < 0.1){ // x == 0
	        return 1000*exp(-1/(delta))*exp(-1/(1-x));
	} else if(x > 0.9){ // x == 1
	        return 1000*exp(-1/x)*exp(-1/(1-(delta)));
	} else {
	        return 1000*exp(-1/x)*exp(-1/(1-x));
	}
*/
}

Doub F4(Doub x){
        return 1/sqrt(x)*cos(x*x)*exp(-x);
}

Doub F4_derule(Doub x, Doub delta){
        return 1/sqrt(x)*cos(x*x)*exp(-x);
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
        std::cout << "Result qpromo: " << result << "\t\t" << std::flush;
	result = qDErule(F1_derule, 0, 1);
        std::cout << "Result DErule: " << result << std::endl;


        result = qmidt(F2, 0, 1);
        std::cout << "Result midt A(h): " << result << "\t";
        result = qtrap(F2, 0, 1);
        std::cout << "Result trap A(h): " << result << "\t";
        result = qsimp(F2, 0, 1);
        std::cout << "Result sim A(h): " << result << "\t";
        Midpnt<Doub(Doub)> q2(F2, 0, 1);
        result = qromo(q2);
        std::cout << "Result qpromo: " << result << "\t\t";

	result = qDErule(F2_derule, 0, 1);
        std::cout << "Result DErule: " << result << std::endl;


        result = qmidt(F3, 0, 1);
        std::cout << "Result midt A(h): " << result << "\t";
        result = qtrap(F3, 0, 1);
        std::cout << "Result trap A(h): " << result << "\t";
        result = qsimp(F3, 0, 1);
        std::cout << "Result sim A(h): " << result << "\t";
        Midpnt<Doub(Doub)> q3(F3, 0, 1);
        result = qromo(q3);
        std::cout << "Result qpromo: " << result << "\t\t";
	result = qDErule(F3_derule, 0, 1);
        std::cout << "Result DErule: " << result << std::endl;



	result = qDErule(F4_derule, 0, 1);
        std::cout << "Result DErule: " << result << "\t" << std::endl;

        result = qmidt(F4, 0, 1);
        std::cout << "Result midt A(h): " << result << "\t" << std::flush;
        result = qtrap(F4, 0, 1);
        std::cout << "Result trap A(h): " << result << "\t";
        result = qsimp(F4, 0, 1);
        std::cout << "Result sim A(h): " << result << "\t";
        Midpnt<Doub(Doub)> q4(F4, 0, 1);
        result = qromo(q4);
        std::cout << "Result qpromo: " << result << std::endl;

}
