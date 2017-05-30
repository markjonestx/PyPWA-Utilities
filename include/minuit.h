#include <string>
#include <matrix.h>

using namespace std;

void mnintr(void (*func)(int*,double*,double*,double*,int*));
void mninit(int, int, int);
void mnseti(string);
int mnparm(int, string, double, double, double, double);
int mnexcm(void (*)(int*,double*,double*,double*,int*), string ,double* , int );
matrix<double> minEmat(int size);
int nMinpars();

extern "C" {
	void fortopen_(const char*, int);
	void fortopeninput_(const char*, int);

}
