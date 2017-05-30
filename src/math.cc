#include <fit.h>
#include <math.h>
#include <errno.h>
extern int errno;
complex<double> errcheck(complex<double>, char*);

complex<double> AbsSq(complex<double> x) {
	return x*conj(x);
}

complex<double> Log(complex<double> x) {
  double y = x.real();
  if (fabs(y) > 0.0001)
	return errcheck(log(x.real()), "log");
  else
    return(0.0);
}

complex<double> Conj(complex<double> x) {
	return conj(x);
}

complex<double> Real(complex<double> x) {
	return real(x);
}

complex<double> Imag(complex<double> x) {
	return imag(x);
}

complex<double> Exp(complex<double> x) {
	return errcheck(exp(x), "exp");
}

complex<double> Sqrt(complex<double> x) {
	return errcheck(sqrt(x), "sqrt");
}

complex<double> Print(complex<double> x) {
  cout << x << endl;
	return x;
}

complex<double> Pow(complex<double> x, complex<double> y) {
	return errcheck(pow(x,y), "exponentiation");
}


complex<double> errcheck(complex<double> d, char* s) {
	if (errno == EDOM) {
		errno = 0;
		execerror(s, "argument out of domain");
	} else if (errno == ERANGE) {
		errno = 0;
		execerror(s, "result out of range");
	}
	return d;
}

