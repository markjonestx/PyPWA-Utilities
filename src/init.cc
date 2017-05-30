#include <fit.h>
#include <fitparse.h>
#include <math.h>

extern complex<double> Log(complex<double>), Log10(complex<double>), Exp(complex<double>), Sqrt(complex<double>), integer(complex<double>), Pow(complex<double>,complex<double>), AbsSq(complex<double>), Conj(complex<double>), Real(complex<double>), Imag(complex<double>);

static struct {
	char *name;
	double cval;
} consts[] = {
	{"PI",	3.14159265358979323846},
	{"E",	2.71828182845904523536},
	{"GAMMA",	0.57721566490153286060},
	{"DEG",	57.29577951308232087680},
	{"PHI",	1.61803398874989484820},
	{"fcn",	0.0},
	{"nevents",	0.0},
	{0,	0}
};

static struct {
	char *name;
	complex<double> (*func)(complex<double>);
} builtins[] = {
	{"absSq", AbsSq},
	{"log", Log},
	{"conj", Conj},
	{"real", Real},
	{"imag", Imag},
	{"sqrt", Sqrt},
	{"print", Print},
	{0,	0}
};

void init() {
	int i;
	Symbol *s;

	for (i = 0; consts[i].name; i++)
		install(consts[i].name, CTYPE, consts[i].cval);
	for (i = 0; builtins[i].name; i++) {
		s = install(builtins[i].name, BLTIN, 0.0);
		s->ptr = builtins[i].func;
	}
}
