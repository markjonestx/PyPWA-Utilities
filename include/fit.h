#include <complex>
#include <ext/hash_map>
#include <stdlib.h>
#include <integral.h>

namespace std {
  using namespace __gnu_cxx;
}

using namespace std;

#ifndef FIT_H
#define FIT_H

typedef struct parameter {
  complex<double> val;
  int real_index;
  int imag_index; 
} parameter;

typedef struct Symbol {
  char *name;
  short type;
  parameter par;
  complex<double> val;
  complex<double> (*ptr)(complex<double>);
  matrix<double> rmatval;
  integral normInt;
  struct Symbol *next;
} Symbol;

Symbol *install(char *s, int t, complex<double> d);
Symbol *lookup(char *s);

void print_stable();
void read_amps(void);
void open_ampfiles(void);
void open_rmatfiles(void);
int update_symbols(void);
int update_amps(void);
int update_amps_from_files(void);
int update_rmat_from_files(void);
void rewind_amps(void);
void rewind_rmat(void);
void read_rmat(void);
int update_rmat(void);
void set_params();
void set_param(string pname, parameter p);
void update_params(double *minval);
void print_params(ofstream& os);


struct eqstr {
    bool operator()(const char* s1, const char* s2) const {
        return strcmp(s1,s2)==0;
    }
};

hash_map<const char*, integral, hash<const char*>, eqstr> get_integrals();
hash_map<const char*, parameter> get_params();
void print_integrals(ofstream& os);

void set_integral(string s, integral ni);
void read_integral(string s);
void read_integral(string s, string file);
integral get_integral(string s);
complex<double> get_integral_val(string s, string wv1, string wv2);
void print_integral(string s);

void init();

typedef struct Datum {
	complex<double> val;
  //        matrix<double> rmatval;
	Symbol *sym;
} Datum;

extern Datum pop();

typedef int (*Inst)();
Inst* code(Inst);
#define STOP (Inst) 0

extern Inst event_prog[];
extern Inst norm_prog[];

int init_event_code();
int init_norm_code();
void execute(Inst*);

extern void eval();
extern void evalptype();
extern void evalrmat();
extern void evalnorm();
extern void evalnormev();
extern void add();
extern void sub();
extern void mul();
extern void divide();
extern void Negate();
extern void power();
extern void assign();
extern void bltin();
extern void varpush();
extern void constpush();
extern void print();

complex<double> Pow(complex<double>, complex<double>);
complex<double> Print( complex<double>);
void execerror(char* s, char* t);

#endif
