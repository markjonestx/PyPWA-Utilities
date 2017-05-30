#include <iostream>
#include <fit.h>
#include <fitparse.h>

#define NSTACK 5000
static Datum stack[NSTACK];
static Datum *stackp;

#define NPROG 500000
Inst event_prog[NPROG];
Inst norm_prog[NPROG];
Inst *prog_start;
Inst *progp;
Inst *pc;

int lstack = 0;
int code_debug = 0;

int init_event_code() {
	stackp = stack;
	progp = event_prog;
	prog_start = progp;
	return 0;
}

int init_norm_code() {
	stackp = stack;
	progp = norm_prog;
	prog_start = progp;
	return 0;
}

void push(Datum d) {
	if (stackp >= &stack[NSTACK])
		execerror("stack overflow", (char*) 0);
	*stackp++ = d;
	lstack++;
}

Datum pop() {
	if (stackp <= stack)
		execerror("stack underflow", (char*) 0);
	lstack--;
	return *--stackp;
}

Inst *code(Inst f) {
	Inst *oprogp = progp;
	if (progp >= &prog_start[NPROG])
		execerror("program too big", (char*) 0);
	*progp++ = f;
	return oprogp;
}

void execute(Inst *p) {
	for (pc = p; *pc != STOP;  )
		(*(*pc++))();
}

void constpush() {
	Datum d;
	if(code_debug) cerr << "executing constpush()" << endl;
	d.val = ((Symbol *)*pc++)->val;
	push(d);
	if(code_debug) cerr << "pushing const " << d.val << " onto stack" << endl;
}

void varpush() {
	Datum d;
	if(code_debug) cerr << "executing varpush()" << endl;
	d.sym = (Symbol *)(*pc++);
	push(d);
	if(code_debug) cerr << "pushing var " << d.sym->name << " onto stack" << endl;
	if(code_debug) cerr << "     sym->val: " << d.sym->val << endl;
	if(code_debug) cerr << "          val: " << d.val << endl;
}

void add() {
	Datum d1, d2;
	if(code_debug) cerr << "executing add()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	d1.val += d2.val;
	push(d1);
	if(code_debug) cerr << "pushed " << d1.val << " onto stack" << endl;
}

void sub() {
	Datum d1, d2;
	if(code_debug) cerr << "executing sub()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	d2.val -= d1.val;
	push(d2);
	if(code_debug) cerr << "pushed " << d2.val << " onto stack" << endl;
}

void mul() {
	Datum d1, d2;
	if(code_debug) cerr << "executing mul()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	d1.val *= d2.val;
	push(d1);
	if(code_debug) cerr << "pushed " << d1.val << " onto stack" << endl;
}

void divide() {
	Datum d1, d2;
	if(code_debug) cerr << "executing divide()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	d2.val /= d1.val;
	push(d2);
	if(code_debug) cerr << "pushed " << d2.val << " onto stack" << endl;
}

void Negate() {
	Datum d1;
	if(code_debug) cerr << "executing Negate()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d1.val = -d1.val;
	push(d1);
	if(code_debug) cerr << "pushed " << d1.val << " onto stack" << endl;
}

void power() {
	Datum d1, d2;
	if(code_debug) cerr << "executing power()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.val << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	d2.val = Pow(d2.val, d1.val);
	push(d2);
	if(code_debug) cerr << "pushed " << d2.val << " onto stack" << endl;
}

void eval() {
	Datum d;
	if(code_debug) cerr << "executing eval()" << endl;
	if(code_debug) print_stable();
	d = pop();
	if(code_debug) cerr << "popped " << d.sym->name << " from stack" << endl;
	if (d.sym->type == UNDEF)
		execerror("undefined variable", d.sym->name);
	d.val = d.sym->val;
	if(code_debug) cerr << d.sym->name << " evaluates to " << d.val << endl;
	push(d);
	if(code_debug) cerr << "pushed " << d.val << " onto stack" << endl;
}

void evalptype() {
	Datum d;
	if(code_debug) cerr << "executing evalptype()" << endl;
	d = pop();
	if(code_debug) cerr << "popped " << d.sym->name << " from stack" << endl;
	if (d.sym->type == UNDEF)
		execerror("undefined variable", d.sym->name);
	d.val = d.sym->par.val;
	if(code_debug) cerr << d.sym->name << " evaluates to " << d.val << endl;
	push(d);
	if(code_debug) cerr << "pushed " << d.val << " onto stack" << endl;
}

void evalnormev() {
	Datum norm_int, d;
	norm_int = pop();
	d.val = get_integral(norm_int.sym->name).nevents();
	push(d);
}

void evalnorm() {
	Datum norm,index1,index2,d;
	if(code_debug) cerr << "executing evalnorm()" << endl;
	norm = pop();
	if(code_debug) cerr << "popped " << norm.sym->name << " from stack" << endl;
	index2 = pop();
	if(code_debug) cerr << "popped " << index2.sym->name << " from stack" << endl;
	index1 = pop();
	if(code_debug) cerr << "popped " << index1.sym->name << " from stack" << endl;
	if ((norm.sym->type == UNDEF))
		execerror("undefined variable", norm.sym->name);
	d.val = get_integral_val(norm.sym->name,index1.sym->name,index2.sym->name);
	if(code_debug) cerr << norm.sym->name << "[" << index1.sym->name << "," << index2.sym->name << "] evaluates to " << d.val << endl;
	push(d);
	if(code_debug) cerr << "pushed " << d.val << " onto stack" << endl;
}

void evalrmat() {
	Datum mat,index1,index2,d;
	if(code_debug) cerr << "executing evalrmat()" << endl;
	int i,j;
	matrix<double> m;
	mat = pop();
	if(code_debug) cerr << "popped " << mat.sym->name << " from stack" << endl;
	index2 = pop();
	if(code_debug) cerr << "popped " << index2.val << " from stack" << endl;
	index1 = pop();
	if(code_debug) cerr << "popped " << index1.val << " from stack" << endl;
	if ((mat.sym->type == UNDEF))
		execerror("undefined variable", mat.sym->name);
	i = (int)index1.val.real();
	j = (int) index2.val.real();
	m = mat.sym->rmatval;
	d.val = mat.sym->rmatval.el((int)index1.val.real(),(int) index2.val.real());
	if(code_debug) cerr << mat.sym->name << "[" << i << "," << j << "] evaluates to " << d.val << endl;
	push(d);
	if(code_debug) cerr << "pushed " << d.val << " onto stack" << endl;
}

void assign() {
	Datum d1, d2;
	if(code_debug) cerr << "executing assign()" << endl;
	d1 = pop();
	if(code_debug) cerr << "popped " << d1.sym->name << " from stack" << endl;
	d2 = pop();
	if(code_debug) cerr << "popped " << d2.val << " from stack" << endl;
	// if (d1.sym->type != VAR && d1.sym->type != UNDEF)
	if (d1.sym->type != CTYPE && d1.sym->type != ATYPE && d1.sym->type != PTYPE)
		execerror("assignment to non-variable", d1.sym->name);
	d1.sym->val = d2.val;
	// d1.sym->type = VAR;
	// push(d2);
}

void print() {
	Datum d;
	if(code_debug) cerr << "executing print()" << endl;
	d = pop();
	if(code_debug) cerr << "popped " << d.val << " from stack" << endl;
	cout << "\t" << d.val << endl;;
}

void bltin() {
	Datum d;
	if(code_debug) cerr << "executing bltin()" << endl;
	d = pop();
	if(code_debug) cerr << "popped " << d.val << " from stack" << endl;
	d.val = (*(complex<double> (*)(complex<double>)) (*pc++)) (d.val);
	push(d);
	if(code_debug) cerr << "pushed " << d.val << " onto stack" << endl;
}
