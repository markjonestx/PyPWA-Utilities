#include <minuit.h>
#include <cstdlib>

extern "C" {
    void mnintr_(...);
    void mninit_(...);
    void mnparm_(...);
    void mnexcm_(...);
    void mnseti_(...);
    void mnemat_(...);
}

// from minuit91.f
//     D/MN7NPR/ MAXINT ,NPAR   ,MAXEXT ,NU
typedef struct {
    int maxInternal;
    int nInternal;
    int maxExternal;
    int nExternal;
} mn7npr_t;

extern "C" {
    mn7npr_t mn7npr_;
}

void mninit(int ird, int iwr, int isav) {
    mninit_(&ird,&iwr,&isav);
}

void mnintr(void (*func)(int*,double*,double*,double*,int*)) {
    mnintr_(func,0);
}

void mnseti(string title) {
    char* ctitle = (char*) malloc (sizeof(char)*title.size() + 1);
    strcpy(ctitle,title.c_str());

    mnseti_(ctitle,strlen(ctitle));
}

int mnparm(int p, string pn, double start, double step, double lbound, double ubound) {

    int ierr;
    char* pname = (char*) malloc (sizeof(char)*pn.size() + 1);
    strcpy(pname,pn.c_str());

    mnparm_(&p,pname,&start,&step,&lbound,&ubound,&ierr,strlen(pname));

    return(ierr);

}

int mnexcm(void (*func)(int*,double*,double*,double*,int*), string command,double* args, int narg) {

    int ierr;
    const char* com = command.c_str();

    mnexcm_(func,com,args,&narg,&ierr,0,strlen(com));

    return(ierr);
}

// from the minuit code:
//        SUBROUTINE MNEMAT(EMAT,NDIM)
//  C ************ DOUBLE PRECISION VERSION *************
//      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
//        DIMENSION EMAT(NDIM,NDIM)
//  CC        Calculates the external error matrix from the internal
//  CC        to be called by user, who must dimension EMAT at (NDIM,NDIM)

matrix<double> minEmat(int size) {
    matrix<double> emat;
    double *memat = (double*) malloc(sizeof(double)*size*size);

    mnemat_(memat,&size);
    emat = matrix<double>(size,size,memat);
    return emat;
}

int nMinpars(){
    return mn7npr_.nExternal;
}

