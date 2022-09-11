#pragma once

#include <complex>
#include <shared/particle.h>


#ifndef YYSTYPE
typedef union {
    int num;
    std::complex<double> *Cnum;
    double Fnum;
    char string[100];
    Decay *Decay;
    particle *Particle;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define    INT    257
# define    FLOAT    258
# define    STRING    261
# define    PARTICLENAME    262
# define    DEBUG    263
# define    CHANNEL    264
# define    MODE    265
# define    MASSDEP    266
# define    HELICITY    267


extern YYSTYPE keylval;
