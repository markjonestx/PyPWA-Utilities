#include <stdio.h>
#include <parseutils.h>
#include <stdlib.h>

extern char* progname;
extern int lineno;

void warning(char *s, char *t)
{
    fprintf(stderr, "%s: %s", progname, s);
    if (t)
        fprintf(stderr, " %s", t);
    fprintf(stderr, " near line %d\n", lineno);
}

void yyerror(char *s)
{
    warning(s, (char *) 0);
    exit(1);
}

void execerror(char *s, char *t)
{
    warning(s, t);
    exit(1);
}

