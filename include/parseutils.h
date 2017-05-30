#ifndef PARSE_UTILS_H
#define PARSE_UTILS_H

int yyparse(void);
void warning(char *s, char *t);
void yyerror(char *s);
void execerror(char *s, char *t);

#endif
