#ifndef BISON_Y_TAB_H
# define BISON_Y_TAB_H

#ifndef YYSTYPE
typedef union {
	Symbol *sym;
	Inst *inst;
	int ival;
} yystype;
# define YYSTYPE yystype
# define YYSTYPE_IS_TRIVIAL 1
#endif
# define	NUMBER	257
# define	VAR	258
# define	BLTIN	259
# define	UNDEF	260
# define	CTYPE	261
# define	RPTYPE	262
# define	PTYPE	263
# define	ATYPE	264
# define	RMTYPE	265
# define	ITYPE	266
# define	STR_INDEX	267
# define	NORMSTATE	268
# define	ELOOPSTATE	269
# define	NEVFIELD	270
# define	INTEGER	271
# define	UNARYMINUS	272


extern YYSTYPE yylval;

#endif /* not BISON_Y_TAB_H */
