/* A Bison parser, made from keyParse.yy
   by GNU bison 1.35.  */

#define YYBISON 1  /* Identify Bison output.  */

#define yyparse keyparse
#define yylex keylex
#define yyerror keyerror
#define yylval keylval
#define yychar keychar
#define yydebug keydebug
#define yynerrs keynerrs
# define    INT 257
# define    FLOAT   258
# define    COMPLEX 259
# define    SIGN    260
# define    STRING  261
# define    PARTICLENAME    262
# define    DEBUG   263
# define    CHANNEL 264
# define    MODE    265
# define    MASSDEP 266
# define    HELICITY    267

#line 1 "keyParse.yy"

#include <iostream>
#include <particle.h>
#include <wave.h>
#include <keyfile.h>
#include <massDep.h>

#define stoi(x) strcmp((x),"+")?-1:1

using namespace std;

int yylex(void);

void yyerror(char *s);

#define YYDEBUG 1

int nwave;
int debug = 0;
string mode;

wave wv;
extern int lineno;
extern particleDataTable PDGtable;
extern event e;
complex<double> amp;
string t_part_init;
particle *t_part_final;


#line 29 "keyParse.yy"
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
#ifndef YYDEBUG
# define YYDEBUG 0
#endif


#define YYFINAL     89
#define YYFLAG      -32768
#define YYNTBASE    27

/* YYTRANSLATE(YYLEX) -- Bison token number corresponding to YYLEX. */
#define YYTRANSLATE(x) ((unsigned)(x) <= 267 ? yytranslate[x] : 40)

/* YYTRANSLATE[YYLEX] -- Bison token number corresponding to YYLEX. */
static const char yytranslate[] =
        {
                0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                20, 21, 16, 14, 19, 15, 2, 2, 26, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 18,
                2, 17, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 24, 2, 25, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 22, 2, 23, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                2, 2, 2, 2, 2, 2, 1, 3, 4, 5,
                6, 7, 8, 9, 10, 11, 12, 13
        };

#if YYDEBUG
static const short yyprhs[] =
        {
                0, 0, 1, 4, 9, 14, 19, 26, 29, 31,
                35, 39, 43, 47, 51, 54, 56, 58, 62, 66,
                72, 80, 88, 99, 106, 113, 115, 118, 124, 126,
                131, 133, 138, 140, 143, 146, 149
        };
static const short yyrhs[] =
        {
                -1, 27, 28, 0, 9, 17, 3, 18, 0, 10,
                17, 7, 18, 0, 11, 17, 7, 18, 0, 7,
                17, 7, 19, 37, 18, 0, 29, 18, 0, 30,
                0, 20, 29, 21, 0, 29, 14, 29, 0, 29,
                15, 29, 0, 4, 16, 29, 0, 39, 16, 29,
                0, 31, 34, 0, 34, 0, 32, 0, 33, 33,
                33, 0, 7, 17, 3, 0, 22, 35, 35, 3,
                23, 0, 22, 35, 35, 7, 17, 3, 23, 0,
                22, 35, 35, 7, 17, 4, 23, 0, 22, 35,
                35, 7, 17, 3, 7, 17, 3, 23, 0, 22,
                35, 35, 3, 3, 23, 0, 22, 35, 35, 35,
                3, 23, 0, 36, 0, 36, 34, 0, 36, 34,
                12, 17, 7, 0, 37, 0, 37, 13, 17, 3,
                0, 38, 0, 38, 24, 3, 25, 0, 8, 0,
                8, 14, 0, 8, 15, 0, 8, 26, 0, 20,
                4, 19, 4, 21, 0
        };

#endif

#if YYDEBUG
/* YYRLINE[YYN] -- source line where rule number YYN was defined. */
static const short yyrline[] =
        {
                0, 62, 63, 66, 69, 72, 75, 79, 94, 98,
                105, 114, 123, 131, 142, 164, 210, 214, 218, 225,
                235, 253, 270, 294, 304, 318, 321, 328, 354, 357,
                363, 366, 372, 376, 380, 384, 390
        };
#endif


#if (YYDEBUG) || defined YYERROR_VERBOSE

/* YYTNAME[TOKEN_NUM] -- String name of the token TOKEN_NUM. */
static const char *const yytname[] =
        {
                "$", "error", "$undefined.", "INT", "FLOAT", "COMPLEX", "SIGN",
                "STRING",
                "PARTICLENAME", "DEBUG", "CHANNEL", "MODE", "MASSDEP",
                "HELICITY",
                "'+'", "'-'", "'*'", "'='", "';'", "','", "'('", "')'", "'{'",
                "'}'",
                "'['", "']'", "'0'", "input", "statement", "waveexp", "wave",
                "resonance", "quantum_nums", "quantum_num", "Decay", "particle",
                "pstate", "particleID", "particleCharge", "complex", 0
        };
#endif

/* YYR1[YYN] -- Symbol number of symbol that rule YYN derives. */
static const short yyr1[] =
        {
                0, 27, 27, 28, 28, 28, 28, 28, 29, 29,
                29, 29, 29, 29, 30, 30, 31, 32, 33, 34,
                34, 34, 34, 34, 34, 35, 35, 35, 36, 36,
                37, 37, 38, 38, 38, 38, 39
        };

/* YYR2[YYN] -- Number of symbols composing right hand side of rule YYN. */
static const short yyr2[] =
        {
                0, 0, 2, 4, 4, 4, 6, 2, 1, 3,
                3, 3, 3, 3, 2, 1, 1, 3, 3, 5,
                7, 7, 10, 6, 6, 1, 2, 5, 1, 4,
                1, 4, 1, 2, 2, 2, 5
        };

/* YYDEFACT[S] -- default rule to reduce with in state S when YYTABLE
   doesn't specify something else to do.  Zero means the default is an
   error. */
static const short yydefact[] =
        {
                1, 0, 0, 0, 0, 0, 0, 0, 0, 2,
                0, 8, 0, 16, 0, 15, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 32, 0, 25, 28, 30,
                0, 0, 7, 14, 0, 0, 12, 18, 0, 0,
                0, 0, 0, 0, 9, 33, 34, 35, 0, 26,
                0, 0, 10, 11, 17, 13, 0, 3, 4, 5,
                0, 0, 0, 0, 0, 0, 0, 0, 36, 0,
                19, 0, 0, 0, 29, 31, 6, 23, 0, 0,
                24, 27, 0, 20, 21, 0, 0, 22, 0, 0
        };

static const short yydefgoto[] =
        {
                1, 9, 10, 11, 12, 13, 14, 15, 26, 27,
                28, 29, 16
        };

static const short yypact[] =
        {
                -32768, 5, -5, -9, 32, 39, 40, 13, 43, -32768,
                26, -32768, 8, -32768, 6, -32768, 42, 25, 45, 56,
                35, 53, 34, 44, 22, -8, 43, 8, 49, 41,
                25, 25, -32768, -32768, 6, 25, -32768, -32768, 47, 46,
                50, 52, 59, 64, -32768, -32768, -32768, -32768, 31, 57,
                54, 69, -32768, -32768, -32768, -32768, 43, -32768, -32768,
                -32768,
                55, -2, 58, 70, 60, 71, 61, 62, -32768, 65,
                -32768, 51, 66, 72, -32768, -32768, -32768, -32768, -4, 67,
                -32768, -32768, 68, -32768, -32768, 75, 73, -32768, 81, -32768
        };

static const short yypgoto[] =
        {
                -32768, -32768, -7, -32768, -32768, -32768, -12, 19, -22,
                -32768,
                27, -32768, -32768
        };


#define YYLAST      96


static const short yytable[] =
        {
                24, 69, 34, 82, 48, 88, 45, 46, 18, 2,
                36, 17, 3, 23, 4, 5, 6, 22, 47, 83,
                23, 70, 54, 52, 53, 7, 63, 8, 55, 2,
                8, 33, 23, 7, 61, 8, 30, 31, 62, 25,
                30, 31, 40, 44, 32, 7, 49, 8, 37, 19,
                17, 25, 38, 42, 78, 79, 20, 21, 35, 39,
                41, 43, 50, 60, 57, 51, 56, 37, 58, 64,
                59, 65, 66, 72, 74, 71, 68, 73, 86, 81,
                76, 89, 0, 67, 0, 85, 75, 0, 77, 80,
                84, 0, 0, 0, 0, 0, 87
        };

static const short yycheck[] =
        {
                7, 3, 14, 7, 26, 0, 14, 15, 17, 4,
                17, 16, 7, 7, 9, 10, 11, 4, 26, 23,
                7, 23, 34, 30, 31, 20, 48, 22, 35, 4,
                22, 12, 7, 20, 3, 22, 14, 15, 7, 8,
                14, 15, 7, 21, 18, 20, 27, 22, 3, 17,
                16, 8, 7, 19, 3, 4, 17, 17, 16, 3,
                7, 17, 13, 4, 18, 24, 19, 3, 18, 12,
                18, 17, 3, 3, 3, 17, 21, 17, 3, 7,
                18, 0, -1, 56, -1, 17, 25, -1, 23, 23,
                23, -1, -1, -1, -1, -1, 23
        };
/* -*-C-*-  Note some compilers choke on comments on `#line' lines.  */
#line 3 "/usr/share/bison/bison.simple"

/* Skeleton output parser for bison,

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002 Free Software
   Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330,
   Boston, MA 02111-1307, USA.  */

/* As a special exception, when this file is copied by Bison into a
   Bison output file, you may use that output file without restriction.
   This special exception was added by the Free Software Foundation
   in version 1.24 of Bison.  */

/* This is the parser code that is written into each bison parser when
   the %semantic_parser declaration is not specified in the grammar.
   It was written by Richard Stallman by simplifying the hairy parser
   used when %semantic_parser is specified.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

#if !defined (yyoverflow) || defined (YYERROR_VERBOSE)

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# if YYSTACK_USE_ALLOCA
#  define YYSTACK_ALLOC alloca
# else
#  ifndef YYSTACK_USE_ALLOCA
#   if defined (alloca) || defined (_ALLOCA_H)
#    define YYSTACK_ALLOC alloca
#   else
#    ifdef __GNUC__
#     define YYSTACK_ALLOC __builtin_alloca
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
/* Pacify GCC's `empty if-body' warning. */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
# else
#  if defined (__STDC__) || defined (__cplusplus)
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   define YYSIZE_T size_t
#  endif
#  define YYSTACK_ALLOC malloc
#  define YYSTACK_FREE free
# endif
#endif /* ! defined (yyoverflow) || defined (YYERROR_VERBOSE) */


#if (!defined (yyoverflow) \
 && (!defined (__cplusplus) \
 || (YYLTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  short yyss;
  YYSTYPE yyvs;
# if YYLSP_NEEDED
  YYLTYPE yyls;
# endif
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAX (sizeof (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# if YYLSP_NEEDED
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE) + sizeof (YYLTYPE))  \
      + 2 * YYSTACK_GAP_MAX)
# else
#  define YYSTACK_BYTES(N) \
     ((N) * (sizeof (short) + sizeof (YYSTYPE))             \
      + YYSTACK_GAP_MAX)
# endif

/* Copy COUNT objects from FROM to TO.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if 1 < __GNUC__
#   define YYCOPY(To, From, Count) \
      __builtin_memcpy (To, From, (Count) * sizeof (*(From)))
#  else
#   define YYCOPY(To, From, Count)      \
      do                    \
    {                   \
      YYSIZE_T yyi;        \
      for (yyi = 0; yyi < (Count); yyi++)   \
        (To)[yyi] = (From)[yyi];        \
    }                   \
      while (0)
#  endif
# endif

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack)                    \
    do                                  \
      {                                 \
    YYSIZE_T yynewbytes;                        \
    YYCOPY (&yyptr->Stack, Stack, yysize);              \
    Stack = &yyptr->Stack;                      \
    yynewbytes = yystacksize * sizeof (*Stack) + YYSTACK_GAP_MAX;   \
    yyptr += yynewbytes / sizeof (*yyptr);              \
      }                                 \
    while (0)

#endif


#if !defined (YYSIZE_T) && defined (__SIZE_TYPE__)
# define YYSIZE_T __SIZE_TYPE__
#endif
#if !defined (YYSIZE_T) && defined (size_t)
# define YYSIZE_T size_t
#endif
#if !defined (YYSIZE_T)
# if defined (__STDC__) || defined (__cplusplus)
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# endif
#endif
#if !defined (YYSIZE_T)
# define YYSIZE_T unsigned int
#endif

#define yyerrok     (yyerrstatus = 0)
#define yyclearin   (yychar = YYEMPTY)
#define YYEMPTY     -2
#define YYEOF       0
#define YYACCEPT    goto yyacceptlab
#define YYABORT     goto yyabortlab
#define YYERROR     goto yyerrlab1
/* Like YYERROR except do call yyerror.  This remains here temporarily
   to ease the transition to the new meaning of YYERROR, for GCC.
   Once GCC version 2 has supplanted version 1, this can go.  */
#define YYFAIL      goto yyerrlab
#define YYRECOVERING()  (!!yyerrstatus)
#define YYBACKUP(Token, Value)                  \
do                              \
  if (yychar == YYEMPTY && yylen == 1)              \
    {                               \
      yychar = (Token);                     \
      yylval = (Value);                     \
      yychar1 = YYTRANSLATE (yychar);               \
      YYPOPSTACK;                       \
      goto yybackup;                        \
    }                               \
  else                              \
    {                               \
      yyerror ("syntax error: cannot back up");         \
      YYERROR;                          \
    }                               \
while (0)

#define YYTERROR    1
#define YYERRCODE   256


/* YYLLOC_DEFAULT -- Compute the default location (before the actions
   are run).

   When YYLLOC_DEFAULT is run, CURRENT is set the location of the
   first token.  By default, to implement support for ranges, extend
   its range to the last symbol.  */

#ifndef YYLLOC_DEFAULT
# define YYLLOC_DEFAULT(Current, Rhs, N)        \
   Current.last_line   = Rhs[N].last_line;  \
   Current.last_column = Rhs[N].last_column;
#endif


/* YYLEX -- calling `yylex' with the right arguments.  */

#if YYPURE
# if YYLSP_NEEDED
#  ifdef YYLEX_PARAM
#   define YYLEX        yylex (&yylval, &yylloc, YYLEX_PARAM)
#  else
#   define YYLEX        yylex (&yylval, &yylloc)
#  endif
# else /* !YYLSP_NEEDED */
#  ifdef YYLEX_PARAM
#   define YYLEX        yylex (&yylval, YYLEX_PARAM)
#  else
#   define YYLEX        yylex (&yylval)
#  endif
# endif /* !YYLSP_NEEDED */
#else /* !YYPURE */
# define YYLEX          yylex ()
#endif /* !YYPURE */


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF

#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */

#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)            \
do {                        \
  if (yydebug)                  \
    YYFPRINTF Args;             \
} while (0)
/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args)
#endif /* !YYDEBUG */

/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   SIZE_MAX < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#if YYMAXDEPTH == 0
# undef YYMAXDEPTH
#endif

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif

#ifdef YYERROR_VERBOSE

# ifndef yystrlen
#  if defined (__GLIBC__) && defined (_STRING_H)
#   define yystrlen strlen
#  else
/* Return the length of YYSTR.  */
static YYSIZE_T
#   if defined (__STDC__) || defined (__cplusplus)
yystrlen (const char *yystr)
#   else
yystrlen (yystr)
     const char *yystr;
#   endif
{
  const char *yys = yystr;

  while (*yys++ != '\0')
    continue;

  return yys - yystr - 1;
}
#  endif
# endif

# ifndef yystpcpy
#  if defined (__GLIBC__) && defined (_STRING_H) && defined (_GNU_SOURCE)
#   define yystpcpy stpcpy
#  else
/* Copy YYSRC to YYDEST, returning the address of the terminating '\0' in
   YYDEST.  */
static char *
#   if defined (__STDC__) || defined (__cplusplus)
yystpcpy (char *yydest, const char *yysrc)
#   else
yystpcpy (yydest, yysrc)
     char *yydest;
     const char *yysrc;
#   endif
{
   char *yyd = yydest;
   const char *yys = yysrc;

  while ((*yyd++ = *yys++) != '\0')
    continue;

  return yyd - 1;
}
#  endif
# endif
#endif

#line 315 "/usr/share/bison/bison.simple"


/* The user can define YYPARSE_PARAM as the name of an argument to be passed
   into yyparse.  The argument should have type void *.
   It should actually point to an object.
   Grammar actions can access the variable by casting it
   to the proper pointer type.  */

#ifdef YYPARSE_PARAM
# if defined (__STDC__) || defined (__cplusplus)
#  define YYPARSE_PARAM_ARG void *YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL
# else
#  define YYPARSE_PARAM_ARG YYPARSE_PARAM
#  define YYPARSE_PARAM_DECL void *YYPARSE_PARAM;
# endif
#else /* !YYPARSE_PARAM */
# define YYPARSE_PARAM_ARG
# define YYPARSE_PARAM_DECL
#endif /* !YYPARSE_PARAM */

/* Prevent warning if -Wstrict-prototypes.  */
#ifdef __GNUC__
# ifdef YYPARSE_PARAM
int yyparse (void *);
# else

int yyparse(void);

# endif
#endif

/* YY_DECL_VARIABLES -- depending whether we use a pure parser,
   variables are global, or local to YYPARSE.  */

#define YY_DECL_NON_LSP_VARIABLES           \
/* The lookahead symbol.  */                \
int yychar;                     \
                            \
/* The semantic value of the lookahead symbol. */   \
YYSTYPE yylval;                     \
                            \
/* Number of parse errors so far.  */           \
int yynerrs;

#if YYLSP_NEEDED
# define YY_DECL_VARIABLES          \
YY_DECL_NON_LSP_VARIABLES           \
                        \
/* Location data for the lookahead symbol.  */  \
YYLTYPE yylloc;
#else
# define YY_DECL_VARIABLES          \
YY_DECL_NON_LSP_VARIABLES
#endif


/* If nonreentrant, generate the variables here. */

#if !YYPURE
YY_DECL_VARIABLES
#endif  /* !YYPURE */

int
yyparse(YYPARSE_PARAM_ARG)
YYPARSE_PARAM_DECL
{
    /* If reentrant, generate the variables here. */
#if YYPURE
    YY_DECL_VARIABLES
#endif  /* !YYPURE */

    int yystate;
    int yyn;
    int yyresult;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus;
    /* Lookahead token as an internal (translated) token number.  */
    int yychar1 = 0;

    /* Three stacks and their tools:
       `yyss': related to states,
       `yyvs': related to semantic values,
       `yyls': related to locations.

       Refer to the stacks thru separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* The state stack. */
    short yyssa[YYINITDEPTH];
    short *yyss = yyssa;
    short *yyssp;

    /* The semantic value stack.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp;

#if YYLSP_NEEDED
    /* The location stack.  */
    YYLTYPE yylsa[YYINITDEPTH];
    YYLTYPE *yyls = yylsa;
    YYLTYPE *yylsp;
#endif

#if YYLSP_NEEDED
# define YYPOPSTACK   (yyvsp--, yyssp--, yylsp--)
#else
# define YYPOPSTACK   (yyvsp--, yyssp--)
#endif

    YYSIZE_T yystacksize = YYINITDEPTH;


    /* The variables used to return semantic value and location from the
       action routines.  */
    YYSTYPE yyval;
#if YYLSP_NEEDED
    YYLTYPE yyloc;
#endif

    /* When reducing, the number of symbols on the RHS of the reduced
       rule. */
    int yylen;

    YYDPRINTF ((stderr, "Starting parse\n"));

    yystate = 0;
    yyerrstatus = 0;
    yynerrs = 0;
    yychar = YYEMPTY;     /* Cause a token to be read.  */

    /* Initialize stack pointers.
       Waste one element of value and location stack
       so that they stay on the same level as the state stack.
       The wasted elements are never initialized.  */

    yyssp = yyss;
    yyvsp = yyvs;
#if YYLSP_NEEDED
    yylsp = yyls;
#endif
    goto yysetstate;

/*------------------------------------------------------------.
| yynewstate -- Push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
    yynewstate:
    /* In all cases, when you get here, the value and location stacks
       have just been pushed. so pushing a state here evens the stacks.
       */
    yyssp++;

    yysetstate:
    *yyssp = yystate;

    if (yyssp >= yyss + yystacksize - 1) {
        /* Get the current used size of the three stacks, in elements.  */
        YYSIZE_T yysize = yyssp - yyss + 1;

#ifdef yyoverflow
        {
      /* Give user a chance to reallocate the stack. Use copies of
         these so that the &'s don't force the real ones into
         memory.  */
      YYSTYPE *yyvs1 = yyvs;
      short *yyss1 = yyss;

      /* Each stack pointer address is followed by the size of the
         data in use in that stack, in bytes.  */
# if YYLSP_NEEDED
      YYLTYPE *yyls1 = yyls;
      /* This used to be a conditional around just the two extra args,
         but that might be undefined if yyoverflow is a macro.  */
      yyoverflow ("parser stack overflow",
              &yyss1, yysize * sizeof (*yyssp),
              &yyvs1, yysize * sizeof (*yyvsp),
              &yyls1, yysize * sizeof (*yylsp),
              &yystacksize);
      yyls = yyls1;
# else
      yyoverflow ("parser stack overflow",
              &yyss1, yysize * sizeof (*yyssp),
              &yyvs1, yysize * sizeof (*yyvsp),
              &yystacksize);
# endif
      yyss = yyss1;
      yyvs = yyvs1;
        }
#else /* no yyoverflow */
# ifndef YYSTACK_RELOCATE
        goto yyoverflowlab;
# else
        /* Extend the stack our own way.  */
        if (yystacksize >= YYMAXDEPTH)
      goto yyoverflowlab;
        yystacksize *= 2;
        if (yystacksize > YYMAXDEPTH)
      yystacksize = YYMAXDEPTH;

        {
      short *yyss1 = yyss;
      union yyalloc *yyptr =
        (union yyalloc *) YYSTACK_ALLOC (YYSTACK_BYTES (yystacksize));
      if (! yyptr)
        goto yyoverflowlab;
      YYSTACK_RELOCATE (yyss);
      YYSTACK_RELOCATE (yyvs);
# if YYLSP_NEEDED
      YYSTACK_RELOCATE (yyls);
# endif
# undef YYSTACK_RELOCATE
      if (yyss1 != yyssa)
        YYSTACK_FREE (yyss1);
        }
# endif
#endif /* no yyoverflow */

        yyssp = yyss + yysize - 1;
        yyvsp = yyvs + yysize - 1;
#if YYLSP_NEEDED
        yylsp = yyls + yysize - 1;
#endif

        YYDPRINTF ((stderr, "Stack size increased to %lu\n",
                (unsigned long int) yystacksize));

        if (yyssp >= yyss + yystacksize - 1)
            YYABORT;
    }

    YYDPRINTF ((stderr, "Entering state %d\n", yystate));

    goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
    yybackup:

/* Do appropriate processing given the current state.  */
/* Read a lookahead token if we need one and don't already have one.  */
/* yyresume: */

    /* First try to decide what to do without reference to lookahead token.  */

    yyn = yypact[yystate];
    if (yyn == YYFLAG)
        goto yydefault;

    /* Not known => get a lookahead token if don't already have one.  */

    /* yychar is either YYEMPTY or YYEOF
       or a valid token in external form.  */

    if (yychar == YYEMPTY) {
        YYDPRINTF ((stderr, "Reading a token: "));
        yychar = YYLEX;
    }

    /* Convert token to internal form (in yychar1) for indexing tables with */

    if (yychar <= 0)      /* This means end of input. */
    {
        yychar1 = 0;
        yychar = YYEOF;       /* Don't call YYLEX any more */

        YYDPRINTF ((stderr, "Now at end of input.\n"));
    } else {
        yychar1 = YYTRANSLATE (yychar);

#if YYDEBUG
        /* We have to keep this `#if YYDEBUG', since we use variables
       which are defined only if `YYDEBUG' is set.  */
        if (yydebug) {
            YYFPRINTF(stderr, "Next token is %d (%s",
                      yychar, yytname[yychar1]);
            /* Give the individual parser a way to print the precise
               meaning of a token, for further debugging info.  */
# ifdef YYPRINT
            YYPRINT (stderr, yychar, yylval);
# endif
            YYFPRINTF(stderr, ")\n");
        }
#endif
    }

    yyn += yychar1;
    if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != yychar1)
        goto yydefault;

    yyn = yytable[yyn];

    /* yyn is what to do for this token type in this state.
       Negative => reduce, -yyn is rule number.
       Positive => shift, yyn is new state.
         New state is final state => don't bother to shift,
         just return success.
       0, or most negative number => error.  */

    if (yyn < 0) {
        if (yyn == YYFLAG)
            goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
    } else if (yyn == 0)
        goto yyerrlab;

    if (yyn == YYFINAL)
        YYACCEPT;

    /* Shift the lookahead token.  */
    YYDPRINTF ((stderr, "Shifting token %d (%s), ",
            yychar, yytname[yychar1]));

    /* Discard the token being shifted unless it is eof.  */
    if (yychar != YYEOF)
        yychar = YYEMPTY;

    *++yyvsp = yylval;
#if YYLSP_NEEDED
    *++yylsp = yylloc;
#endif

    /* Count tokens shifted since error; after three, turn off error
       status.  */
    if (yyerrstatus)
        yyerrstatus--;

    yystate = yyn;
    goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
    yydefault:
    yyn = yydefact[yystate];
    if (yyn == 0)
        goto yyerrlab;
    goto yyreduce;


/*-----------------------------.
| yyreduce -- Do a reduction.  |
`-----------------------------*/
    yyreduce:
    /* yyn is the number of a rule to reduce with.  */
    yylen = yyr2[yyn];

    /* If YYLEN is nonzero, implement the default value of the action:
       `$$ = $1'.

       Otherwise, the following line sets YYVAL to the semantic value of
       the lookahead token.  This behavior is undocumented and Bison
       users should not rely upon it.  Assigning to YYVAL
       unconditionally makes the parser a bit smaller, and it avoids a
       GCC warning that YYVAL may be used uninitialized.  */
    yyval = yyvsp[1 - yylen];

#if YYLSP_NEEDED
    /* Similarly for the default location.  Let the user run additional
       commands if for instance locations are ranges.  */
    yyloc = yylsp[1-yylen];
    YYLLOC_DEFAULT (yyloc, (yylsp - yylen), yylen);
#endif

#if YYDEBUG
    /* We have to keep this `#if YYDEBUG', since we use variables which
       are defined only if `YYDEBUG' is set.  */
    if (yydebug) {
        int yyi;

        YYFPRINTF(stderr, "Reducing via rule %d (line %d), ",
                  yyn, yyrline[yyn]);

        /* Print the symbols being reduced, and their result.  */
        for (yyi = yyprhs[yyn]; yyrhs[yyi] > 0; yyi++)
            YYFPRINTF(stderr, "%s ", yytname[yyrhs[yyi]]);
        YYFPRINTF(stderr, " -> %s\n", yytname[yyr1[yyn]]);
    }
#endif

    switch (yyn) {

        case 3:
#line 66 "keyParse.yy"
        {
            debug = yyvsp[-1].num;
        }
            break;
        case 4:
#line 69 "keyParse.yy"
        {
            wv.channel(yyvsp[-1].string);
        }
            break;
        case 5:
#line 72 "keyParse.yy"
        {
            mode = yyvsp[-1].string;
        }
            break;
        case 6:
#line 75 "keyParse.yy"
        {
            t_part_init = yyvsp[-3].string;
            t_part_final = yyvsp[-1].Particle;
        }
            break;
        case 7:
#line 79 "keyParse.yy"
        {
            if (mode == "binary") {
                cout.write((char *) yyvsp[-1].Cnum, sizeof(complex<double>));
            } else {
                cout << "Mass = " << ~(wv.get4P()) << "\t";
                if (wv.channel() == "t") {
                    cout << "t = " << (e.beam().get4P() - wv.get4P()).lenSq()
                         << "\t";
                }
                cout << "Amp = " << *yyvsp[-1].Cnum << endl;
            }
            delete yyvsp[-1].Cnum;
        }
            break;
        case 8:
#line 94 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(*yyvsp[0].Cnum);
            delete (yyvsp[0].Cnum);
        }
            break;
        case 9:
#line 98 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(*yyvsp[-1].Cnum);
            if (debug) {
                cout << " ( " << *yyvsp[-1].Cnum << " ) = " << *yyval.Cnum
                     << endl;
            }
            delete (yyvsp[-1].Cnum);
        }
            break;
        case 10:
#line 106 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum + *yyvsp[0].Cnum);
            if (debug) {
                cout << *yyvsp[-2].Cnum << " + " << *yyvsp[0].Cnum << " = "
                     << *yyval.Cnum << endl;
            }
            delete (yyvsp[-2].Cnum);
            delete (yyvsp[0].Cnum);
        }
            break;
        case 11:
#line 115 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum - *yyvsp[0].Cnum);
            if (debug) {
                cout << *yyvsp[-2].Cnum << " - " << *yyvsp[0].Cnum << " = "
                     << *yyval.Cnum << endl;
            }
            delete (yyvsp[-2].Cnum);
            delete (yyvsp[0].Cnum);
        }
            break;
        case 12:
#line 124 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(yyvsp[-2].Fnum * *yyvsp[0].Cnum);
            if (debug) {
                cout << yyvsp[-2].Fnum << " * " << *yyvsp[0].Cnum << " = "
                     << *yyval.Cnum << endl;
            }
            delete (yyvsp[0].Cnum);
        }
            break;
        case 13:
#line 132 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(*yyvsp[-2].Cnum * *yyvsp[0].Cnum);
            if (debug) {
                cout << *yyvsp[-2].Cnum << " * " << *yyvsp[0].Cnum << " = "
                     << *yyval.Cnum << endl;
            }
            delete (yyvsp[-2].Cnum);
            delete (yyvsp[0].Cnum);
        }
            break;
        case 14:
#line 142 "keyParse.yy"
        {
            wv.setDecay(*yyvsp[0].Decay);
            delete yyvsp[0].Decay;
            if (debug) {
                cout << "@@Found a wave" << endl;
                wv.print();
                cout << "@@Filling wave" << endl;
            }
            wv.fill(e, debug);
            if (debug) {
                cout << "@@Wave before boosts" << endl;
                wv.print();
            }
            wv.setupFrames(debug);
            if (debug) {
                cout << "@@Wave after boosts" << endl;
                wv.print();
            }
            amp = wv.decayAmp(debug);
            yyval.Cnum = new complex<double>(amp);
            nwave++;
        }
            break;
        case 15:
#line 164 "keyParse.yy"
        {
            wv.setDecay(*yyvsp[0].Decay);
            delete yyvsp[0].Decay;
            if (debug) {
                cout << "@@Found a wave" << endl;
                wv.print();
                cout << "@@Filling wave" << endl;
            }
            wv.fill(e, debug);
            if (debug) {
                cout << "@@Wave before boosts" << endl;
                wv.print();
            }
            wv.setupFrames(debug);
            if (debug) {
                cout << "@@Wave after boosts" << endl;
                wv.print();
            }
            if (debug) {
                cout << "This should compute Decay amplitude expt wave" << endl;
            }
            double t = 0.0;
            fourVec t_init(0.0, threeVec(0.0, 0.0, 0.0));
            if (t_part_init == "beam") {
                t_init = wv.getBeam();
            } else if (t_part_init == "target") {
                t_init = wv.getTarget();
            } else {
                cerr << "unknown initial t specifier: " << t_part_init << endl;
                abort();
            }
            t = (t_init - *wv.get4P(t_part_final, debug)).lenSq();
            if (debug) {
                cout << "calulating amplitude with t = " << t << endl;
            }
            delete t_part_final;

            wv.setT(t);
            amp = wv.decayAmp(debug);
            yyval.Cnum = new complex<double>(amp);
            nwave++;
        }
            break;
        case 16:
#line 210 "keyParse.yy"
        {
        }
            break;
        case 17:
#line 214 "keyParse.yy"
        {
        }
            break;
        case 18:
#line 218 "keyParse.yy"
        {
            if (!strcmp(yyvsp[-2].string, "J")) wv.setJ(yyvsp[0].num);
            if (!strcmp(yyvsp[-2].string, "M")) wv.setM(yyvsp[0].num);
            if (!strcmp(yyvsp[-2].string, "P")) wv.setP(yyvsp[0].num);
        }
            break;
        case 19:
#line 225 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-3].Particle);
            d->addChild(*yyvsp[-2].Particle);
            delete yyvsp[-3].Particle;
            delete yyvsp[-2].Particle;
            d->setL(yyvsp[-1].num);
            d->calculateS();
            yyval.Decay = d;
        }
            break;
        case 20:
#line 235 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-5].Particle);
            d->addChild(*yyvsp[-4].Particle);
            delete yyvsp[-5].Particle;
            delete yyvsp[-4].Particle;
            if (!strcmp(yyvsp[-3].string, "l")) {
                d->setL(yyvsp[-1].num);
                d->calculateS();
            } else {
                cerr << "unexpected field at line " << lineno << endl;
                cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                cerr << "expected \'l\'" << endl;
                exit(1);
            }
            yyval.Decay = d;
        }
            break;
        case 21:
#line 253 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-5].Particle);
            d->addChild(*yyvsp[-4].Particle);
            delete yyvsp[-5].Particle;
            delete yyvsp[-4].Particle;
            if (!strcmp(yyvsp[-3].string, "b")) {
                wv.setSlope(yyvsp[-1].Fnum);
            } else {
                cerr << "unexpected field at line " << lineno << endl;
                cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                cerr << "expected \'b\'" << endl;
                exit(1);
            }
            yyval.Decay = d;
        }
            break;
        case 22:
#line 270 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-8].Particle);
            d->addChild(*yyvsp[-7].Particle);
            delete yyvsp[-8].Particle;
            delete yyvsp[-7].Particle;
            if (!strcmp(yyvsp[-6].string, "l")) {
                d->setL(yyvsp[-4].num);
            } else {
                cerr << "expecting \'l\' at line " << lineno << endl;
                cerr << "found \'" << yyvsp[-6].string << "\'" << endl;
                exit(1);
            }
            if (!strcmp(yyvsp[-3].string, "s")) {
                d->setS(yyvsp[-1].num);
            } else {
                cerr << "expecting \'l\' at line " << lineno << endl;
                cerr << "found \'" << yyvsp[-3].string << "\'" << endl;
                exit(1);
            }
            yyval.Decay = d;
        }
            break;
        case 23:
#line 294 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-4].Particle);
            d->addChild(*yyvsp[-3].Particle);
            delete yyvsp[-4].Particle;
            delete yyvsp[-3].Particle;
            d->setL(yyvsp[-2].num);
            d->setS(yyvsp[-1].num);
            yyval.Decay = d;
        }
            break;
        case 24:
#line 304 "keyParse.yy"
        {
            Decay * d = new Decay;
            d->addChild(*yyvsp[-4].Particle);
            d->addChild(*yyvsp[-3].Particle);
            d->addChild(*yyvsp[-2].Particle);
            delete yyvsp[-4].Particle;
            delete yyvsp[-3].Particle;
            delete yyvsp[-2].Particle;
            d->setL(yyvsp[-1].num);
            d->calculateS();
            yyval.Decay = d;
        }
            break;
        case 25:
#line 318 "keyParse.yy"
        {
            yyval.Particle = yyvsp[0].Particle;
        }
            break;
        case 26:
#line 321 "keyParse.yy"
        {
            yyvsp[-1].Particle->setDecay(*yyvsp[0].Decay);
            massDep *bw = new breitWigner();
            yyvsp[-1].Particle->setMassDep(bw);
            delete yyvsp[0].Decay;
            yyval.Particle = yyvsp[-1].Particle;
        }
            break;
        case 27:
#line 328 "keyParse.yy"
        {
            yyvsp[-4].Particle->setDecay(*yyvsp[-3].Decay);
            massDep *md;
            if (!strcmp(yyvsp[0].string, "bw")) {
                md = new breitWigner();
            } else if (!strcmp(yyvsp[0].string, "amp")) {
                md = new AMP_M();
            } else if (!strcmp(yyvsp[0].string, "amp_ves")) {
                md = new AMP_ves();
            } else if (!strcmp(yyvsp[0].string, "flat")) {
                md = new flat();
            } else {
                cerr << "unknown mass dependence: " << yyvsp[0].string;
                cerr << " at line " << lineno << endl;
                exit(1);
            }
            yyvsp[-4].Particle->setMassDep(md);
            delete yyvsp[-3].Decay;
            yyval.Particle = yyvsp[-4].Particle;
        }
            break;
        case 28:
#line 354 "keyParse.yy"
        {
            yyval.Particle = yyvsp[0].Particle;
        }
            break;
        case 29:
#line 357 "keyParse.yy"
        {
            yyvsp[-3].Particle->addHelicity(yyvsp[0].num);
            yyval.Particle = yyvsp[-3].Particle;
        }
            break;
        case 30:
#line 363 "keyParse.yy"
        {
            yyval.Particle = yyvsp[0].Particle;
        }
            break;
        case 31:
#line 366 "keyParse.yy"
        {
            yyvsp[-3].Particle->Index(yyvsp[-1].num);
            yyval.Particle = yyvsp[-3].Particle;
        }
            break;
        case 32:
#line 372 "keyParse.yy"
        {
            particle *p = new particle(PDGtable.get(yyvsp[0].string), 0);
            yyval.Particle = p;
        }
            break;
        case 33:
#line 376 "keyParse.yy"
        {
            particle *p = new particle(PDGtable.get(yyvsp[-1].string), +1);
            yyval.Particle = p;
        }
            break;
        case 34:
#line 380 "keyParse.yy"
        {
            particle *p = new particle(PDGtable.get(yyvsp[-1].string), -1);
            yyval.Particle = p;
        }
            break;
        case 35:
#line 384 "keyParse.yy"
        {
            particle *p = new particle(PDGtable.get(yyvsp[-1].string), 0);
            yyval.Particle = p;
        }
            break;
        case 36:
#line 390 "keyParse.yy"
        {
            yyval.Cnum = new complex<double>(yyvsp[-3].Fnum, yyvsp[-1].Fnum);
        }
            break;
    }

#line 705 "/usr/share/bison/bison.simple"


    yyvsp -= yylen;
    yyssp -= yylen;
#if YYLSP_NEEDED
    yylsp -= yylen;
#endif

#if YYDEBUG
    if (yydebug) {
        short *yyssp1 = yyss - 1;
        YYFPRINTF(stderr, "state stack now");
        while (yyssp1 != yyssp)
            YYFPRINTF(stderr, " %d", *++yyssp1);
        YYFPRINTF(stderr, "\n");
    }
#endif

    *++yyvsp = yyval;
#if YYLSP_NEEDED
    *++yylsp = yyloc;
#endif

    /* Now `shift' the result of the reduction.  Determine what state
       that goes to, based on the state we popped back to and the rule
       number reduced by.  */

    yyn = yyr1[yyn];

    yystate = yypgoto[yyn - YYNTBASE] + *yyssp;
    if (yystate >= 0 && yystate <= YYLAST && yycheck[yystate] == *yyssp)
        yystate = yytable[yystate];
    else
        yystate = yydefgoto[yyn - YYNTBASE];

    goto yynewstate;


/*------------------------------------.
| yyerrlab -- here on detecting error |
`------------------------------------*/
    yyerrlab:
    /* If not already recovering from an error, report this error.  */
    if (!yyerrstatus) {
        ++yynerrs;

#ifdef YYERROR_VERBOSE
        yyn = yypact[yystate];

        if (yyn > YYFLAG && yyn < YYLAST)
      {
        YYSIZE_T yysize = 0;
        char *yymsg;
        int yyx, yycount;

        yycount = 0;
        /* Start YYX at -YYN if negative to avoid negative indexes in
           YYCHECK.  */
        for (yyx = yyn < 0 ? -yyn : 0;
             yyx < (int) (sizeof (yytname) / sizeof (char *)); yyx++)
          if (yycheck[yyx + yyn] == yyx)
            yysize += yystrlen (yytname[yyx]) + 15, yycount++;
        yysize += yystrlen ("parse error, unexpected ") + 1;
        yysize += yystrlen (yytname[YYTRANSLATE (yychar)]);
        yymsg = (char *) YYSTACK_ALLOC (yysize);
        if (yymsg != 0)
          {
            char *yyp = yystpcpy (yymsg, "parse error, unexpected ");
            yyp = yystpcpy (yyp, yytname[YYTRANSLATE (yychar)]);

            if (yycount < 5)
          {
            yycount = 0;
            for (yyx = yyn < 0 ? -yyn : 0;
                 yyx < (int) (sizeof (yytname) / sizeof (char *));
                 yyx++)
              if (yycheck[yyx + yyn] == yyx)
                {
              const char *yyq = ! yycount ? ", expecting " : " or ";
              yyp = yystpcpy (yyp, yyq);
              yyp = yystpcpy (yyp, yytname[yyx]);
              yycount++;
                }
          }
            yyerror (yymsg);
            YYSTACK_FREE (yymsg);
          }
        else
          yyerror ("parse error; also virtual memory exhausted");
      }
        else
#endif /* defined (YYERROR_VERBOSE) */
        yyerror("parse error");
    }
    goto yyerrlab1;


/*--------------------------------------------------.
| yyerrlab1 -- error raised explicitly by an action |
`--------------------------------------------------*/
    yyerrlab1:
    if (yyerrstatus == 3) {
        /* If just tried and failed to reuse lookahead token after an
       error, discard it.  */

        /* return failure if at end of input */
        if (yychar == YYEOF)
            YYABORT;
        YYDPRINTF ((stderr, "Discarding token %d (%s).\n",
                yychar, yytname[yychar1]));
        yychar = YYEMPTY;
    }

    /* Else will try to reuse lookahead token after shifting the error
       token.  */

    yyerrstatus = 3;      /* Each real token shifted decrements this */

    goto yyerrhandle;


/*-------------------------------------------------------------------.
| yyerrdefault -- current state does not do anything special for the |
| error token.                                                       |
`-------------------------------------------------------------------*/
    yyerrdefault:
#if 0
    /* This is wrong; only states that explicitly want error tokens
       should shift them.  */

    /* If its default is to accept any token, ok.  Otherwise pop it.  */
    yyn = yydefact[yystate];
    if (yyn)
      goto yydefault;
#endif


/*---------------------------------------------------------------.
| yyerrpop -- pop the current state because it cannot handle the |
| error token                                                    |
`---------------------------------------------------------------*/
    yyerrpop:
    if (yyssp == yyss)
        YYABORT;
    yyvsp--;
    yystate = *--yyssp;
#if YYLSP_NEEDED
    yylsp--;
#endif

#if YYDEBUG
    if (yydebug) {
        short *yyssp1 = yyss - 1;
        YYFPRINTF(stderr, "Error: state stack now");
        while (yyssp1 != yyssp)
            YYFPRINTF(stderr, " %d", *++yyssp1);
        YYFPRINTF(stderr, "\n");
    }
#endif

/*--------------.
| yyerrhandle.  |
`--------------*/
    yyerrhandle:
    yyn = yypact[yystate];
    if (yyn == YYFLAG)
        goto yyerrdefault;

    yyn += YYTERROR;
    if (yyn < 0 || yyn > YYLAST || yycheck[yyn] != YYTERROR)
        goto yyerrdefault;

    yyn = yytable[yyn];
    if (yyn < 0) {
        if (yyn == YYFLAG)
            goto yyerrpop;
        yyn = -yyn;
        goto yyreduce;
    } else if (yyn == 0)
        goto yyerrpop;

    if (yyn == YYFINAL)
        YYACCEPT;

    YYDPRINTF ((stderr, "Shifting error token, "));

    *++yyvsp = yylval;
#if YYLSP_NEEDED
    *++yylsp = yylloc;
#endif

    yystate = yyn;
    goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
    yyacceptlab:
    yyresult = 0;
    goto yyreturn;

/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
    yyabortlab:
    yyresult = 1;
    goto yyreturn;

/*---------------------------------------------.
| yyoverflowab -- parser overflow comes here.  |
`---------------------------------------------*/
    yyoverflowlab:
    yyerror("parser stack overflow");
    yyresult = 2;
    /* Fall through.  */

    yyreturn:
#ifndef yyoverflow
    if (yyss != yyssa)
        YYSTACK_FREE (yyss);
#endif
    return yyresult;
}

#line 394 "keyParse.yy"


