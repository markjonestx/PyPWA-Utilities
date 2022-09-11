/*
 * ntypes.h
 *
*/

#pragma once

/*
static const char sccsid_ntypesH[] = "@(#)ntypes.h\t4.1\tCreated 10 Dec 1995 18:48:47, \tcompiled "__DATE__;
*/

/*
 * This include file defines the basic types
 * used to define the data structures written on E852 data tapes
 * and elsewhere.
 *
 * int8, uint8    -  8 bit wide MSB integers (more commonly known as 'char')
 * int16, uint16  - 16 bit wide MSB integers (usually 'short')
 * int32, uint32  - 32 bit wide MSB integers (usually 'int')
 * int64, uint64  - 64 bit wide MSB integers (on some systems 'long', on others 'long long')
 *
 * float32   - 32 bit wide IEEE-754 floating point numbers (usually C 'float',  FORTRAN 'REAL')
 * float64   - 64 bit wide IEEE-754 floating point numbers (usually C 'double', FORTRAN 'DOUBLE PRECISION')
 * 
 * MSB means "big-endian", like SGI-MIPS, RS6000 and others.
 * LSB means "little-endian", like Intel x86, WindowsNT and others.
 *
 * Read 'man vector' for more details on IEEE-754 floating point numbers.
 *
*/

#ifdef sgi

/*
 * SGI-only definitions 
*/

typedef int                   int32;
typedef unsigned int         uint32;
typedef long long             int64;
typedef unsigned long long   uint64;

#else

/*
 * these definitions work for rs6000 under AIX 3.2.5, untested for others
 * The IBM xlc compiler seem to understand 'long long', but it is not documented, so...
*/

typedef int int32;
typedef unsigned int uint32;
typedef long long int64;
typedef unsigned long long uint64;

#endif

typedef float float32;

/* These are the 'float' structures ... */

typedef struct {
    float32 x, y, z;
} vector3_t;


/* ... and their 'double' siblings */



/*
 * Kludge to get rid of the Particle_t without breaking everybody's code
*/
