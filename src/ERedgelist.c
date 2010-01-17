#include <R.h> 
#include <Rinternals.h> 
#include <R_ext/Rdynload.h> 
#include <Rmath.h>
#include <math.h> 
#include <Rdefines.h>

/**************************
ERedgelist

Generates an Erdos-Renyi graph with a given N
and p.  Based on code given in the article:

network: a Package for Managing Relational Data in R
Journal of Statistical Software
Butts, 2008

***************************/

/* Notes:
(1) Was getting integer overflow in calculation of w
for values of p on the order of 1e-9, so changed w
to type long long, which seemed to fix the problem
*/

SEXP ERedgelist (SEXP N, SEXP p) 
  { 
  int v = 1, count = 0; 
  long long w = -1;
  double r; 
  SEXP tailhead;
  PROTECT_INDEX pind;
  
  GetRNGstate();  /* R function enabling uniform RNG */ 
  
  PROTECT(N = coerceVector(N,INTSXP));
  PROTECT(p = coerceVector(p,REALSXP));
  
  PROTECT_WITH_INDEX(tailhead = allocVector(INTSXP, count),&pind); /*Allocate head/tail*/
  while(v < INTEGER(N)[0])
    {
    r = unif_rand();
    w += 1+ (long long)floor(log(1.0 - r) / log(1.0 - REAL(p)[0]));
    while((w >= v) && (v < INTEGER(N)[0]))
      { 
      w -= v; 
      v++; 
      }
    if (v < INTEGER(N)[0])
      {
      count += 2;
      SET_LENGTH(tailhead,count);
      INTEGER(tailhead)[count-1] = v + 1;
      INTEGER(tailhead)[count-2] = w + 1;
      REPROTECT(tailhead,pind);
      } 
    } 
  PutRNGstate(); /* Must be called after GetRNGstate before returning to R */
  UNPROTECT(3);
  return tailhead;
  }
