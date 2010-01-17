
/*
 * File: ERedgelist.h
 * 
 * 
 */


 
 
/* Notes:
(1) Was getting integer overflow in calculation of w
for values of p on the order of 1e-9, so changed w
to type long long, which seemed to fix the problem
*/

#ifndef EREDGELIST_H
#define EREDGELIST_H

/* Function to produce an ER graph in edgelist format */
SEXP ERedgelist (SEXP N, SEXP p); 

#endif
