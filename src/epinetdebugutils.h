
/*
 * File epinetdebugutils.h
 *
 * Part of the `epinet' R package
 *
 * Debugging utilities for epinet package
 */ 

#ifndef EPINETDEBUGUTILS_H
#define EPINETDEBUGUTILS_H

#include "NetworkFunctions.h" 
#include "NetworkEpiMCMCSEIR.h"

/* Functions to calculate various "pressures" from scratch */
double CalcAii(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *rectimes);
double CalcA(Network *nwp, double *exptimes, double *inftimes, double *rectimes);
double CalcB(double *exptimes, double *inftimes, int m);
double CalcBln(double *exptimes, double *inftimes, int m);
double CalcC(double *inftimes, double *rectimes, int m);
double CalcCln(double *inftimes, double *rectimes, int m);

/* Functions to count number of edges in the graph*/
void TreeCount(TreeNode *edges, int x, int *countptr);
int CountEdges(Network *nwp);

#endif
