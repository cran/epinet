
/*
 * File epinetdebugutils.c
 *
 * Part of the `epinet' R package
 *
 * Debugging utilities for epinet package
 */ 

#include <R.h> 
#include <Rinternals.h> 
#include <R_ext/Rdynload.h> 
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <math.h> 
#include <Rdefines.h>
#include "NetworkEpiMCMCSEIR.h"
#include "NetworkFunctions.h"


/**********************
CalcAii

Calculates the infectious pressure associated with one node.
Uses a recursive algorithm.
Just used as a debugging function and is only called by CalcA().
************************/

double CalcAii(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *rectimes)
{
	Vertex priorinf, latterinf;
	
	if (x != 0) 
	{
		priorinf = MININDEX(inftimes,orig,(edges+x)->value);
 		latterinf = MAXINDEX(inftimes,orig,(edges+x)->value);

		return( MAX(MIN(exptimes[latterinf], rectimes[priorinf]) - inftimes[priorinf] , 0 )
			+ CalcAii(edges, orig, (edges+x)->left, exptimes, inftimes, rectimes)
			+ CalcAii(edges, orig, (edges+x)->right, exptimes, inftimes, rectimes) ); 		
	}
	else return (0);
}


/******************************
CalcA

Calculates the total amount of infectious pressure applied over the course 
of the epidemic.  This differs from AdjustA in that this calculation is done 
from scratch, rather than by adjusting the current value.  Uses a recursive 
algorithm, similar to AdjustA.  Just used as a debugging function to verify A is correct.
Only need to calculate infectious pressure coming from infected nodes; no infectious
pressure comes from the susceptibles.
*****************************/

double CalcA(Network *nwp, double *exptimes, double *inftimes, double *rectimes)
{
	double tempA = 0;
	
	for (Vertex i = 1; i<=(nwp->ninfnodes); i++)
		if (((nwp->outedges)+i)->value != 0) tempA += CalcAii(nwp->outedges,i,i,exptimes,inftimes,rectimes);		/* Calculate contribution to A from all edges */ 												

	return(tempA);
}


/************************
CalcB

Calculates the total amount of transition pressure
applied over the course of the epidemic.  Used
only as a debugging check.
*************************/

double CalcB(double *exptimes, double *inftimes, int m)
{
	double tempB = 0;
	
	for (int i=1; i<=m; i++)
		tempB += inftimes[i] - exptimes[i];	
	
return(tempB);	
}


/************************
CalcBln

Calculates the "log sum" (Bln) of the total amount of transition pressure
applied over the course of the epidemic.  Used
only as a debugging check.
*************************/

double CalcBln(double *exptimes, double *inftimes, int m)
{
	double tempBln = 0;
	
	for (int i=1; i<=m; i++)
		tempBln += log(inftimes[i] - exptimes[i]);	
	
return(tempBln);	
}


/************************
CalcC

Calculates the total amount of removal pressure
applied over the course of the epidemic.  Used
only as a debugging check.
*************************/

double CalcC(double *inftimes, double *rectimes, int m)
{
	double tempC = 0;
	
	for (int i=1; i<=m; i++)
		tempC += rectimes[i] - inftimes[i];	
	
return(tempC);	
}


/************************
CalcCln

Calculates the "log sum" (Cln) of the total amount of removal pressure
applied over the course of the epidemic.  Used
only as a debugging check.
*************************/

double CalcCln(double *inftimes, double *rectimes, int m)
{
	double tempCln = 0;
	
	for (int i=1; i<=m; i++)
		tempCln += log(rectimes[i] - inftimes[i]);	
	
return(tempCln);	
}


/*****************
 void TreeCount
 
 Diagnostic routine that counts the edges in the tree rooted
 at edges[x].  Only called by CountEdgesRecursive.
 *****************/

void TreeCount(TreeNode *edges, int x, int *countptr) 
{
	if (x != 0) 
	{
		TreeCount(edges, (edges+x)->left, countptr);
		(*countptr)++;
		TreeCount(edges, (edges+x)->right, countptr);
	}
}


/*****************
 int CountEdges
 
 Diagnostic routine that counts the number of edges
 in a network.  Uses a recursive algorithm.
 *****************/

int CountEdges(Network *nwp) 
{
	int i, count = 0, *countptr = &count;
	
	for(i=1; i<=(nwp->N); i++) 
		if (((nwp->outedges) + i) -> value != 0)
			TreeCount(nwp->outedges, i, countptr);
	
	return(count);
}
