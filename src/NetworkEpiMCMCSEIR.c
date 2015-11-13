
/*
 *  File NetworkEpiMCMCSEIR.c
 *
 *  Part of the `epinet' R package.
 *  This file contains the routines to perform
 *  inference on the parameters in the network and
 *  epidemic models.
 */


#include <R.h> 
#include <Rinternals.h> 
#include <R_ext/Rdynload.h> 
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <math.h> 
#include <Rdefines.h>
#include <time.h>
#include "NetworkEpiMCMCSEIR.h"
#include "epinetdebugutils.h"


/********************
InitializeTransTree

Initializes the transmission tree, given the input exposure, infection, 
and recovery times.  For each infected node, sample uniformly among 
the vertices that could possibly be its parents.  Called only at the beginning
of the MCMC routine.  Also finds the initial exposed.
*******************/

int InitializeTransTree(Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, int m, int N, Vertex *initexp, double *A)
{
	double maxrand, posspar, lowexptime = exptimes[1];
	Vertex i, j, currpar, lowexp = 1, numinits = 0;

	*A = 0;				/* A holds the infectious pressure for the epidemic */
	transtree[0] = 0;				/* 0th element isn't used - just set equal to 0 */
	for (i=1;i<=m;i++)
	{
		/* Find parent for (infected) Vertex i */
		maxrand = 0; posspar = -1; currpar = -999;
		for (j=1;j < m+1;j++)
			if ((inftimes[j] < exptimes[i]) && (exptimes[i] < rectimes[j]))	/* In order for j to possibly be the parent of i, i must've been infected during j's infectious period */
			{
				posspar = unif_rand();
				if (posspar > maxrand)				
				{
					currpar = j;
					maxrand = posspar;
				}
			}
		transtree[i] = currpar;
		
		if (currpar == -999) numinits++;	/* numinits counts the number of initial infecteds -- should only be one */
		
		/* Add infectious pressure associated with this infection to A 
		   Since we're starting out assuming the only edges are those  
		   required due to the tree, this is the only source of infectious pressure */
		if (transtree[i] != -999) *A += exptimes[i] - inftimes[transtree[i]]; 		
	}
	
	/* Find the initial exposed */
	for (i = 2; i<=m; i++)
		if (exptimes[i] < lowexptime)
		{
			lowexptime = exptimes[i];
			lowexp = i;
		}
	*initexp = lowexp;	
	
	for (i = (m+1); i <= N; i++)	/* Set transmission tree entries for susceptibles to -77; this is the code for non-infected (i.e., susceptible) */
		transtree[i] = -77;
	
	return(numinits);
}


/**************************
 InitializeNetworkFromTree

Initializes the contact network, given the transmission tree.
Just adds the edges associated with the transmission tree,
by calling NetworkInitialize(). 
**************************/

Network InitializeNetworkFromTree(Vertex *transtree, int m, Vertex N)
{
	Vertex head[m-1], tail[m-1];
	int i, count = 0;
	
	for(i=1; i<=m; i++)				/* For each node that isn't the initial exposed, add an edge corresponding to its exposure */
		if (transtree[i] != -999)
		{
			head[count] = i;
			tail[count] = transtree[i];
			count++;
		}
	
	return(NetworkInitialize(head,tail,m-1,m,N));	/* Initialize the network with only the edges corresponding to the trans tree */
}


/******************
ProposedInitExp

Generates a proposal for the identity of the initial exposed node (Kappa).
The proposal is generated uniformly randomly from all of the children of the current 
initial exposed.  This proposal is used in a Hastings update step. This uses a similar 
algorithm as the update for the transmission tree (see DrawTransTree)
******************/

Vertex ProposedInitExp(Vertex initexp, Vertex *transtree, int m)
{
	double maxrand = 0, posspar;
	Vertex propinitexp = initexp;
	
	for (Vertex i=1; i<=m; i++)			
		if (transtree[i] == initexp)
		{
			posspar = unif_rand();
			if (posspar > maxrand)
			{
				maxrand = posspar;
				propinitexp = i;
			}
		}	
	return(propinitexp);	
}


/************************
ProposedInitialExptime

Generates a proposal for the exposure time of the initial exposed (E_Kappa)
The proposal is an exponential RV subtracted from the minimum possible exposure time
Theta is a tuning parameter.  In the case of the SEIR model, the maximum possible  
exposure time is simply the transition (from exposed to infectious) time of the initial exposed.
************************/

double ProposedInitialExptime(Vertex init, double *inftimes, double delta)
{	
	return(inftimes[init] - rexp(1/delta));	
}


/************************
ProposedExptime

Generates a proposal for the exposure time for any vertex except
the initial exposed.  The proposal is generated uniformly from the
interval of possible exposure times for the node.  The proposal is
used in a Hastings update step.
*************************/

double ProposedExptime(Vertex j, Vertex *transtree, double *inftimes)
{
	double maxtime = inftimes[j], mintime = inftimes[transtree[j]];
	
	return(unif_rand() * (maxtime - mintime) + mintime);	/* Sample uniformly on (mintime, maxtime) */
}


/************************
ProposedInftime

Generates a proposal for the transition time for a vertex. 
The proposal is generated uniformly from the
interval of possible transition times for the node.  The proposal is
used in a Hastings update step.
*************************/

double ProposedInftime(Vertex j, Vertex *transtree, double *exptimes, double *rectimes, int m)
{
	double mintime = exptimes[j], maxtime = rectimes[j];	/* Contruct the window of possible trans times */
	
	for (Vertex i=1; i<=m; i++)	/* mintime and maxtime are the lower and upper bounds that we're sampling inftimes[j] uniformly from */
		if (transtree[i] == j) maxtime = MIN(maxtime, exptimes[i]);
				
	return(unif_rand() * (maxtime - mintime) + mintime);	/* Sample uniformly on (mintime, maxtime) */
}


/****************************
DrawTransTree

Generates a transmission tree (P) from its full conditional
distribution.  For each infected node, sample 
among the vertices that could possibly be its parents, possibly
giving the putative parent node ("probparentprior") more weight 
(a multiplier specified by "probparentmult").
Node with parent = -999 is the initial exposed.  Only need to
update the first m nodes; these are the only ones that are actually
part of the transmission tree, since these are the ones that were
infected during the course of the epidemic.
******************************/

void DrawTransTree(Vertex *transtree, Network *nwp, double *exptimes, double *inftimes, double *rectimes, int m, int *probparentprior, int *probparentmult)
{
	double maxrand;
	Vertex currpar;

	for (Vertex i=1; i<=m; i++)	/* Find the parent of each Vertex i; the vector of parents comprises the trans. tree */
	{
		maxrand = 0; currpar = -999;		
		if (((nwp->inedges)+i)->value != 0) DrawParent(nwp->inedges,i,i,exptimes,inftimes,rectimes,&maxrand,&currpar,probparentprior[i],*probparentmult);
		if (((nwp->outedges)+i)->value != 0) DrawParent(nwp->outedges,i,i,exptimes,inftimes,rectimes,&maxrand,&currpar,probparentprior[i],*probparentmult);
		transtree[i] = currpar;
	}
}


/****************************
DrawParent

Find the parent of a particular node (orig).  Sample from 
all of the node's potential parents.  Uses a recursive algorithm.
Tried using the distribution of the max order stat., but this
method is actually slower due to the pow() function. 
******************************/

void DrawParent(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *rectimes, double *maxrand, Vertex *currpar, int priorparentnode, int probmult)
{
	double posspar, temprand;
	
	if (x != 0) 
	{
		DrawParent(edges, orig, (edges+x)->left, exptimes, inftimes, rectimes, maxrand, currpar, priorparentnode, probmult);
		if ( (inftimes[(edges+x)->value] < exptimes[orig]) && (exptimes[orig] < rectimes[(edges+x)->value])  )
		{
			posspar = unif_rand();		/* For each Vertex that could be the parent of the current node, we generate U; the vertex with the largest U is chosen as the parent  */
			if (((edges+x) -> value) == priorparentnode)			/* Give extra weight to the most likely parent node */
				for (int count = 2; count <= probmult; count++)
				{
					temprand = unif_rand();
					posspar = MAX(posspar,temprand);
				}
				
			if (posspar > *maxrand)		/* If no Vertices could possibly be the parent, then the parent  will stay at -999, which indicates that it's the initial infected */		
			{
				*currpar = (edges+x) -> value;
				*maxrand = posspar;
			}
		}
		DrawParent(edges, orig, (edges+x)->right, exptimes, inftimes, rectimes, maxrand, currpar, priorparentnode, probmult);
	}
}


/*************************
CalcEdgeProb

Calculates the probability of an edge between two nodes,
given the current values of the parameters, and the 
dyadic covariates.
*************************/

double CalcEdgeProb(int dyadcovindex, double *dyadcovs, double *eta, int etapars, int N)
{
	double tempsum = 0;
	int k, totaldyads = (N)*(N-1)/2;
	
	for (k = 0; k < etapars; k++)
		tempsum += eta[k] * dyadcovs[k*totaldyads + dyadcovindex];
		
	return( exp(tempsum) / ( 1 + exp(tempsum) ) );
}


/**************************
DrawGraph

Sample the contact network (G) from its full conditional, given all of the other 
parameters (including P).  Changing G also changes the total infectious 
pressure applied during the epidemic (A), so this is modified as well.
We cycle through the dyads in a random order -- testing indicated that this
runs about 15% quicker than updating the dyads in numerical order.
******************************/

void DrawGraph(Network *nwp, Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, double currbeta, double *dyadcovs, int *dyadindex1, int *dyadindex2, double *eta, int etapars, double *A)
{
	double u, v, pij; 
	Vertex i, j, k, priorinf, latterinf, dyadcovindex;

	for (k=1; k<=((nwp->N)*((nwp->N)-1)/2); k++)
	{	
		i = dyadindex1[k]; 
		j = dyadindex2[k];

		dyadcovindex =  GetDyadIndex(i, j, (nwp->N));
		
		priorinf = MININDEX(inftimes,i,j);		/* Set priorinf to the node (i or j) with the earlier infection time, and latterinf to the node with the later infection time */
		latterinf = MAXINDEX(inftimes,i,j);		/* If the two times are equal, these macros will set priorinf to j and latterinf to i */
		
		pij = CalcEdgeProb(dyadcovindex, dyadcovs, eta, etapars, nwp->N); 		/* Calculate probability of edge between nodes i and j) */
		
		u = exp(-currbeta*(MAX(MIN(rectimes[priorinf], exptimes[latterinf])  - inftimes[priorinf] , 0)));	/* Construct the probability of an edge existing */
		v = u*pij / (1 - pij + u*pij);

		if ((unif_rand() < v) || (transtree[latterinf] == priorinf))	/* Randomly determine if the edge exists - the edge will exist with prob. u, unless it's part of the trans tree in which case it always exists */
			*A += AddEdgeToTrees(i, j, nwp) * MAX(MIN( exptimes[latterinf], rectimes[priorinf] ) - inftimes[priorinf] , 0);	/* If necessary, add the edge and the associated infectious pressure */
		else
			*A -= DeleteEdgeFromTrees(i, j, nwp) * MAX(MIN( exptimes[latterinf], rectimes[priorinf] ) - inftimes[priorinf] , 0);	/* If necessary, delete the edge and the associated infectious pressure */						
	}	
}


/***********************
GetRandomOrder

Used to determine the (uniformly random) order in which
to update the exposure/infection times.  Returns a vector with the random 
order of the infecteds.  When updating the exposure times, we exclude
the initial exposed, because it's updated separately.  Also use this
routine to set the random order to cycle through the dyads.
************************/

void GetRandomOrder (Vertex *order, Vertex initexp, int includeinit, int m)
{
	Vertex i, j, temp;
	int count = 1;
	
	for (i=1; i<=m; i++)	/* First set up an ordered vector containing all nodes except possibly initial exposed */
		if ((i != initexp) || (includeinit == 1))
		{
			order[count] = i;
			count++;
		}

	for(i=m-(1-includeinit); i>1; i--)
	{								/* Knuth shuffle for generating a random permutation vector */
		j = floor(unif_rand() * i) + 1;
		temp = order[j];
		order[j] = order[i];
		order[i] = temp;
	}
}


/***************************
 GetDyadIndex
  
 Used to find out which group a particular dyad belongs to.  Takes the two node indices
 of the individuals in the dyad and returns the index in the dyad group array corresponding
 to the dyad.  If the first index is a and the second index is b, then index for the dyad
 is [N*(a-1)] + [(b-a)-1] - [a*(a-1)/2].  
 To get the group that this dyad belongs to, we then have to look up the dyad in the dyad 
 group array, using this index.
 ***************************/

int GetDyadIndex(int nodeindex1, int nodeindex2, int N)
{
	int a = MIN(nodeindex1, nodeindex2), b = MAX(nodeindex1, nodeindex2);
	
	return((N * (a - 1)) + ((b - a) - 1) - (a * (a - 1) / 2));
}


/************************
 CreateRandomDyadOrder

 When updating the graph (G), we update each dyad individually.
 This function creates a random ordering of dyads that we use to
 cycle through all of the dyads.  (We don't want to update in order
 because doing so will create potentially very unbalanced tree structures
 -- we store the edges for a node in binary trees.)  dyadindex1 will
 contain the first node in each dyad (the smaller number), and dyadindex2 
 will contain the second node (larger number) in the dyad.
**************************/
 
void CreateRandomDyadOrder(int *dyadindex1, int *dyadindex2, int N)
{
	int i, j, totaldyads = (N)*(N-1)/2;	/* There are N-choose-2 dyads */
	
	GetRandomOrder(dyadindex2, 0, 1, totaldyads);	/* Create a random permutation of the numbers 1..N-choose-2 */

	/* Convert the numbers 1 .. N-choose-2 to dyad indices */
	for (i = 1; i <= totaldyads; i++) dyadindex1[i] = 1;	

	for (j = 1; j <= (N-2); j++)	
		for (i = 1; i <= totaldyads; i++)
			if ( (dyadindex2[i] > (N-j)) && (dyadindex1[i] >= j) )
			{
				dyadindex1[i]++;
				dyadindex2[i] -= (N-j);
			}

	for (i = 1; i <= totaldyads; i++) dyadindex2[i] += dyadindex1[i];
}


/****************************
LogLikelihood

Calculates the log-likelihood function, excluding the indicator of whether
or not the transmission tree is legal.
****************************/

double LogLikelihood(double beta, double thetai, double ki ,double thetae, double ke, int m, double A, double B, double Bln, double C, double Cln)
{
	double llkd = 0;
	
	/* Add contribution due to recovery process */
	llkd -= m*(ki*log(thetai) + Rf_lgammafn(ki));
	llkd += (ki - 1) * Cln - C/thetai;
	
	/* Add contribution due to transition (exposed to infected) process */
	llkd -= m*(ke*log(thetae) + Rf_lgammafn(ke));
	llkd += (ke - 1) * Bln - B/thetae;

	/* Add contribution due to exposure process */
	llkd += (m-1)*log(beta) - beta*A;

	return(llkd);	
}


/***********************
IsTreeLegal

Calculates the part of the likelihood corresponding
to the indicator function of the transmission tree is legal.
Returns a value of 1 if the proposed transmission
tree, exposure, and infection times are legal, and a value of 0 otherwise.
**************************/

int IsTreeLegal(double *exptimes, double *inftimes, double *rectimes, int *transtree, Network *nwp, int m)
{	
	Vertex i;
	int j = 0, initialexp = 1;
		
	/* Check if there's only one initial exposed (parent = - 999) */
	for (i=1; i<=m; i++) 	
		if (transtree[i] == -999) 
		{
			j++;
			initialexp = i;
		}	
	if (j > 1) return(0);
	
	/* Make sure initial infected has lowest exposure time */
	for (i=1; i<=m; i++)
		if (i != initialexp)
			if (exptimes[i] < exptimes[initialexp])
				return(0);
			
	/* Check each infection (other than initial) to verify it's indeed possible, given the tree and times */
	for (i=1; i<=m; i++)
		if (i != initialexp)
			if ( (inftimes[transtree[i]] > exptimes[i]) || (exptimes[i] > rectimes[transtree[i]]) || ( !EDGE_EXISTS(transtree[i],i,nwp) ) ) 
				return(0);

	/* Make sure E_j <= I_j <= R_j for each node j */		
	for (i=1; i<=m; i++)
		if ( (exptimes[i] > inftimes[i]) || (inftimes[i] > rectimes[i]) )
			return(0);		
			
	/* If we've made it to this point, tree is legal */
	return(1);	
}


/**********************
AdjustAiiExpTime

Calculates the change in infectious pressure (A) associated with a 
proposed change in the infection time for one node.
Uses a recursive algorithm.
************************/

void AdjustAiiExpTime(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *pexptimes, double *inftimes, double *rectimes, double *propA)
{
	Vertex priorinf, latterinf;
	
	if (x != 0) 
	{
		AdjustAiiExpTime(edges, orig, (edges+x)->left, exptimes, pexptimes, inftimes, rectimes, propA);

		priorinf = MININDEX(inftimes,orig,(edges+x)->value);
 		latterinf = MAXINDEX(inftimes,orig,(edges+x)->value);
		*propA -= MAX(MIN(exptimes[latterinf], rectimes[priorinf]) - inftimes[priorinf],0); /* subtract out previous contribution to propA, based on current inf time */
		*propA += MAX(MIN(pexptimes[latterinf], rectimes[priorinf]) - inftimes[priorinf],0); /* add in new contribution to propA, based on proposed inf time */
		
		AdjustAiiExpTime(edges, orig, (edges+x)->right, exptimes, pexptimes, inftimes, rectimes, propA);
	}

}


/******************************
AdjustABExpTime

Calculates the changes in infectious (A) and transition (B, Bln) 
 pressures associated with a proposed change in the infection 
 time for one node.  Uses a recursive algorithm twice - once 
 for the inedges and once for the outedges.
*****************************/

void AdjustABExpTime(Network *nwp, Edge orig, double *exptimes, double *pexptimes, double *inftimes, double *rectimes, double *propA, double *propB, double *propBln)
{
	/* Adjust A */
	if (((nwp->outedges)+orig)->value != 0) AdjustAiiExpTime(nwp->outedges,orig,orig,exptimes,pexptimes,inftimes,rectimes,propA);
	if (((nwp->inedges)+orig)->value != 0) AdjustAiiExpTime(nwp->inedges,orig,orig,exptimes,pexptimes,inftimes,rectimes,propA);															

	/* Adjust B */
	*propB += exptimes[orig] - pexptimes[orig];
	*propBln += log(inftimes[orig] - pexptimes[orig]) - log(inftimes[orig] - exptimes[orig]);	
}


/**********************
AdjustAiiInfTime

Calculates the change in infectious pressure (A) associated with a 
proposed change in the transition time for one node.
Uses a recursive algorithm.
************************/

void AdjustAiiInfTime(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *pinftimes, double *rectimes, double *propA)
{
	Vertex priorinf, latterinf;
	
	if (x != 0) 
	{
		AdjustAiiInfTime(edges, orig, (edges+x)->left, exptimes, inftimes, pinftimes, rectimes, propA);

		priorinf = MININDEX(inftimes,orig,(edges+x)->value);
 		latterinf = MAXINDEX(inftimes,orig,(edges+x)->value);
		*propA -= MAX(MIN(exptimes[latterinf], rectimes[priorinf]) - inftimes[priorinf],0); /* subtract out previous contribution to propA, based on current inf time */

		priorinf = MININDEX(pinftimes,orig,(edges+x)->value);
 		latterinf = MAXINDEX(pinftimes,orig,(edges+x)->value);
		*propA += MAX(MIN(exptimes[latterinf], rectimes[priorinf]) - pinftimes[priorinf],0); /* add in new contribution to propA, based on proposed inf time */
		
		AdjustAiiInfTime(edges, orig, (edges+x)->right, exptimes, inftimes, pinftimes, rectimes, propA);
	}

}


/******************************
AdjustABCInfTime

Calculates the change in infectious (A), transition (B, Bln), and removal (C, Cln) pressures
 associated with a proposed change in the transition time for one node.  
Uses a recursive algorithm twice - once for the inedges and once for the outedges.
*****************************/

void AdjustABCInfTime(Network *nwp, Edge orig, double *exptimes, double *inftimes, double *pinftimes, double *rectimes, 
	double *propA, double *propB, double *propBln, double *propC, double *propCln)
{
	/* Adjust A */
	if (((nwp->outedges)+orig)->value != 0) AdjustAiiInfTime(nwp->outedges,orig,orig,exptimes,inftimes,pinftimes,rectimes,propA);
	if (((nwp->inedges)+orig)->value != 0) AdjustAiiInfTime(nwp->inedges,orig,orig,exptimes,inftimes,pinftimes,rectimes,propA);															
	
	/* Adjust B and C, and their log versions */
	*propB += pinftimes[orig] - inftimes[orig];
	*propBln += log(pinftimes[orig] - exptimes[orig]) - log(inftimes[orig] - exptimes[orig]);
	*propC += inftimes[orig] - pinftimes[orig];
	*propCln += log(rectimes[orig] - pinftimes[orig]) - log(rectimes[orig] - inftimes[orig]);
}


/*********************
AdjustAClnKappa

 Calculates the change in infectious (A) and log-removal (Cln) pressures
 associated with a change in the initial exposed (kappa). 
*********************/

void AdjustAClnKappa(Network *nwp, Edge initexp, Edge pinitexp, double *exptimes,  double *pexptimes, 
	double *inftimes, double *pinftimes, double *rectimes, double *propA, double *propCln)
{
	/* Adjust A */
	if (((nwp->outedges)+initexp)->value != 0) AdjustAiiExpTime(nwp->outedges,initexp,initexp,exptimes,pexptimes,inftimes,rectimes,propA);
	if (((nwp->inedges)+initexp)->value != 0) AdjustAiiExpTime(nwp->inedges,initexp,initexp,exptimes,pexptimes,inftimes,rectimes,propA);																
	
	if (((nwp->outedges)+pinitexp)->value != 0) AdjustAiiExpTime(nwp->outedges,pinitexp,pinitexp,exptimes,pexptimes,inftimes,rectimes,propA);
	if (((nwp->inedges)+pinitexp)->value != 0) AdjustAiiExpTime(nwp->inedges,pinitexp,pinitexp,exptimes,pexptimes,inftimes,rectimes,propA);															

	if (((nwp->outedges)+initexp)->value != 0) AdjustAiiInfTime(nwp->outedges,initexp,initexp,pexptimes,inftimes,pinftimes,rectimes,propA);
	if (((nwp->inedges)+initexp)->value != 0) AdjustAiiInfTime(nwp->inedges,initexp,initexp,pexptimes,inftimes,pinftimes,rectimes,propA);															
	
	if (((nwp->outedges)+pinitexp)->value != 0) AdjustAiiInfTime(nwp->outedges,pinitexp,pinitexp,pexptimes,inftimes,pinftimes,rectimes,propA);
	if (((nwp->inedges)+pinitexp)->value != 0) AdjustAiiInfTime(nwp->inedges,pinitexp,pinitexp,pexptimes,inftimes,pinftimes,rectimes,propA);															
	
	/* Adjust Cln */
	*propCln += log(rectimes[initexp] - pinftimes[initexp]) - log(rectimes[initexp] - inftimes[initexp]);
	*propCln += log(rectimes[pinitexp] - pinftimes[pinitexp]) - log(rectimes[pinitexp] - inftimes[pinitexp]);	
}


/************************
LogkPrior

Calculates the log of the density of the prior distribution for either of the k 
parameters (or really any other parameter), evaluated at the argument x. 
The prior distribution parameter(s) are passed in as well.  
*************************/

double LogkPrior(double x, double *kprior, int *kpriordist)
{
	double logdens = 0;
	
	if (*kpriordist == 1)	/* Gamma prior */
		logdens = Rf_dgamma(x,kprior[0],kprior[1],1); /* Calculate log-density of prior distribution, evaluated at x */
	else if (*kpriordist == 0)	/* Uniform prior */
		logdens = 0;		/* We've already checked to make sure the proposal is in the bounds, so nothing to do here */
	else if (*kpriordist == 2)	/* Flat prior */
		logdens = 0;		/* Nothing to do here */
	else if (*kpriordist == 3)  /* Normal prior */
		logdens = dnorm(x,kprior[0],kprior[1],1);
	
	return(logdens);
}


/**************
 LogEdgeCalc
 
 Recursive routine to calculate the numerator of the log-density
 of the graph, since this only depends on the dyads in which an
 edge is present.  Only called by LogGraphCalc().
 **************/

void LogEdgeCalc(TreeNode *edges, int orig, int x, double *eta, int etapars, double *dyadcovs, int totaldyads, double *logdensptr, int N)
{
	int q, dyadindex;
	
	if (x != 0) 
	{
		LogEdgeCalc(edges, orig, (edges+x)->left, eta, etapars, dyadcovs, totaldyads, logdensptr, N);
		
		dyadindex = GetDyadIndex(orig, (edges+x)->value, N);
		for (q = 0; q < etapars; q++)
			(*logdensptr) += eta[q] * dyadcovs[q*totaldyads + dyadindex];				
		
		LogEdgeCalc(edges, orig, (edges+x)->right, eta, etapars, dyadcovs, totaldyads, logdensptr, N);
	}
}


/************************
LogGraphCalc

Calculates the log of the ERGM pdf, given the current state
 of the graph and the eta parameters.  Uses a recursive algorithm
 to traverse the graph for the numerator; cycles through all dyads
 to calculate the denominator.
*************************/

double LogGraphCalc(Network *nwp, double *eta, int etapars, double *dyadcovs)
{
	int i, j, q, dyadindex, totaldyads = (int) (nwp->N) * ( (nwp->N) - 1 ) / 2;
	double tempsum, logdens = 0, *logdensptr = &logdens;
	
	/* Calculate numerator */
	for(i=1; i<=(nwp->N); i++) 
		if (((nwp->outedges) + i) -> value != 0)
			LogEdgeCalc(nwp->outedges, i, i, eta, etapars, dyadcovs, totaldyads, logdensptr, (nwp->N));
	
	/* Calculate denominator */
	for (i = 1; i < (nwp->N); i++)
	{
		for (j = i + 1; j <= (nwp->N); j++)
		{			
			tempsum = 0;
			dyadindex = GetDyadIndex(i, j, (nwp->N));
			for (q = 0; q < etapars; q++)
				tempsum += eta[q] * dyadcovs[q*totaldyads + dyadindex];

			logdens -= log(1 + exp(tempsum));
		}	
	}
	return(logdens);
}


/* **************************** */
/* MAIN MCMC FUNCTION */
/* **************************** */

/*************************
epigraphmcmcc

Main MCMC function to produce a sample of the parameters.  
************************/

void epigraphmcmcc (double *etime, double *itime, double *rtime, int *etapars, double *dyadcovs, int *nsamp, int *thinning,  double *bprior, double *tiprior, double *teprior, double *etaprior, double *kiprior,
	double *keprior, int *ninf, int *initN, double *initbeta, double *initthetai, double *initki, double *initthetae, double *initke, double *initeta,  int *bpriordist, 
	int *tipriordist, int *tepriordist, int *kipriordist, int *kepriordist, int *etapriordist, double *etapropsd, int *accept, int *propose, double *allkd, double *abeta, double *athetai, double *aki, double *athetae, double *ake ,double *aeta, 
	int *ainit, double *ainitexptime, double *aexptimes, double *ainftimes, int *atranstree, int *extrathinning, int *inferEtimes, int *inferItimes, int *parentprior, int *probparentmult,
	int *verbose, int *burnin, int *numsamp, int *numsamptimes)
{ 

	GetRNGstate();  /* R function enabling uniform RNG */ 
	
	/* VARIABLE DECLARATIONS */
	
	Network gr, *nwp=&gr;								/* This is the network that holds the edges between individuals in the population */
	int m = *ninf, N = *initN; 							/* m = # of infecteds (considered fixed and known); note that this variable holds the same constant as nwp->ninfnodes, just have it here for convenience.  Similar story for N. */
	double llkd = 0, pllkd = 0, currlogprior, proplogprior;					/* Variables to hold current and proposed values of the log-likelihood and log-priors */
	Vertex *transtree = (Vertex *) malloc((N+1) * sizeof(Vertex));	/* transtree[i] holds the label of the node that infected node i in the current infection tree... */
	Vertex *ptranstree = (Vertex *) malloc((N+1) * sizeof(Vertex));		/* ... ptranstree is the proposal version */
	Vertex *orderrand = (Vertex *) malloc((m+1) * sizeof(Vertex));		/* Update the inf/exp times in random order - orderrand holds this random order (only need to update the infecteds) */
	int *dyadindex1 = (int *) malloc(((N)*(N-1)+1) * sizeof(int));		/* Holds the first node of each dyad - this is used to cycle through the dyads in random order */
	int *dyadindex2 = (int *) malloc(((N)*(N-1)+1) * sizeof(int));		/* Holds the second node of each dyad - this is used to cycle through the dyads in random order */
	int ok = 1, zz, i, j, currcount, propcount;					/* dummy & counter variables */
	long iter;
	time_t last_t, curr_t;
	double currbeta = *initbeta, currthetai = *initthetai, currthetae = *initthetae; /* Starting values for some of the parameters... */
	double currki = *initki, currke = *initke;	/*... starting values for more parameters */
	double *eta = (double *) malloc((*etapars) * sizeof(double));		/* Current values of eta parameters (one for each dyadic parameter) */
	double *propeta = (double *) malloc((*etapars) * sizeof(double));		/* Proposed values of eta parameters (one for each dyadic parameter) */
	double A, propA, B = 0, propB, Bln = 0, propBln, C = 0, propC, Cln = 0, propCln;				/* Variables to keep track of infectious, transition and removal pressures, and their log versions */
	double *exptimes = (double *) malloc((N+1) * sizeof(double));				/* Use indices 1..N and leave 0 unused to be consistent with Network object and R (1..m are the infecteds) */ 
	double *pexptimes = (double *) malloc((N+1) * sizeof(double));			/* These hold the current and proposed exposure times */
	double *inftimes = (double *) malloc((N+1) * sizeof(double));				
	double *pinftimes = (double *) malloc((N+1) * sizeof(double));			/* These hold the current and proposed infection times */
	double *rectimes = (double *) malloc((N+1) * sizeof(double));			/* These hold the recovery times, considered fixed and known */
	int *probparentprior = (int *) malloc((m+1) * sizeof(int));			/* These hold the prior belief about most likely parent for each node */
	Vertex initexp, pinitexp;						/* These hold the current and proposed identity of the initial exposed (Kappa)	*/
	double propthetai, propthetae, propbeta, propki, propke;	/* Variables to hold proposal values for parameters */
	double delta = 0.5, uniwidth = 5;	// Tuning parameters - maybe will pass them in or figure out a way to automatically set them - these values seem to work reasonably well...
	double bwidth = (bprior[1] - bprior[0])/uniwidth, tiwidth = (tiprior[1] - tiprior[0])/uniwidth, tewidth = (teprior[1] - teprior[0])/uniwidth;
	double kiwidth = (kiprior[1] - kiprior[0])/uniwidth, kewidth = (keprior[1] - keprior[0])/uniwidth;	// Variables used in case of uniform priors
	int maxmove = 11;	// Sets the number of updates in each sweep of the algorithm
    
	/* INITIALIZE VARIABLES */
	
	if (*verbose == 1) Rprintf("Initializing variables for MCMC. \n");
			
	for (i = 0; i < N; i++) 
	{
		exptimes[i+1] = etime[i];		/* Copy initial exposure times (since they might change during the procedure)	*/
		inftimes[i+1] = itime[i];	/* Copy initial infection times (since they might change during the procedure)	*/
		rectimes[i+1] = rtime[i];	/* Copy recovery times just for indexing consistency  */
		if (i < m) probparentprior[i+1] = parentprior[i];	/* Copy prior beliefs about parents just for indexing consistency  */
	}

	for (i = (m+1); i <= N; i++)	/* Set permanent values of proposed I and E times for susceptibles -- these should never change */
	{
		pexptimes[i] = exptimes[i];
		pinftimes[i] = inftimes[i];
	}
	
	for (i = 0; i < (*etapars); i++)	/* Eta parameters are indexed 0 .. (etapars - 1) */
		eta[i] = initeta[i];
	
 	if (InitializeTransTree(transtree, exptimes, inftimes, rectimes, m, N, &initexp, &A) != 1) 	/* Initializes transmission tree, determines the initial infected, and calculates the amount of infectionus pressure (A) due to the transmission tree */
	{
		Rprintf("Faulty exposure/infection/recovery time data.  Cannot build legal initial transmission tree.  Aborting routine. \n");
		return;
	}
		
	gr = InitializeNetworkFromTree(transtree, m, N);	/* Initialize network to only contain the edges from the transmission tree */		
	if (!IsTreeLegal(exptimes,inftimes,rectimes,transtree,nwp,m)) 		/* If initial transmission tree isn't legal, then the set of infection times and recovery times are somehow faulty, so we print an error and abort */
	{
		Rprintf("Faulty exposure/infection/recovery time data.  Cannot build legal initial transmission tree.  Aborting routine. \n");
		return;
	}

	for (i=1; i <= m; i++)	/* Calculate initial transition and removal pressures and their log versions -- only the infecteds add to these pressures */
	{
		C += rectimes[i] - inftimes[i];
		Cln += log(rectimes[i] - inftimes[i]);
		B += inftimes[i] - exptimes[i];
		Bln += log(inftimes[i] - exptimes[i]);
	}

	CreateRandomDyadOrder(dyadindex1, dyadindex2, N);		/* Create random order to update dyads with updating the graph -- will use the same order throughout */
	
	DrawGraph(nwp, transtree, exptimes, inftimes, rectimes, currbeta, dyadcovs, dyadindex1, dyadindex2, eta, *etapars, &A);  /* Draw initial graph from its full conditional dist. */	
	
	/* BEGIN MAIN MCMC LOOP */
	
	if (*verbose == 1) Rprintf("Beginning MCMC. \n");
	
	for (iter = - ( (long) (*burnin) * (maxmove) ); iter < ( (long) (*nsamp) * (maxmove) ); iter++)
	{
		if (iter == 0) last_t = time(NULL); // Start clock after burn-in is complete
		zz = labs(iter % (maxmove));		// Cyclical updates - can change this line to update different sets of parameters or change the order of updates
		if (zz <= 7) propose[zz] ++;		// Not currently keeping track of propose / accept info for E/I times or Kappa
		switch(zz)		/* Choose which item (parameter) to update */
		{
			
		case 0:	/* Update Transmission Tree (P)	*/
			
			DrawTransTree(transtree, nwp, exptimes, inftimes, rectimes, m, probparentprior, probparentmult);	/* transtree will have the updated transmission tree */
			accept[zz]++;
			break;

		case 1:	/* Update eta parameters */
						
			/* Doing a block update on eta */
				
			ok = 1;
			for (i = 0; i < (*etapars); i++)	
				propeta[i] = eta[i] + norm_rand() * etapropsd[i];

			pllkd = LogGraphCalc(nwp, propeta, *etapars, dyadcovs);
			llkd = LogGraphCalc(nwp, eta, *etapars, dyadcovs);
			proplogprior = 0;
			currlogprior = 0;
			for (i = 0; i < (*etapars); i++)
				proplogprior += LogkPrior(propeta[i], &etaprior[2*i], &etapriordist[i]);
			for (i = 0; i < (*etapars); i++)
				currlogprior += LogkPrior(eta[i], &etaprior[2*i], &etapriordist[i]);

			if ( log(unif_rand()) > (pllkd + proplogprior - llkd - currlogprior) )
				ok = 0;
								
			if (ok == 1)
			{
				for (i = 0; i < (*etapars); i++)
					eta[i] = propeta[i];						
				accept[zz]++;					
			}
			
			break;
			
		case 2:	/* Update graph (G) */
		
			DrawGraph(nwp, transtree, exptimes, inftimes, rectimes, currbeta, dyadcovs, dyadindex1, dyadindex2, eta, *etapars,  &A);
			accept[zz]++;
			break;						

		case 3:	/* Update Beta	*/
			
			if (*bpriordist == 0)	/* Uniform prior for Beta -- use MH random walk update */
			{
				ok = 1;			
				propbeta = currbeta + unif_rand()*bwidth - bwidth/2;
				if ((propbeta > bprior[1]) || (propbeta < bprior[0])) ok = 0;
				if (ok == 1)
				{
					llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	
					pllkd = LogLikelihood(propbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	
					if (log(unif_rand()) > pllkd - llkd) ok = 0;				
				}
				if (ok == 1)
				{
					accept[zz]++;
					currbeta = propbeta;
				}
			} else if(*bpriordist == 1)	/* Gamma (conjugate) prior for Beta -- use Gibbs update */
			{
				currbeta = rgamma( bprior[0] + m - 1, 1 / ( ( 1 / bprior[1] ) + A ) );
				accept[zz]++;
			}	
			break;
						
		case 4:	/* Update Theta_I */
		
			if (*tipriordist == 0)	/* Uniform prior for Theta_I -- use MH random walk update */
			{
				ok = 1;			
				propthetai = currthetai + unif_rand()*tiwidth - tiwidth/2;
				if ((propthetai > tiprior[1]) || (propthetai < tiprior[0])) ok = 0;
				if (ok == 1)
				{
					llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	
					pllkd = LogLikelihood(currbeta,propthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	
					if (log(unif_rand()) > pllkd - llkd) ok = 0;				
				}
				if (ok == 1)
				{
					accept[zz]++;
					currthetai = propthetai;
				}
			} else if(*tipriordist == 1)	/* Inverse Gamma (conjugate) prior for Theta_I -- use Gibbs update */
			{
				currthetai = 1/rgamma( tiprior[0] + currki*m, 1/ (tiprior[1] + C) );			
				accept[zz]++;
			} 	
			break;
			
		case 5:  /* Update k_I parameter for removal process */	
			
			ok = 1;
			propki = currki + unif_rand()*kiwidth - kiwidth/2;	/* Random walk MH update */

			if (*kipriordist == 1)	/* Check to see if we've proposed a value outside the prior bounds before proceeding further */
			{
				if (propki < 0) 
					ok = 0;
			} else if (*kipriordist == 0)
			{
				if ( (propki < kiprior[0] ) || ( propki > kiprior[1] ) ) 
					ok = 0;
			}

			if (ok == 1)
			{
				currlogprior = LogkPrior(currki, kiprior, kipriordist);		/* calculate current prior for k_I */
				proplogprior = LogkPrior(propki, kiprior, kipriordist);		/* calculate proposed prior for k_I */
				llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* calculate current log likelihood */
				pllkd = LogLikelihood(currbeta,currthetai, propki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* calculate proposed log likelihood */
				if (log(unif_rand()) > (pllkd + proplogprior - llkd - currlogprior) ) ok = 0;
			}
			if (ok == 1)
			{
				currki = propki;
				accept[zz]++;
			}
			
			break;

		case 6:		/* Update theta_E */

			if (*tepriordist == 0)	/* Uniform prior for Theta_E -- use MH Random Walk update */
			{
				ok = 1;			
				propthetae = currthetae + unif_rand()*tewidth - tewidth/2;
				if ((propthetae > teprior[1]) || (propthetae < teprior[0])) ok = 0;
				if (ok == 1)
				{
					llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	
					pllkd = LogLikelihood(currbeta,currthetai, currki, propthetae, currke, m, A, B, Bln, C, Cln); 	
					if (log(unif_rand()) > pllkd - llkd) ok = 0;				
				}
				if (ok == 1)
				{
					accept[zz]++;
					currthetae = propthetae;
				}
			} else if(*tepriordist == 1)	/* Inverse Gamma (conjugate) prior for Theta_E -- use Gibbs update */
			{
				currthetae = 1/rgamma( teprior[0] + currke*m, 1/ (teprior[1] + B) );			
				accept[zz]++;
			} 	
			
			break;
			
		case 7:		/* Update k_E */

			ok = 1;
			propke = currke + unif_rand()*kewidth - kewidth/2;		/* MH Random walk update */

			if (*kepriordist == 1)	/* Check to see if we've proposed a value outside the prior bounds before proceeding further */
			{
				if (propke < 0) 
					ok = 0;
			} else if (*kepriordist == 0)
			{
				if ( (propke < keprior[0] ) || ( propke > keprior[1] ) ) 
					ok = 0;
			}

			if (ok == 1)
			{
				currlogprior = LogkPrior(currke, keprior, kepriordist);		/* calculate current prior for k_E */
				proplogprior = LogkPrior(propke, keprior, kepriordist);		/* calculate proposed prior for k_E */
				llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* calculate current log likelihood */
				pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, propke, m, A, B, Bln, C, Cln); 	/* calculate proposed log likelihood */
				if (log(unif_rand()) > (pllkd + proplogprior - llkd - currlogprior) ) ok = 0;
			}
			if (ok == 1)
			{
				currke = propke;
				accept[zz]++;
			}
			
			break;
			
		case 8:		/* Update infection times (I) */

			if (*inferItimes == 0) break;			/* If the infection times are known and fixed, then there's nothing to do here. */
			
			for (i=1; i <=m; i++) pinftimes[i] = inftimes[i];		/* Initialize proposed infection times - will start with current infection times and update one at a time*/
			
			llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* Calculate current log likelihood */

			GetRandomOrder(orderrand, initexp, 1, m);	/* Orderrand contains the (random) order to update infection times  */
			
			for (i=1; i <= m; i++)	/* Cycle through the (infected) nodes  */
			{
				ok = 1;
				j = orderrand[i];	/* Update infection time for Vertex j */
				
				pinftimes[j] = ProposedInftime(j,transtree, exptimes, rectimes, m);	/* Get proposed infection time for vertex j */
				if (IsTreeLegal(exptimes,pinftimes,rectimes,transtree,nwp,m) == 0) ok = 0;	/* If proposed time doesn't result in a legal tree, no need to calculate the rest of the likelihood */
								
				if (ok == 1) 		/* We've proposed a legal infection time, now check to see if we accept it */
				{
					propA = A; propB = B; propBln = Bln; propC = C; propCln = Cln;
					AdjustABCInfTime(nwp, j, exptimes, inftimes, pinftimes, rectimes, &propA, &propB, &propBln, &propC, &propCln);
					pllkd = LogLikelihood(currbeta, currthetai, currki, currthetae, currke, m, propA, propB, propBln, propC, propCln);
					if (log(unif_rand()) > (pllkd - llkd) ) ok = 0;				
				}
				
				if (ok == 1) 		/* Accepted */
				{		
					A = propA; B = propB; Bln = propBln; C = propC; Cln = propCln;
					inftimes[j] = pinftimes[j];							
					llkd = pllkd;
				}
				else pinftimes[j] = inftimes[j];	/* If not accepted, need to roll back proposal before proposing next infection time*/
			}
						
			break;
			
		case 9:	/* Update Exposure Times (E) */

			if (*inferEtimes == 0) break;			/* If the exposure times are known and fixed, then there's nothing to do here. */

			for (i=1; i <=m; i++) pexptimes[i] = exptimes[i];		/* Initialize proposed exposure times - will start with current exposure times and update one at a time */
			
			llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* Calculate current log likelihood */

			GetRandomOrder(orderrand, initexp, 0, m);	/* Orderrand contains the (random) order to update exposure times except for initial exposed */
			
			for (i=1; i<m; i++)	/* Cycle through the exposeds, other than the initial exposed */
			{
				ok = 1;
				j = orderrand[i];	/* Update exposure time for Vertex j */
				pexptimes[j] = ProposedExptime(j, transtree, inftimes);	/* Get proposed exposure time for vertex j */
				if (IsTreeLegal(pexptimes, inftimes, rectimes, transtree, nwp, m) == 0) ok = 0;	/* If proposed time doesn't result in a legal tree, no need to calculate the rest of the likelihood */
				
				if (ok == 1) 		/* We've proposed a legal time, now check to see if we accept it */
				{
					propA = A; propB = B; propBln = Bln;
					AdjustABExpTime(nwp, j, exptimes, pexptimes, inftimes, rectimes, &propA, &propB, &propBln);
					pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, propA, propB, propBln, C, Cln);
					if (log(unif_rand()) > (pllkd - llkd) ) ok = 0;
				}
				
				if (ok == 1) 		/* Accepted */
				{		
					A = propA; B = propB; Bln = propBln;
					exptimes[j] = pexptimes[j];							
					llkd = pllkd;
				}
				else pexptimes[j] = exptimes[j];	/* If not accepted, need to roll back proposal before proposing next exposure time*/
			
			}
						
			/* Now update exposure time for orig exposed */
			ok = 1;
			pexptimes[initexp] = ProposedInitialExptime(initexp, inftimes, delta);
			if ((pexptimes[initexp] > 0) || ( !IsTreeLegal(pexptimes,inftimes, rectimes,transtree,nwp,m) )) ok = 0; /* Initial exposure time must be in (-inf, 0)  - if it's not, no need to go further */

			if (ok == 1) 	/* We've proposed a legal time, now check to see if we accept it */
			{
				propB = B; propBln = Bln;
				AdjustABExpTime(nwp, initexp, exptimes, pexptimes, inftimes, rectimes, &propA, &propB, &propBln);
				pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, propB, propBln, C, Cln);	/* Note: proposing a new exposure time for the initial exposed doesn't impact A */
				if ( log(unif_rand()) > pllkd  - llkd  - (delta) * (pexptimes[initexp] - exptimes[initexp]) ) ok = 0;	
			}
			
			if (ok == 1)	/* Accepted */
			{
				exptimes[initexp] = pexptimes[initexp];
				B = propB; Bln = propBln;
			}
			
			break;
			
		case 10:		/* Update which vertex is the initial exposed (Kappa)  */
			
			ok = 1;
				
			if ( (*inferItimes == 0) || (*inferEtimes == 0) ) break;			/* If the infection or exposure times are known and fixed, then there's nothing to do here. */
			
			for (i = 1; i <= m; i++) 
			{								/* Initialize proposal vectors for exposure times, infection times, and trans tree */
				pexptimes[i] = exptimes[i];		
				pinftimes[i] = inftimes[i];
				ptranstree[i] = transtree[i];
			}

			pinitexp = ProposedInitExp(initexp, transtree, m);		/* Get proposed value for vertex to be initial exposed */
			
			pexptimes[initexp] = exptimes[pinitexp];		/* Swap proposal exposure times */
			pexptimes[pinitexp] = exptimes[initexp];
			
			pinftimes[initexp] = inftimes[pinitexp];		/* Swap proposal infection times */
			pinftimes[pinitexp] = inftimes[initexp];			
			
			ptranstree[pinitexp] = -999;			/* Swap proposal direction of initial infection  */
			ptranstree[initexp] = pinitexp;
							
			if ( !IsTreeLegal(pexptimes,pinftimes,rectimes,ptranstree,nwp,m) ) ok = 0;		/* If the proposal produces an illegal tree, no need to calculate the likelihood */		
			
			if (ok == 1) 	/* If tree is legal, see if we accept proposal */
			{
				propA = A; propCln = Cln;

				AdjustAClnKappa(nwp, initexp, pinitexp, exptimes, pexptimes, inftimes, pinftimes, rectimes, &propA, &propCln);
				
				currcount = 0; propcount = 0;
				for (i=1; i<=m; i++)
				{							/* currcount and propcount count the number of infections by the current and proposed initial infectives, respectively - this is needed for Hastings ratio below */
					if (transtree[i] == initexp) currcount++;		
					if (ptranstree[i] == pinitexp) propcount++;
				}
				
				llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* calculate current log likelihood */
				pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, propA, B, Bln,C, propCln);	/* calculate proposed log likelihood */
				
				if (log(unif_rand()) > pllkd + log(currcount) - llkd - log(propcount)) ok = 0;
			}

			if (ok == 1)	/* Accepted proposal, so update values */
			{
				exptimes[initexp] = pexptimes[initexp];
				exptimes[pinitexp] = pexptimes[pinitexp];
				inftimes[initexp] = pinftimes[initexp];
				inftimes[pinitexp] = pinftimes[pinitexp];
				transtree[initexp] = ptranstree[initexp];
				transtree[pinitexp] = ptranstree[pinitexp];				
				initexp = pinitexp;
				A = propA; Cln = propCln;
			}

			break;			
			
		}

		if (*verbose == 1)
			if ( ( ( ( (long) iter ) % ( (long) maxmove * (*nsamp) / 100 ) ) == 0 ) && ( iter >= 0 ) )
			{
				curr_t = time(NULL);
				if (iter == 0)
				{
					if (*burnin > 0) Rprintf("Burn-in complete. \n");	
				}
				else Rprintf("%d of %d MCMC iterations complete (%ld secs/%d iterations). \n", (iter / maxmove), (*nsamp) , (long) difftime(curr_t, last_t) , (int) ( (*nsamp) / 100 ) );
				last_t = curr_t;
			}
		
		if ( (iter % ( (maxmove) * (*thinning) ) == 0) && (iter >= 0) && (floor(iter / ( (maxmove) * (*thinning) ) ) < (*numsamp) ) )
		{  
			j = floor(iter / ( (maxmove) * (*thinning) ) );				/* 	Every thinning sweeps (updates of each parameter), record parameter values */
			athetai[j] = currthetai;				/*	These are storage vectors that are initialized in R and passed (empty) into C	*/		
			abeta[j] = currbeta;					/*	They're passed back to R containing the sample values */
			aki[j] = currki;
			athetae[j] = currthetae;
			ake[j] = currke;
			ainitexptime[j] = exptimes[initexp];
			ainit[j] = initexp;
			allkd[j] = (ok == 1) ? pllkd : llkd;
            
			for (i = 0; i < (*etapars); i++)
				aeta[j*(*etapars)+i] = eta[i];						/* Record eta parameter values */
						
			R_CheckUserInterrupt();			/* Occasionally check for an R user interrupt */
			
			if (*extrathinning > 0)			/* Record exposure and infection time, as well as the transmission tree, if desired */
				if  ( (iter % ( (maxmove) * (*thinning) * (*extrathinning)) == 0) && (floor(iter / ( (maxmove) * (*thinning) * (*extrathinning)) ) < (*numsamptimes) ) )
				{
					j = floor(iter / ( (maxmove) * (*thinning) * (*extrathinning)) );
					for (i=1; i<=m; i++)
					{
						aexptimes[j*m+i - 1] = exptimes[i];
						ainftimes[j*m+i - 1] = inftimes[i];
						atranstree[j*m+i - 1] = transtree[i];
					}					
				}
		}
	}	/* End MCMC loop */
			
	if (*verbose == 1) Rprintf("MCMC complete.\n");
	
	/* Return memory used before leaving */
	NetworkDestroy(nwp);
	free(transtree); free(ptranstree);
	free(exptimes); free(pexptimes);
	free(inftimes); free(pinftimes);
	free(rectimes); free(orderrand);
	free(dyadindex1); free(dyadindex2);
	free(probparentprior);
	free(eta); free(propeta);
	
	PutRNGstate(); /* Must be called after GetRNGstate before returning to R */

}
