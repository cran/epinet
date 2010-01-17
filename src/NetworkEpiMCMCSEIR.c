
/*
 *  File NetworkEpiMCMCSEIR.c
 *
 *  
 */


#include <R.h> 
#include <Rinternals.h> 
#include <R_ext/Rdynload.h> 
#include <R_ext/Utils.h>
#include <Rmath.h>
#include <math.h> 
#include <Rdefines.h>
#include "NetworkEpiMCMCSEIR.h"
//#include "epinetdebugutils.h"



/********************
InitializeTransTree

Initializes the transmission tree, given the input exposure, infection, 
and recovery times.  For each infected node, sample uniformly among 
the vertices that could possibly be its parents.  Called only at the beginning
of the MCMC routine.  Also finds the initial exposed.
*******************/

int InitializeTransTree(Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, int m, Vertex *initexp, double *A)
{

	double maxrand, posspar, lowexptime = exptimes[1];
	Vertex i, j, currpar, lowexp = 1, numinits = 0;

	*A = 0;		
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
		
		if (currpar == -999) numinits++;
		
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
	
	return(numinits);
}


/**************************
 InitializeNetworkFromTree

Initializes the I-I contact nework, given the transmission tree.
Just adds the edges associated with the transmission tree,
by calling NetworkInitialize(). Sets the number of I-S edges
to 0 for each infected, and the total number of S-S edges to 0.
**************************/

Network InitializeNetworkFromTree(Vertex *transtree, int m, Vertex nsus, Vertex N)
{
	Vertex *head = (Vertex *) malloc((m-1) * sizeof(Vertex));
	Vertex *tail = (Vertex *) malloc((m-1) * sizeof(Vertex));
	int i, count = 0;
	
	for(i=1; i<=m; i++)				/* For each node that isn't the initial exposed, add an edge corresponding to its exposure */
		if (transtree[i] != -999)
		{
			head[count] = i;
			tail[count] = transtree[i];
			count++;
		}
								
	return(NetworkInitialize(head,tail,m-1,m,nsus,N));	/* Initialize the network with only the edges corresponding to the trans tree */
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

Generates a proposal for the exposure time of the initial exposed (I_Kappa)
The proposal is an exponential RV subtracted from the minimum possible exposure time
Theta is a tuning parameter.  In the case of the SEIR model, the minimum possible  
exposure time is simply the transition (from exposed to infectious) time of the initial exposed.
************************/

double ProposedInitialExptime(Vertex init, double *inftimes, double delta)
{	
	double mintime = inftimes[init];
	double offset = rexp(1/delta);

	return(mintime - offset);	
}


/************************
ProposedExptime

Generates a proposal for the exposure time for any vertex except
the initial exposed.  The proposal is generated uniformly from the
interval of possible exposure times for the node.  The proposal is
used in a Hastings update step.
*************************/

double ProposedExptime(Vertex j,Vertex *transtree,double *inftimes)
{
	double mintime, maxtime;
	
	maxtime = inftimes[j];    /* Contruct the window of possible exposure times */
	mintime = inftimes[transtree[j]];
				
	return(unif_rand() * (maxtime - mintime) + mintime);	/* Sample uniformly on (mintime, maxtime) */
}


/************************
ProposedInftime

Generates a proposal for the transition time for a vertex. 
The proposal is generated uniformly from the
interval of possible transition times for the node.  The proposal is
used in a Hastings update step.
*************************/

double ProposedInftime(Vertex j, Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, int m)
{
	double mintime, maxtime;
	
	mintime = exptimes[j];    /* Contruct the window of possible trans times */
	maxtime = rectimes[j];
	for (Vertex i=1; i<=m; i++)	/* mintime and maxtime are the lower and upper bounds that we're sampling inftimes[j] uniformly from    */
		if (transtree[i] == j) maxtime = MIN(maxtime, exptimes[i]);
				
	return(unif_rand() * (maxtime - mintime) + mintime);	/* Sample uniformly on (mintime, maxtime) */
}


/****************************
DrawTransTree

Generates a transmission tree (P) from its full conditional
distribution.  For each infected node, sample uniformly 
among the vertices that could possibly be its parents.
Node with parent = -999 is the initial exposed 
******************************/

void DrawTransTree(Vertex *transtree, Network *nwp, double *exptimes, double *inftimes, double *rectimes, int m, int *probparentprior, int *probparentmult)
{
	double maxrand, posspar;
	Vertex currpar;

	for (Vertex i=1; i<=m; i++)	/* Find the parent of each Vertex i; the vector of parents comprises the trans. tree */
	{
		maxrand = 0; posspar = -1; currpar = -999;
		
		if (((nwp->inedges)+i)->value != 0) DrawParent(nwp->inedges,i,i,exptimes,inftimes,rectimes,&maxrand,&currpar,probparentprior[i],*probparentmult);
		if (((nwp->outedges)+i)->value != 0) DrawParent(nwp->outedges,i,i,exptimes,inftimes,rectimes,&maxrand,&currpar,probparentprior[i],*probparentmult);
		transtree[i] = currpar;
	}
}


/****************************
DrawParent

Find the parent of a particular node (orig).  Sample from 
all of the node's potential parents.  Uses a recursive algorithm.
******************************/

void DrawParent(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *rectimes, double *maxrand, Vertex *currpar, int priorparentnode, int probmult)
{
	double posspar, temprand;
	/* double expo; */
	
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
				/* {
				temprand = unif_rand();
				expo = pow( temprand , ( 1.0 / (probmult - 1) ) ); 
				posspar = MAX( posspar, expo );		// Can use the distribution of the max order stat. instead of actually drawing each uniform RV
				} */									// This is perhaps more efficient
				
			if (posspar > *maxrand)		/* If no Vertices could possibly be the parent, then the parent  will stay at -999, which indicates that it's the initial infected */		
			{
				*currpar = (edges+x) -> value;
				*maxrand = posspar;
			}
		}
		DrawParent(edges, orig, (edges+x)->right, exptimes, inftimes, rectimes, maxrand, currpar, priorparentnode, probmult);
	}

}


/**************************
DrawGraph

Sample the contact network (G) from its full conditional, given all of the other 
parameters (including P).  Also sets the number of I-S edges for each infected, 
and the total number of S-S edges.  Changing G also changes the total infectious 
pressure applied during the epidemic (A), so this is modified as well.
******************************/

void DrawGraph(Network *nwp, Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, double currbeta, double p, int m, double *A)
{
	double u, timeinfected; 						
	Vertex i,j, priorinf, latterinf;
			
	for (j=1; j<=m; j++)	/* Cycle through each dyad - make sure that we consider each dyad (a,b) where a became infectious before b  */
		for(i=1; i<j; i++)
		{			
			priorinf = MININDEX(inftimes,i,j);		/* Set priorinf to the node (i or j) with the earlier infection time, and latterinf to the node with the later infection time */
			latterinf = MAXINDEX(inftimes,i,j);
			
			u = exp(-currbeta*(MAX(MIN(rectimes[priorinf], exptimes[latterinf])  - inftimes[priorinf] , 0)));	/* Construct the probability of an I-I edge existing */
			u = u*p / (1 - p + u*p);

			if ((unif_rand() < u) || (transtree[latterinf] == priorinf))	/* Randomly determine if the edge exists - the edge will exist with prob. u, unless it's part of the trans tree in which case it always exists */
			{
				if ( !EDGE_EXISTS( priorinf, latterinf, nwp ) )
				{
					ToggleEdge( priorinf, latterinf, nwp );	/* Edge doesn't exist, but we need it, so add it */
					*A += MAX(MIN( exptimes[latterinf], rectimes[priorinf] ) - inftimes[priorinf] , 0);	/* Add in the associated infectious pressure */							
				}
			} else
			{
				if ( EDGE_EXISTS(priorinf,latterinf,nwp) )
				{
					ToggleEdge( priorinf, latterinf, nwp );	/* Edge exists, but we don't want it, so delete it */
					*A -= MAX(MIN( exptimes[latterinf], rectimes[priorinf] ) - inftimes[priorinf] , 0);	/* Subtract out the associated infectious pressure */		
				}
			}
		}

	nwp->totaledges = 0;  /* Reset edge count:  Total number of edges = # of I-S edges + # of S-S edges + # of I-I edges */
	
	/* For each infected, set number of edges to susceptibles, i.e., # of I-S edges */
	for (i=1; i<=m; i++) 
	{
		timeinfected = rectimes[i] - inftimes[i]; /* This is the amount of time that Vertex i was infected for */
		*A -= nwp->nisedges[i] * timeinfected;	/* Subtract out old I-S infectious pressure based on old number of I-S edges */
		
		u = exp( -currbeta * timeinfected );	/* Construct the probability of an I-S edge existing */
		u = u * p / ( 1 - p + u*p );
		
		nwp->nisedges[i] = rbinom( nwp->nsus , u );	/* Generate number of I-S edges for (infected) Vertex i  */
		
		*A += nwp->nisedges[i] * timeinfected;	/* Add in new I-S infectious pressure based on new number of I-S edges */
		nwp->totaledges += nwp->nisedges[i];
	}

	nwp->nssedges = rbinom( nwp->nsus * ( nwp->nsus - 1 ) / 2 , p );		/* Get total number of edges between susceptibles */
	
	nwp->totaledges += nwp->nssedges + nwp->niiedges;
}


/***********************
GetRandomOrder

Used to determine the (uniformly random) order in which
to update the exposure/infection times.  Returns a vector with the random 
order of the infecteds.  When updating the exposure times, we exclude
the initial exposed, because it's updated separately 
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


/****************************
LogLikelihood

Calculates the likelihood function, excluding the indicator of whether
or not the transmission tree is legal
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
	double lowexptime;
		
	/* Check if there's only one initial exposed (parent = - 999) */
	for (i=1; i<=m; i++) 	
		if (transtree[i] == -999) 
		{
			j++;
			initialexp = i;
		}	
	if (j > 1) return(0);
	
	/* Make sure initial infected has lowest exposure time */
	lowexptime = exptimes[initialexp];
	for (i=1; i<=m; i++)
		if (i != initialexp)
			if (exptimes[i] < exptimes[initialexp])
				return(0);
			
	/* Check each infection (other than initial) to verify it's indeed possible, given the tree and times */
	for (i=1; i<=m; i++)
		if (i != initialexp)
			if ( (inftimes[transtree[i]] > exptimes[i]) || (exptimes[i] > rectimes[transtree[i]]) || ( !EDGE_EXISTS(transtree[i],i,nwp) ) ) 
				return(0);

	/* Make sure E_j < I_j < R_j for each node j */		
	for (i=1; i<=m; i++)
		if ( (exptimes[i] > inftimes[i]) || (inftimes[i] > rectimes[i]) )
			return(0);		
			
	/* If we've made it to this point, tree is legal */
	return(1);	
}


/**********************
AdjustAiiExpTime

Calculates the change in infectious pressure associated with a 
proposed change in the infection time for one node,due to I-I edges.
Uses a recursive algorithm, and is based on the function InOrderTreeWalk()
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

Calculates the change in infectious pressure associated with a 
proposed change in the infection time for one node. For the change 
due to I-I edges, use a recursive algorithm twice - once for the inedges 
and once for the outedges.  There is no change due to I-S edges.
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

Calculates the change in infectious pressure associated with a 
proposed change in the transition time for one node,due to I-I edges.
Uses a recursive algorithm, and is based on the function InOrderTreeWalk()
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

Calculates the change in infectious pressure associated with a proposed 
change in the transition time for one node. For the change due to I-I edges, 
use a recursive algorithm twice - once for the inedges and once for the outedges.
The change due to I-S edges can be calculated directly.
*****************************/

void AdjustABCInfTime(Network *nwp, Edge orig, double *exptimes,  double *inftimes, double *pinftimes, double *rectimes, 
	double *propA, double *propB, double *propBln, double *propC, double *propCln)
{
	if (((nwp->outedges)+orig)->value != 0) AdjustAiiInfTime(nwp->outedges,orig,orig,exptimes,inftimes,pinftimes,rectimes,propA);
	if (((nwp->inedges)+orig)->value != 0) AdjustAiiInfTime(nwp->inedges,orig,orig,exptimes,inftimes,pinftimes,rectimes,propA);															
	*propA += (inftimes[orig] - pinftimes[orig]) * nwp->nisedges[orig];		/* Adjust A for I-S edges */ 
	
	/* Adjust B and C */
	*propB += pinftimes[orig] - inftimes[orig];
	*propBln += log(pinftimes[orig] - exptimes[orig]) - log(inftimes[orig] - exptimes[orig]);
	*propC += inftimes[orig] - pinftimes[orig];
	*propCln += log(rectimes[orig] - pinftimes[orig]) - log(rectimes[orig] - inftimes[orig]);
}


/*********************
Adjust AClnKappa

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
	*propA += (inftimes[initexp] - pinftimes[initexp]) * nwp->nisedges[initexp];
	
	if (((nwp->outedges)+pinitexp)->value != 0) AdjustAiiInfTime(nwp->outedges,pinitexp,pinitexp,pexptimes,inftimes,pinftimes,rectimes,propA);
	if (((nwp->inedges)+pinitexp)->value != 0) AdjustAiiInfTime(nwp->inedges,pinitexp,pinitexp,pexptimes,inftimes,pinftimes,rectimes,propA);															
	*propA += (inftimes[pinitexp] - pinftimes[pinitexp]) * nwp->nisedges[pinitexp];
	
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
		
	return(logdens);
}



/* **************************** */
/* MAIN MCMC FUNCTION */
/* **************************** */

/*************************
epigraphmcmcc

Main MCMC function to produce a sample of the parameters.  
************************/

void epigraphmcmcc (double *etime, double *itime, double *rtime, int *nsamp, int *thinning,  double *bprior, double *tiprior, double *teprior, double *pprior, double *kiprior,
	double *keprior, int *ninf, int *initN, double *initbeta, double *initthetai, double *initki, double *initthetae, double *initke, double *initp,  int *bpriordist, 
	int *tipriordist, int *tepriordist, int *kipriordist, int *kepriordist, int *accept, int *propose, double *abeta, double *athetai, double *aki, double *athetae, double *ake ,double *ap, 
	int *ainit, double *ainitexptime, double *aexptimes, double *ainftimes, int *atranstree, int *extrathinning, int *inferEtimes, int *inferItimes, int *parentprior, int *probparentmult)
{ 

	GetRNGstate();  /* R function enabling uniform RNG */ 
	
	/* VARIABLE DECLARATIONS */
	
	Network gr, *nwp=&gr;
	int m = *ninf; 							// m = # of infecteds (considered fixed and known); note that this variable holds the same constant as nwp->ninfnodes, just have it here for convenience 
	double llkd, pllkd, currlogprior, proplogprior;					// Variables to hold current and proposed values of the log-likelihoodand priors
	Vertex *transtree = (Vertex *) malloc((m+1) * sizeof(Vertex));	// transtree[i] holds the label of the node that infected node i in the current infection tree...
	Vertex *ptranstree = (Vertex *) malloc((m+1) * sizeof(Vertex));		// ... ptranstree is the proposal version
	Vertex *orderrand = (Vertex *) malloc((m+1) * sizeof(Vertex));		// Update the inf/exp times in random order - orderrand holds this random order
	int iter, ok = 1, zz, i, j, currcount, propcount;					// dummy & counter variables
	double currbeta = *initbeta, currthetai = *initthetai, currthetae = *initthetae; // Starting values for the parameters...
	double p = *initp, currki = *initki, currke = *initke;	//... starting values for more parameters
	double A, propA, B = 0, propB, Bln = 0, propBln, C = 0, propC, Cln = 0, propCln;				// Variables to keep track of infectious, transition and removal pressure
	double *exptimes = (double *) malloc((m+1) * sizeof(double));				// Use indices 1..m and leave 0 unused to be consistent with Network object and R
	double *pexptimes = (double *) malloc((m+1) * sizeof(double));			// These hold the current and proposed exposure times
	double *inftimes = (double *) malloc((m+1) * sizeof(double));				
	double *pinftimes = (double *) malloc((m+1) * sizeof(double));			// These hold the current and proposed infection times
	double *rectimes = (double *) malloc((m+1) * sizeof(double));			// These hold the recovery times
	int *probparentprior = (int *) malloc((m+1) * sizeof(int));			// These hold the prior belief about most likely parent nodes
	Vertex initexp, pinitexp;						// These hold the current and proposed initial exposeds (Kappa)	
	double propthetai, propthetae, propbeta, propki, propke;	// Variables to hold proposal valuesfor parameters
	double delta = 0.5, uniwidth = 5;	// Tuning parameters - maybe will pass them in or figure out a way to automatically set them - these values seem to work reasonably well...
	double bwidth = (bprior[1] - bprior[0])/uniwidth, tiwidth = (tiprior[1] - tiprior[0])/uniwidth, tewidth = (teprior[1] - teprior[0])/uniwidth;
	double kiwidth = (kiprior[1] - kiprior[0])/uniwidth, kewidth = (keprior[1] - keprior[0])/uniwidth;	// Variables used in case of uniform priors
	int maxmove = 11;	// Sets the number of updates in each sweep of the algorithm
	
	/* INITIALIZE VARIABLES */

	for (i = 0; i < m; i++) 
	{
		exptimes[i+1] = etime[i];		/* Copy initial exposure times(since they might change during the procedure)	*/
		inftimes[i+1] = itime[i];	/* Copy initial infection times(since they might change during the procedure)	*/
		rectimes[i+1] = rtime[i];	/* Copy recovery times just for indexing consistency  */
		probparentprior[i+1] = parentprior[i];	/* Copy prior beliefs about parents just for indexing consistency  */
	}
		
 	if (InitializeTransTree(transtree, exptimes, inftimes, rectimes, m, &initexp, &A) != 1) 	/* Initializes transmission tree, determines the initial infected, and calculates the amount of  infection pressure (A) due to the transmission tree */
	{
		Rprintf("Faulty exposure/infection/recovery time data.  Cannot build legal initial transmission tree.  Aborting routine. \n");
		return;
	}
	
	gr = InitializeNetworkFromTree(transtree, m, *initN-m, *initN);	/* Initialize network to only contain the edges from the transmission tree */		
	
	if (!IsTreeLegal(exptimes,inftimes,rectimes,transtree,nwp,m)) 		/* If initial transmission tree isn't legal, then the set of infection times and recovery times are somehow faulty, so we print an error and abort */
	{
		Rprintf("Faulty exposure/infection/recovery time data.  Cannot build legal initial transmission tree.  Aborting routine. \n");
		return;
	}

	for (i=1; i <= m; i++)	/* Calculate initial transition and removal pressures and their "log" versions */
	{
		C += rectimes[i] - inftimes[i];
		Cln += log(rectimes[i] - inftimes[i]);
		B += inftimes[i] - exptimes[i];
		Bln += log(inftimes[i] - exptimes[i]);
	}
	
	DrawGraph(nwp, transtree, exptimes, inftimes, rectimes, currbeta, p, m, &A);  /* Draw initial graph from its full conditional dist. */
		
	/* BEGIN MAIN MCMC LOOP */
	
	for (iter = 0; iter < ( (*nsamp) * (maxmove) ); iter++)
	{
		
		zz = iter % (maxmove);		/* Cyclical updates - can change this line to update different sets of parameters or change the order of updates */		
		if (zz <= 7) propose[zz] ++;		/* Not currently keeping track of propose / accept info for E/I times or Kappa */
		switch(zz)		/* Choose which item (parameter) to update */
		{
			
		case 0:	/* Update Transmission Tree (P)	*/
			
			DrawTransTree(transtree, nwp, exptimes, inftimes, rectimes, m, probparentprior, probparentmult);	/* transtree will have the updated transmission tree */
			accept[zz]++;
			break;

		case 1:	/* Update p */
						
			p = rbeta( nwp->totaledges + pprior[0], ( nwp->N ) * ( nwp->N - 1 ) / 2 - nwp->totaledges + pprior[1] ); /* Beta is the conjugate prior for p - use Gibbs sampler to update p */
			accept[zz]++;
			break;
			
		case 2:	/* Update graph (G) */
		
			DrawGraph(nwp, transtree, exptimes, inftimes, rectimes, currbeta, p, m,  &A);
			accept[zz]++;
			break;						

		case 3:	/* Update Beta	*/
			
			if (*bpriordist == 0)	/* Uniform prior for Beta */
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
			} else if(*bpriordist == 1)	/* Gamma (conjugate) prior for Beta */
			{
				currbeta = rgamma( bprior[0] + m - 1, 1 / ( ( 1 / bprior[1] ) + A ) );
				accept[zz]++;
			}	
			break;
						
		case 4:	/* Update Theta_I */
		
			if (*tipriordist == 0)	/* Uniform prior for Theta_I */
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
			} else if(*tipriordist == 1)	/* Inverse Gamma (conjugate) prior for Theta_I */
			{
				currthetai = 1/rgamma( tiprior[0] + currki*m, 1/ (tiprior[1] + C) );			
				accept[zz]++;
			} 	
			break;
			
		case 5:  /* Update k_I parameter for removal process */	
			
			ok = 1;
			propki = currki + unif_rand()*kiwidth - kiwidth/2;

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

			if (*tepriordist == 0)	/* Uniform prior for Theta_E */
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
			} else if(*tepriordist == 1)	/* Inverse Gamma (conjugate) prior for Theta_E */
			{
				currthetae = 1/rgamma( teprior[0] + currke*m, 1/ (teprior[1] + B) );			
				accept[zz]++;
			} 	
			
			break;
			
		case 7:		/* Update k_E */

			ok = 1;
			propke = currke + unif_rand()*kewidth - kewidth/2;

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
			
		case 8:		/* Update infection times */

			if (*inferItimes == 0) break;			/* If the infection times are known and fixed, then there's nothing to do here. */
			
			for (i=1; i <=m; i++) pinftimes[i] = inftimes[i];		/* Initialize proposed infection times - will start with current infection times and update one at a time*/
			
			llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* Calculate current log likelihood */

			GetRandomOrder(orderrand, initexp, 1, m);	/* Orderrand contains the (random) order to update infection times  */
			
			for (i=1; i <= m; i++)	/* Cycle through the nodes  */
			{
				ok = 1;
				j = orderrand[i];	/* Update infection time for Vertex j */
				
				pinftimes[j] = ProposedInftime(j,transtree, exptimes, inftimes, rectimes, m);	/* Get proposed infection time for vertex j */
				if (IsTreeLegal(exptimes,pinftimes,rectimes,transtree,nwp,m) == 0) ok = 0;	/* If proposed times don't result in a legal tree, no need to calculate the rest of the likelihood */
								
				if (ok == 1) 		/* We've proposed a legal infection time, now check to see if we accept it */
				{
					propA = A; propB = B; propBln = Bln; propC = C; propCln = Cln;
					AdjustABCInfTime(nwp, j, exptimes, inftimes, pinftimes, rectimes, &propA, &propB, &propBln, &propC, &propCln);
					pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, propA, propB, propBln, propC, propCln);
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

			for (i=1; i <=m; i++) pexptimes[i] = exptimes[i];		/* Initialize proposed exposure times - will start with current exposure times and update one at a time*/
			
			llkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, A, B, Bln, C, Cln); 	/* Calculate current log likelihood */

			GetRandomOrder(orderrand, initexp, 0, m);	/* Orderrand contains the (random) order to update exposure times except for initial exposed */
			
			for (i=1; i<m; i++)	/* Cycle through the exposeds, other than the initial exposed */
			{
				ok = 1;
				j = orderrand[i];	/* Update exposure time for Vertex j */
				pexptimes[j] = ProposedExptime(j,transtree,inftimes);	/* Get proposed exposure time for vertex j */
				if (IsTreeLegal(pexptimes,inftimes,rectimes,transtree,nwp,m) == 0) ok = 0;	/* If proposed times don't result in a legal tree, no need to calculate the rest of the likelihood */
				
				if (ok == 1) 		/* We've proposed a legal time, now check to see if we accept it */
				{
					propA = A; propB = B; propBln = Bln; propC = C; propCln = Cln;
					AdjustABExpTime(nwp, j, exptimes, pexptimes, inftimes, rectimes, &propA, &propB, &propBln);
					pllkd = LogLikelihood(currbeta,currthetai, currki, currthetae, currke, m, propA, propB, propBln, propC, propCln);
					if (log(unif_rand()) > (pllkd - llkd) ) ok = 0;
				}
				
				if (ok == 1) 		/* Accepted */
				{		
					A = propA; B = propB; Bln = propBln; C = propC; Cln = propCln;
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
			
			ptranstree[pinitexp] = -999;			/* Swap proposal direction of infection  */
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
		
		if (iter % ( (maxmove) * (*thinning) ) == 0)		
		{  
			j = floor(iter / ( (maxmove) * (*thinning) ) );				/* 	Every thinning sweeps (updates of each parameter), record parameter values */
			athetai[j] = currthetai;				/*	These are storage vectors that are initialized in R and passed (empty) into C	*/		
			abeta[j] = currbeta;					/*	They're passed back to R containing the sample values */
			aki[j] = currki;
			athetae[j] = currthetae;
			ake[j] = currke;
			ap[j] = p;
			ainitexptime[j] = exptimes[initexp];
			ainit[j] = initexp;
			
			R_CheckUserInterrupt();			/* Occasionally check for an R user interrupt */
			
			if (*extrathinning > 0)			/* Record exposure and infection time, as well as the transmission tree, if desired */
				if  (j % (*extrathinning) == 0)				
				{
					j = floor(j / (*extrathinning) );
					for (i=1; i<=m; i++)
					{
						aexptimes[j*m+i - 1] = exptimes[i];
						ainftimes[j*m+i - 1] = inftimes[i];
						atranstree[j*m+i - 1] = transtree[i];
					}
				}
		}

	}	/* End MCMC loop */
			
	/* Return memory used before leaving */
	NetworkDestroy(nwp);
	free(transtree);
	free(ptranstree);
	free(exptimes);
	free(pexptimes);
	free(rectimes);
	
	PutRNGstate(); /* Must be called after GetRNGstate before returning to R */

}