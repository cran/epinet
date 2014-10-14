
/*
 *  File NetworkEpiMCMCSEIR.h
 *
 *  Part of the `epinet' R package.
 */

 
#ifndef NETWORKEPIMCMCSEIR_H
#define NETWORKEPIMCMCSEIR_H

#include "NetworkFunctions.h" 

/* 
Some Definitions:

A is the total time that susceptibles are exposed to infecteds, i.e., total number 
of person-to-person units of infectious pressure applied during the epidemic 

B is the total time spent (by all individuals) in the exposed state, i.e., total number of
units of transition (from Exposed to Infective) pressure applied during the epidemic.

C is the total time spent (by all individuals) in the infected state, i.e., the total number 
of units of removal pressure applied during the epidemic.

All of these "pressures" (as well as the log-versions of B and C) are maintained throughout
the MCMC routine; whenever any variable changes that would affect these pressures, their values are
updated.  This is done using the "Adjustxxx()" functions.
 
transtree is a vector holding the (current) transmission tree.  transtree[i] contains
the identity of the parent of i.  transtree[0] isn't used and is set to 0.  The node j with 
transtree[j] = -999 is the initial exposed.  Any node k with a value of transtree[k] = -77
is a susceptible, i.e., an individual that was never infected during the course of the
epidemic, and hence has no parent in the transmission tree
*/


/* Initialization Functions */
int InitializeTransTree(Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, int m, int N, Vertex *initexp, double *A);
Network InitializeNetworkFromTree(Vertex *transtree, int m, Vertex N);

/* Functions to Generate Proposal Values */ 
Vertex ProposedInitExp(Vertex initexp, Vertex *transtree, int m);	/* Generate proposed identity of initial exposed */
double ProposedInitialExptime(Vertex init, double *inftimes, double delta); /* Generate proposed exposure time for initial exposed */
double ProposedExptime(Vertex j,Vertex *transtree,double *inftimes); /* Generate proposed exposure time for infecteds other than the initial exposed */
double ProposedInftime(Vertex j, Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, int m);	/* Generate proposed transition time for one individual */
void DrawTransTree(Vertex *transtree, Network *nwp, double *exptimes, double *inftimes, double *rectimes, int m, int *probparentprior, int *probparentmult); /* Sample a transmission tree (P) from its full conditional */
void DrawParent(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *rectimes, double *maxrand, Vertex *currpar, int priorparentnode, int probmult);
double CalcEdgeProb(int dyadcovindex, double *dyadcovs, double *eta, int etapars, int N);
void DrawGraph(Network *nwp, Vertex *transtree, double *exptimes, double *inftimes, double *rectimes, double currbeta, double *dyadcovs, int *dyadindex1, int *dyadindex2, double *eta, int etapars, double *A); /* Sample I-I network from its full conditional */

/* Utility Functions */
void GetRandomOrder(Vertex *order, Vertex initexp, int includeinit, int m);
int GetDyadIndex(int nodeindex1, int nodeindex2, int N);
void CreateRandomDyadOrder(int *dyadindex1, int *dyadindex2, int N);
void AdjustAiiExpTime(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *pexptimes, double *inftimes, double *rectimes, double *propA);
void AdjustABExpTime(Network *nwp, Edge orig, double *exptimes, double *pexptimes, double *inftimes, double *rectimes, double *propA, double *propB, double *propBln);
void AdjustAiiInfTime(TreeNode *edges, Edge orig, Edge x, double *exptimes, double *inftimes, double *pinftimes, double *rectimes, double *propA);
void AdjustABCInfTime(Network *nwp, Edge orig, double *exptimes,  double *inftimes, double *pinftimes, double *rectimes, double *propA, double *propB, double *propBln, double *propC, double *propCln);
void AdjustAClnKappa(Network *nwp, Edge initexp, Edge pinitexp, double *exptimes,  double *pexptimes, double *inftimes, double *pinftimes, double *rectimes, double *propA, double *propCln);

/* Likelihood, Priors and Related Functions */
double LogLikelihood(double beta, double thetai, double ki ,double thetae, double ke, int m, double A, double B, double Bln, double C, double Cln);
int IsTreeLegal(double *exptimes, double *inftime, double *rectimes, int *transtree, Network *nwp, int m);
double LogkPrior(double x, double *kprior, int *kpriordist);
void LogEdgeCalc(TreeNode *edges, int orig, int x, double *eta, int etapars, double *dyadcovs, int totaldyads, double *logdensptr, int N);
double LogGraphCalc(Network *nwp, double *eta, int etapars, double *dyadcovs);

/* MAIN MCMC FUNCTION - This function is accessed from the R wrapper function ... */
void epigraphmcmcc (double *etime, double *itime, double *rtime, int *etapars, double *dyadcovs, int *nsamp, int *thinning,  double *bprior, double *tiprior, double *teprior, double *etaprior, double *kiprior,
	double *keprior, int *ninf, int *initN, double *initbeta, double *initthetai, double *initki, double *initthetae, double *initke, double *initeta,  int *bpriordist, 
	int *tipriordist, int *tepriordist, int *kipriordist, int *kepriordist, int *etapriordist, double *etapropsd, int *accept, int *propose, double *allkd, double *abeta, double *athetai, double *aki, double *athetae, double *ake ,double *aeta, 
	int *ainit, double *ainitexptime, double *aexptimes, double *ainftimes, int *atranstree, int *extrathinning, int *inferEtimes, int *inferItimes, int *parentprior, int *probparentmult,
	int *verbose); 

#endif
