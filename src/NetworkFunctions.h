/*
 *  File NetworkFunctions.h
 *  Part of the `epinet' R package
 *
 *  Based on 'statnet' project software (http://statnetproject.org). 
 *  For license and citation information see http://statnetproject.org/attribution 
 *  
 *  In particular, based on code from the "ergm" R package
 * 
 *  This file contains the definition of the Network structure, and
 *  also the associated functions to manipulate this structure.  This
 *  structure is a modified version of the Network structure from the
 *  ergm package (part of the 'statnet' project software), modified for
 *  the epinet package.
 */

#ifndef NETWORKFUNCTIONS_H
#define NETWORKFUNCTIONS_H

#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#define MIN(a,b) ((a)<(b) ? (a) : (b))
#define MAX(a,b) ((a)<(b) ? (b) : (a))
#define MININDEX(x,a,b) ((x[(a)])<(x[(b)]) ? (a) : (b))
#define MAXINDEX(x,a,b) ((x[(a)])<(x[(b)]) ? (b) : (a))
#define EDGE_EXISTS(a,b,nwp) (EdgetreeSearch(MIN(a,b),MAX(a,b),nwp->outedges)!=0?1:0)

typedef int Vertex;
typedef int Edge;

/*  TreeNode is a binary tree structure, which is how the edgelists 
    are stored.  The root of the tree for vertex i will be inedges[i]
    or outedges[i].  inedges[0] and outedges[0] are unused, since the
    index 0 will indicate no link.  
*/

typedef struct TreeNodestruct {
  Vertex value;      /*  the vertex at the other end of this edge  */
  Edge parent;   /*  parent of this node in the tree (0 for root) */
  Edge left;     /*  left child (0 if none)  */
  Edge right;    /*  right child (0 if none) */
} TreeNode;


/* Network is a structure containing all essential elements
   of a given network; it is based on the Network data
   structure from the "ergm" package (see above), with 
   deletion of unnecessary things, and a few new fields 
   that are specific to this application.  Unlike previous
   versions, we are storing ALL individuals in the network,
   rather than just the infecteds.  By convention, we assume
   that the individuals labelled 1..m (i.e., 1..ninfnodes) are
   infected, and the rest are susceptible (never infected)

   Some of the fields in a Network structure are:
   inedges and outedges are arrays of TreeNode that are used to 
     store all of the incoming and outgoing edges, respectively. 
   next_inedge and next_outedge are continually updated to give
     the smallest index of an edge object not being used.  
   outdegree[] and indegree[] are continually updated to give
     the appropriate degree values for each vertex.  These should
     point to Vertex-vectors of length ninfnodes.  
   ninfnodes is the number of infected individuals (also referred
    to sometimes as "m")
   N is the total number of members in the population (infecteds
	+ susceptibles)
   numgroups holds the number of dyadic groups; this is defined
	by the data and should not change
   numdyads holds the number of dyads in each dyadic group; this
	is defined by the data and should not change
   numedges keeps the number of edges (ties) in each dyadic group;
	this is updated as the graph changes
 */
     
typedef struct Networkstruct {
  TreeNode *inedges;
  TreeNode *outedges;
  Vertex ninfnodes;
  Vertex N;
  Edge next_inedge;
  Edge next_outedge;
  Vertex *indegree;
  Vertex *outdegree;
  Edge maxedges;
} Network;


/* Initialization and destruction. */
Network NetworkInitialize(Vertex *heads, Vertex *tails, Edge edgecount, Vertex ninfnodes, Vertex N);
void NetworkDestroy(Network *nwp);

/* Accessors. */
Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges);
Edge EdgetreeSuccessor (TreeNode *edges, Edge x);
Edge EdgetreeMinimum (TreeNode *edges, Edge x);

/* Modifiers. */
int ToggleEdge (Vertex head, Vertex tail, Network *nwp);
int AddEdgeToTrees(Vertex head, Vertex tail, Network *nwp);
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge next_edge);
void UpdateNextedge (TreeNode *edges, Edge *nextedge, Network *nwp);
int DeleteEdgeFromTrees(Vertex head, Vertex tail, Network *nwp);
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges, Edge *next_edge);

/* Utility functions. */
void printedge(Edge e, TreeNode *edges);
void InOrderTreeWalk(TreeNode *edges, Edge x);
void NetworkEdgeList(Network *nwp);
void ShuffleEdges(Vertex *heads, Vertex *tails, Edge edgecount);

#endif
