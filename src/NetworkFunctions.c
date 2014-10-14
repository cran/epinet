/*
 *  File NetworkFunctions.c
 *  Part of the `epinet' R package
 *
 *  Based on 'statnet' project software (http://statnetproject.org). 
 *  For license and citation information see http://statnetproject.org/attribution 
 *  
 *  In particular, based on code from the "ergm" R package, modified
 *  for use in the epinet package.  
 *
 *  This file contains the code to manipulate the Network structure.
 *  This structure is used to hold the information about the
 *  relationships between members of a population across which
 *  an epidemic has spread.  This version implements the "Group Model",
 *  so that it also contains some information about the number of
 *  edges and dyads within each group.
 * 
 */

#include "NetworkFunctions.h"
#include "NetworkEpiMCMCSEIR.h"

/*******************
 Network NetworkInitialize

 Initialize, construct binary tree version of network.  Note
 that the 0th TreeNode in the array is unused and should 
 have all its values set to zero.
*******************/

Network NetworkInitialize(Vertex *heads, Vertex *tails, Edge edgecount, 
			  Vertex ninfnodes, Vertex N) {
  Network nw;
	int i;

  nw.next_inedge = nw.next_outedge = (Edge)N+1;
  nw.outdegree = (Vertex *) calloc((N+1),sizeof(Vertex));
  nw.indegree  = (Vertex *) calloc((N+1),sizeof(Vertex));
  nw.maxedges = MAX(edgecount,1)+N+2; /* Maybe larger than needed? */
  nw.inedges = (TreeNode *) calloc(nw.maxedges,sizeof(TreeNode));
  nw.outedges = (TreeNode *) calloc(nw.maxedges,sizeof(TreeNode));

  /*Configure a Network*/

  /* Will add edges one-by-one, so leave nw.numedges as all zeros for now */
	  
  nw.ninfnodes = ninfnodes;
  nw.N = N;

  ShuffleEdges(heads,tails,edgecount); /* shuffle to avoid worst-case performance */
	
  for(i = 0; i < edgecount; i++) {
    Vertex h=heads[i], t=tails[i];
    if (h > t) 
      AddEdgeToTrees(t,h,&nw); /* Undir edges always have head < tail */ 
    else 
      AddEdgeToTrees(h,t,&nw);
  }  
  return nw;
}

/*******************
 void NetworkDestroy
*******************/
void NetworkDestroy(Network *nwp) {
  free (nwp->indegree);
  free (nwp->outdegree);
  free (nwp->inedges);
  free (nwp->outedges);
}

/*****************
 Edge EdgetreeSearch

 Check to see if there's a TreeNode with value b 
 in the tree rooted at edges[a].  Return i such that 
 edges[i] is that TreeNode, or 0 if none.
*****************/
Edge EdgetreeSearch (Vertex a, Vertex b, TreeNode *edges) {
  TreeNode *es;
  Edge e = a;
  Vertex v;

  es = edges + e;
  v = es->value;
  while(e != 0 && b != v)
    {
      e = (b<v)?  es->left : es->right;
      es = edges + e;
      v = es->value;
    }
  return e;
}

/*****************
 Edge EdgetreeSuccessor

 Return the index of the TreeNode with the smallest value
 greater than edges[x].value in the same edge tree, or 0
 if none.  This is used by (for instance)
 the DeleteHalfedgeFromTree function.
*****************/
Edge EdgetreeSuccessor (TreeNode *edges, Edge x) {
  TreeNode *ptr;
  Edge y;

  if ((y=(ptr=edges+x)->right) != 0) 
    return EdgetreeMinimum (edges, y);
  while ((y=ptr->parent)!=0 && x==(ptr=edges+y)->right) 
    x=y;
  return y; 
}   

/*****************
 Edge EdgetreeMinimum

 Return the index of the TreeNode with the
 smallest value in the subtree rooted at x
*****************/
Edge EdgetreeMinimum (TreeNode *edges, Edge x) {
  Edge y;

  while ((y=(edges+x)->left) != 0)
    x=y;
  return x;
}

/*****************
 Edge ToggleEdge

 Toggle an edge:  Set it to the opposite of its current
 value.  Return 1 if edge added, 0 if deleted.
*****************/
int ToggleEdge (Vertex head, Vertex tail, Network *nwp) 
{
  if (head > tail) {
    Vertex temp;
    temp = head; /*  Make sure head<tail always for undirected edges */
    head = tail;
    tail = temp;
  }
  if (AddEdgeToTrees(head,tail,nwp))
    return 1;
  else 
    return 1 - DeleteEdgeFromTrees(head,tail,nwp);
}

/*****************
 Edge AddEdgeToTrees

 Add an edge from head to tail after checking to see
 if it's legal. Return 1 if edge added, 0 otherwise.  Since each
 "edge" should be added to both the list of outedges and the list of 
 inedges, this actually involves two calls to AddHalfedgeToTree (hence
 "Trees" instead of "Tree" in the name of this function).  Also increments
 the number of edges in the group corresponding to the given dyad.
*****************/
int AddEdgeToTrees(Vertex head, Vertex tail, Network *nwp){
  if (EdgetreeSearch(head, tail, nwp->outedges) == 0) {
    AddHalfedgeToTree(head, tail, nwp->outedges, nwp->next_outedge);
    AddHalfedgeToTree(tail, head, nwp->inedges, nwp->next_inedge);
    ++nwp->outdegree[head];
    ++nwp->indegree[tail];
    UpdateNextedge (nwp->inedges, &(nwp->next_inedge), nwp); 
    UpdateNextedge (nwp->outedges, &(nwp->next_outedge), nwp);
    return 1;
  }
  return 0;
}

/*****************
 void AddHalfedgeToTree:  Only called by AddEdgeToTrees
*****************/
void AddHalfedgeToTree (Vertex a, Vertex b, TreeNode *edges, Edge next_edge){
  TreeNode *eptr = edges+a, *newnode;
  Edge e;
  if (eptr->value==0) { /* This is the first edge for vertex a. */
    eptr->value=b;
    return;
  }
  (newnode = edges + next_edge)->value=b;  
  newnode->left = newnode->right = 0;
  /* Now find the parent of this new edge */
  for (e=a; e!=0; e=(b < (eptr=edges+e)->value) ? eptr->left : eptr->right);
  newnode->parent=eptr-edges;  /* Point from the new edge to the parent... */
  if (b < eptr->value)  /* ...and have the parent point back. */
    eptr->left=next_edge; 
  else
    eptr->right=next_edge;
}

/*****************
void UpdateNextedge
*****************/
void UpdateNextedge (TreeNode *edges, Edge *nextedge, Network *nwp) {
  int mult=2;
  /*TreeNode *tmp_in, *tmp_out; */
  
  while (++*nextedge < nwp->maxedges) {
    if (edges[*nextedge].value==0) return;
  }
  /* Reached end of allocated memory;  back to start and recheck for "holes" */
  for (*nextedge = (Edge)nwp->N+1; *nextedge < nwp->maxedges; ++*nextedge) {
    if (edges[*nextedge].value==0) return;
  }
  /* There are no "holes" left, so this network overflows mem allocation */
  nwp->maxedges *= mult;
  nwp->inedges = (TreeNode *) realloc(nwp->inedges, 
                                      sizeof(TreeNode) * nwp->maxedges);
  memset(nwp->inedges+*nextedge,0,sizeof(TreeNode) * (nwp->maxedges-*nextedge));
  nwp->outedges = (TreeNode *) realloc(nwp->outedges, 
                                       sizeof(TreeNode) * nwp->maxedges);
  memset(nwp->outedges+*nextedge,0,sizeof(TreeNode) * (nwp->maxedges-*nextedge));
}

/*****************
 int DeleteEdgeFromTrees

 Find and delete the edge from head to tail.  
 Return 1 if successful, 0 otherwise.  As with AddEdgeToTrees, this must
 be done once for outedges and once for inedges.  Decrements the number 
 of edges in the given group.
*****************/
int DeleteEdgeFromTrees(Vertex head, Vertex tail, Network *nwp){
  if (DeleteHalfedgeFromTree(head, tail, nwp->outedges,&(nwp->next_outedge))&&
      DeleteHalfedgeFromTree(tail, head, nwp->inedges, &(nwp->next_inedge))) {
    --nwp->outdegree[head];
    --nwp->indegree[tail];
    return 1;
  }
  return 0;
}

/*****************
 int DeleteHalfedgeFromTree

 Delete the TreeNode with value b from the tree rooted at edges[a].
 Return 0 if no such TreeNode exists, 1 otherwise.  Also update the
 value of *next_edge appropriately.
*****************/
int DeleteHalfedgeFromTree(Vertex a, Vertex b, TreeNode *edges,
		     Edge *next_edge){
  Edge x, z, root=(Edge)a;
  TreeNode *xptr, *zptr, *ptr;

  if ((z=EdgetreeSearch(a, b, edges))==0)  /* z is the current TreeNode. */
    return 0;  /* This edge doesn't exist, so return 0 */
  /* First, determine which node to splice out; this is z.  If the current
     z has two children, then we'll actually splice out its successor. */
  if ((zptr=edges+z)->left != 0 && zptr->right != 0) {
    z=EdgetreeSuccessor(edges, z);  
    zptr->value = (ptr=edges+z)->value;
    zptr=ptr;
  }
  /* Set x to the child of z (there is at most one). */
  if ((x=zptr->left) == 0)
    x = zptr->right;
  /* Splice out node z */
  if (z == root) {
    zptr->value = (xptr=edges+x)->value;
    if (x != 0) {
      if ((zptr->left=xptr->left) != 0)
	(edges+zptr->left)->parent = z;
      if ((zptr->right=xptr->right) != 0)
	(edges+zptr->right)->parent = z;
      zptr=edges+(z=x);
    }  else 
      return 1;
  } else {
    if (x != 0)
      (xptr=edges+x)->parent = zptr->parent;
    if (z==(ptr=(edges+zptr->parent))->left)
      ptr->left = x;
    else 
      ptr->right = x;
  }  
  /* Clear z node, update *next_edge if necessary. */
  zptr->value=0;
  if (z < *next_edge)
    *next_edge=z;
  return 1;
}

/*****************
 void printedge

 Diagnostic routine that prints out the contents
 of the specified TreeNode (used for debugging).  
*****************/
void printedge(Edge e, TreeNode *edges){
  Rprintf("Edge structure [%d]:\n",e);
  Rprintf("\t.value=%d\n",edges[e].value);
  Rprintf("\t.parent=%d\n",edges[e].parent);
  Rprintf("\t.left=%d\n",edges[e].left);
  Rprintf("\t.right=%d\n",edges[e].right);
}

/*****************
 void NetworkEdgeList

 Print an edgelist for a network
*****************/
void NetworkEdgeList(Network *nwp) {
  Vertex i;
  for(i=1; i<=(nwp->N); i++) {
    Rprintf("Node %d:\n  ", i);
    InOrderTreeWalk(nwp->outedges, i);
    Rprintf("\n");
  }
}

/*****************
 void InOrderTreeWalk

 Diagnostic routine that prints the nodes in the tree rooted
 at edges[x], in increasing order of their values.
*****************/
void InOrderTreeWalk(TreeNode *edges, Edge x) {
  if (x != 0) {
    InOrderTreeWalk(edges, (edges+x)->left);
    //printedge(x, edges);
    Rprintf(" %d ",(edges+x)->value); 
    InOrderTreeWalk(edges, (edges+x)->right);
  }
}

/*****************
 void ShuffleEdges
 
 Knuth shuffle to give a random permutation of the
 list of edges provided.  This is done to avoid
 worst-case performance.  Done prior to adding these
 edges to the trees when building (initializing) the network.
*****************/

void ShuffleEdges(Vertex *heads, Vertex *tails, Edge edgecount){
  for(Edge i = edgecount; i > 0; i--) {
    Edge j = (double) i * unif_rand();
    Vertex h = heads[j];
    Vertex t = tails[j];
    heads[j] = heads[i-1];
    tails[j] = tails[i-1];
    heads[i-1] = h;
    tails[i-1] = t;
  }
}
