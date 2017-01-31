/*
 * The routine functions for 2-D Barnes-Hut Trees
 *
 * Copyright (c) 2014, Zhirong Yang (Aalto University)
 * All rights reserved.
 *
 */

#ifndef BARNES_HUT
#define BARNES_HUT

typedef struct tagNode {
    int id;
    double weight;
} Node;

typedef struct tagQuadTree{
    Node *node;
    struct tagQuadTree **children;
    int childrenLength;
    int childCount;
    double position[2];
    double weight;
    double coef; /*reserved for application use*/
    double extcoef[4]; /*reserved for application extended use*/
    double minPos[2];
    double maxPos[2];
} QuadTree;

/* externally callable functions */
QuadTree* buildQuadTree(double *Y, double *weights, int n); /* weights==NULL would use 1.0 for each node weight*/
void destroyQuadTree(QuadTree *tree);
double getTreeWidth(QuadTree *tree);

/* externally callable debugging utilities */
void dumpTree(QuadTree *tree, int depth);
int getNumberOfBranchNodes(QuadTree *tree);
void getBoundaries(QuadTree *tree, double *bnds, int *pind, int nn);
void getTreeNodeInfo(QuadTree *tree, double *bnds, double *pos, double *depths, int *pind, int nn, int depth);

/* internal functions */
QuadTree* createQuadTree(Node *node, double position[], double minPos[], double maxPos[]);
void addNode(QuadTree *tree, Node *newNode, double newPos[], int depth);
void addNode2(QuadTree *tree, Node *newNode, double newPos[], int depth);

#endif