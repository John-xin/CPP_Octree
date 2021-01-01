#ifndef COCTREE_H
#define COCTREE_H

// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>

#include <utility>   // pair

#include "cFeature.h"
#include "cGeomData.h"
#include "cOctNode.h"


// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cOctree {
public:
    static const int MAX_OCTREE_LEVELS = 4;
    static const int MIN_OCTREE_LEVELS = 3;
    static const int MAX_OCTNODE_FEATS = 1;
    int branchOffsets[8][3];
    int depth;
    cOctNode root;
    vector<cOctNode*> allNodesList;
    vector<cOctNode*> bNodesList;
    vector<cOctNode*> nonBNodesList;
    vector<cOctNode*> leafNodesList;
    vector<cOctNode*> tmpNodesList;
    int num_allNodes, num_intlNodes, num_extNodes, num_bNodes;

    //int numOfNodesTested=0;

    cOctree();
    ~cOctree();
    int countAllNodes(cOctNode *node);
    int countIntlNodes(cOctNode *node);
    int countExtNodes(cOctNode *node);
    int countBNodes(cOctNode *node);
    void showOctreeNodes();


    void setup_root();
    double getSizeRoot();//get length of root cube
    vector<double> getPositionRoot();//find low/upp of geom -> get root cube center position
    vector<double> setCustomPosition(double x,double y,double z);
    vector<cOctNode*> getLeafNodeByPt(vector<double> pt, cOctNode* node);
    void addLeafNodeByPt(vector<double> pt,cOctNode* node);
    cOctNode* getNodeFromId(string nodeId);
    cOctNode* findBranchById(string nodeId, cOctNode *node);
    //void findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList);

    vector<Intersection> findRayIntersect(cLine &ray);
    set<int> getListPolysToCheck(cLine &ray);
    void getPolysToCheck(cOctNode *node, cLine &ray, set<int> &intTestPolys);
    //void getNodesToCheck(cOctNode &node, cLine &ray, vector<pair<cOctNode*,double> > &nodeList);
    vector<int> findRayIntersects(vector<cLine> &rayList);
    vector<int> findRayIntersectsSorted(vector<cLine> &rayList);		
    vector<cOctNode*> getNodesFromLabel(int polyLabel);	
    vector<cOctNode*> getSortedNodesToCheck(cLine &ray);

    void splitNode(cOctNode* node);
    void splitNodeByLevel(int level, cOctNode *node);
    void splitOctreeByMinLevel(cOctNode *node);
    void splitNodeByFeaturePt(cOctNode *node);
    void splitOctreeByFeaturePt(cOctNode *node);
    void getOctreeDepth(cOctNode *node);
    void balanceOctree(cOctNode *node);
    void splitNodeByLevelDiff(int lvlDiff,cOctNode *node);
    void splitNodeById(string node_id);
    void splitExtNodeNbr(cOctNode* node);
    void update_allNodesList(cOctNode* node);
    void splitNodeByPhyName(string phyName, int level, cOctNode* node);
    void setup_boundaryNode(cOctNode* node);
    void setup_interiorNode(cOctNode* node);
    int isInteriorNode(cOctNode* node);
    void setup_nbrNodesState(cOctNode* node);

    void setup_leafNodesList(cOctNode* node);
    void setup_leafNodesNbr();
    void setLeafNodeNbr(cOctNode* node);
    void delExtNodes();
    void delExtNodes2(cOctNode* &node);

    void addGapNodes();

};

//vector<cTri> geoFFacesList; //geom feature faces list
//vector<cFeaturePt> geoFPtsList; //geom feature pts list
//vector<cFeatureEdge> geoFEdgesList; //geom feature edges list

#endif
