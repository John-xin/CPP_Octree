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
#include "cFace.h"
#include "cBoundary.h"
#include "cTri.h"
#include "cLine.h"
#include "cGeomData.h"
#include "cOctNode.h"


// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cOctree {
public:

    static const int MAX_OCTREE_LEVELS = 3;
    static const int MIN_OCTREE_LEVELS = 1;
    static const int MAX_OCTNODE_FEATS = 2;
    int branchOffsets[8][3];
    int depth;
    cOctNode root;

    //octree setup depends on geomData - run geomData read STL first
    //simple geom data
    vector<vector<double> >& geoPts3DList=cGeomData::getInstance()->pts3DList;
    vector<vector<int> >& geoTrisList=cGeomData::getInstance()->trisList;
    vector<string>& geoPhyNamesList=cGeomData::getInstance()->phyNamesList;
    int numOfGeoPts;
    int numOfGeoTris;
    int numOfGeoPhyNames;

    //geom features
    vector<cTri*> geoFFacesList; //geom feature faces list
    vector<cFeaturePt*> geoFPtsList; //geom feature pts list
    vector<cFeatureEdge*> geoFEdgesList; //geom feature edges list
    vector<double> ptInGeom;

    //msh data
    vector<vector<double>> mshPtsList;
    int mshPtIndxCount=0;
    int mshVolIndxCount = 0;
    //vector<cFace> mshExtFacesList;
    vector<cFace*> mshIntlFacesList;
    vector<cFace*> mshBFacesList;
    vector<cFace*> mshAllFacesList;
    vector<cBoundary> mshBdsList;
    int num_mshBFaces, num_mshIntlFaces;
    //vector<int> eleVolList;
    

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
    void outputNodeName(cOctNode *node);
    void outputMshPts(cOctNode *node);
    void outputMshFaces(cOctNode *node);


    void defineBody(vector<double> _ptInGeom);
    void buildOctree();

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


    void setup_geoFFacesList();
    void setup_geoFPtsList();
    void setup_geoFEdgesList();
    int numFFaces();

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

    void setup_mshPtList();
    void addMshPtsOfNode(cOctNode* node);
    bool isPtSame(vector<double>& pt1, vector<double>& pt2);
    void setup_nodeMshFaces();

    void setup_boundaryMshFace();
    void setup_intlMshFace();
    void findMshIntlFaces();
    void put2List(int flag,cFace* currentMshFace,cFace* listMshFace);

    void setup_nodePhyName();
    void setup_mshBdsList();


    void saveAsGmsh(const char* _fileName);
    void saveAsOFMesh();
    void saveAsOFMeshPts();
    void saveAsOFMeshFaces();
    void saveAsOFMeshNeis();
    void saveAsOFMeshOwns();
    void saveAsOFMeshBds();
};

//vector<cTri> geoFFacesList; //geom feature faces list
//vector<cFeaturePt> geoFPtsList; //geom feature pts list
//vector<cFeatureEdge> geoFEdgesList; //geom feature edges list

#endif
