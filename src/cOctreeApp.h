#ifndef COCTREEAPP_H
#define COCTREEAPP_H
// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include "cFeature.h"
#include "cGeomData.h"
#include "cOctNode.h"


// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cOctree;
class cOctreeUtil;

class cOctreeApp {
public:
    static const int MAX_OCTREE_LEVELS = 3;
    static const int MIN_OCTREE_LEVELS = 1;
    static const int MAX_OCTNODE_FEATS = 2;
    //octree setup depends on geomData - run geomData read STL first
    //simple geom data
    cGeomData* geom;
    cOctree* octree;
    cOctreeUtil* util;
    cOctreeApp();
    ~cOctreeApp();

    //msh data
    vector<vector<double>> mshPtsList;
    int mshPtIndxCount = 0;
    int mshVolIndxCount = 0;
    //vector<cFace> mshExtFacesList;
    vector<cFace*> mshIntlFacesList;
    vector<cFace*> mshBFacesList;
    vector<cFace*> mshAllFacesList;
    vector<cBoundary> mshBdsList;
    int num_mshBFaces, num_mshIntlFaces;
    //vector<int> eleVolList;

    void defineBody(vector<double> _ptInGeom);
    void buildOctree();

    void setup_mshPtList();
    void addMshPtsOfNode(cOctNode* node);
    bool isPtSame(vector<double>& pt1, vector<double>& pt2);
    void setup_nodeMshFaces();

    void setup_boundaryMshFace();
    void setup_intlMshFace();
    void findMshIntlFaces();
    void put2List(int flag, cFace* currentMshFace, cFace* listMshFace);

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

#endif