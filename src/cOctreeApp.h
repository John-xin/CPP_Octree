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
class cOctreeOFMesh;

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
    cOctreeOFMesh* ofMesh;
    cOctreeApp();
    ~cOctreeApp();

    void defineBody(vector<double> _ptInGeom);
    void buildOctree();
    void buildOFMesh();
    
    void saveAsGmsh(const char* _fileName);
    void saveAsOFMesh();
};

#endif