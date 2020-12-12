#ifndef CFEATURE_H
#define CFEATURE_H

#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm> // find, sort
#include <utility>   // pair

//#include "cOctNode.h"

// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cOctNode;

class cFeaturePt {
public:

    int indx;
    cFeaturePt();
    cFeaturePt(int _indx);
    ~cFeaturePt();
    bool isInNode(cOctNode *node);
    bool isNewInList(vector<cFeaturePt*> _geoFPtsList);
};

class cFeatureEdge {
public:

    int indx;
    vector<vector<double> > pts3D_db2;
    vector<int> ptIndx_int2;
    vector<double> low, upp;
    cFeatureEdge();
    cFeatureEdge(int _indx, vector<vector<double> > _pts3D_db2);
    cFeatureEdge(int _indx, vector<int> _ptIndx_int2);
    ~cFeatureEdge();
    bool isInNode(cOctNode *node);
    void getLowerVert();
    void getUpperVert();
    void getV();
};



#endif
