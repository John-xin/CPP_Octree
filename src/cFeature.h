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
#include "commonFunc.h"
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

//******************************************************************
//******************************************************************
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

//******************************************************************
//******************************************************************

class cLine { //cFeatureFace would use cLine to define ray
public:
    vector<double> p0,p1,dir;
    cLine();
    cLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~cLine();
    void getDir();
    void getP1();
};

// ------------------------------------------------------

class cFeatureFace { //acts as feature faces of geom data
public:
    double D;
    int indx; //geom face indx
    vector<vector<double> > pts3D_db3; //geom pt3D
    vector<int> ptIndx_int3; //geom pt Indx
    int phyNameIndx;
    string phyName;
    vector<double> N;
    vector<double> low, upp;

    cFeatureFace();
    cFeatureFace(int _indx, vector<vector<double> > _pts3D_db3);
    cFeatureFace(int _indx, int _phyNameIndx, vector<int> _ptIndx_int3);
    ~cFeatureFace();

    bool isInNode(cOctNode *node);
    bool isPointInTri(vector<double> &p);
    bool rayPlaneIntersectPoint(cLine &ray, bool entryOnly);
    bool rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s);
    void getN();
    void getD();
    void getV();
    void getLowerVert();
    void getUpperVert();
    double getAngle(cFeatureFace* surface);
    bool isPointSame(int& vertLabel1, int& vertLabel2);
    bool isAdjacent(cFeatureFace* surface);
    vector<int> getAdjacentEdge(cFeatureFace* surface);
};


#endif
