#ifndef CTRI_H
#define CTRI_H

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

#include "cLine.h"
#include "commonFunc.h"

// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cOctNode;

class cTri { //cTri acts as feature faces of geom data
public:
    double D;
    int indx; //geom tri indx
    vector<vector<double> > pts3D_db3; //geom pt3D
    vector<int> ptIndx_int3; //geom pt Indx
    int phyNameIndx;
    string phyName;
    vector<double> N;
    vector<double> low, upp;

    cTri();
    cTri(int _indx, vector<vector<double> > _pts3D_db3);
    cTri(int _indx, int _phyNameIndx, vector<int> _ptIndx_int3);
    ~cTri();

    bool isInNode(cOctNode *node);
    bool isPointInTri(vector<double> &p);
    bool rayPlaneIntersectPoint(cLine &ray, bool entryOnly);
    bool rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s);
    void getN();
    void getD();
    void getV();
    void getLowerVert();
    void getUpperVert();
    double getAngle(cTri* surface);
    bool isPointSame(int& vertLabel1, int& vertLabel2);
    bool isAdjacent(cTri* surface);
    vector<int> getAdjacentEdge(cTri* surface);
};



#endif
