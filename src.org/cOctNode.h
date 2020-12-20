#ifndef COCTNODE_H
#define COCTNODE_H

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


#include "cFace.h"
#include "cLine.h"



// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace std;

class cFeaturePt;
class cFeatureEdge;
class cTri;

class cOctNode {
public:

    //static const int MAX_OCTNODE_OBJECTS  = 1;
    static const int NUM_BRANCHES_OCTNODE = 8;
      
    double size; //node cube length
    int level;
    string nid;
    vector<double> position; //assigned from constructor
    vector<cOctNode> branches;
    cOctNode* parent;
    vector<vector<cOctNode*>> nbr;

    vector<double> low, upp;
    int mshVolIndx;
    int state;//0-boundary node ; -1 - exterior node; 1 - internal node

    vector<cTri*> geoFFacesList; //geom feature faces list
    vector<cFeaturePt*> geoFPtsList; //geom feature pts list
    vector<cFeatureEdge*> geoFEdgesList; //geom feature edges list

    vector<int> geoFFacesIndxList; //geam tri
//    vector<int> geoFPtsIndxList; //feature pts
//    vector<int> geoFEdgesIndxList; //feature geoFEdgesIndxList

    vector<vector<double>> mshPts3DList;
    vector<int> mshPtsIndxList;
    vector<bool> mshPtsRepeatedList;
    vector<cFace> mshFacesList;//element faces


    cOctNode();
    cOctNode(int _level, string _nid, vector<double> _position, double _size,
    		vector<cFeaturePt*> _geoFPtsList,
    		vector<cFeatureEdge*> _geoFEdgesList,
    		vector<cTri*> _geoFFacesList,
    		cOctNode* _parent);
    ~cOctNode();
    bool isLeafNode();
    int numOfGeoFFaces();
    int numOfGeoFPts();
    int numOfGeoFEdges();
    //void addGeoFFacesIndx(int _indx);//geom tri
    //void addGeoFPtsIndx(int _indx);//add feature geoFPtsIndxList
    void setGeoFFacesIndxList();
    void addNode(int _level, string _nid, vector<double> _position, double _size,
    			vector<cFeaturePt*> _geoFPtsList,
    			vector<cFeatureEdge*> _geoFEdgesList,
    			vector<cTri*> _geoFFacesList,
    			cOctNode* _parent);
    void getLowUppVerts();//calc low and upp by position
    void calMshPts3D();
    void calMshFaces();
    bool boxRayIntersect(cLine &ray);
    bool sphereRayIntersect(cLine &ray);
    void update_MshFaces_PtIndxList();
    void updateMshFaces_mshVolIndx(int _mshVolIndx);
    void removeExtraFeats();
    bool isPtInNode(vector<double> pt);
};




// 
#endif
