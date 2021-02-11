#ifndef OCTREEUTIL_H
#define OCTREEUTIL_H

#include <iostream>
#include <vector>
#include <queue>
#include <deque>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm> // find, sort
#include <utility>   // pair
#include <fstream>


using namespace std;

class cOctree;
class cOctNode;
class cFeatureFace;
class cFeaturePt;

class cOctreeUtil {
public:
    cOctreeUtil();
    ~cOctreeUtil();

    cOctree* octree;

    void setup_octree(cOctree* _octree);

    //test program
    void outputNodeName(cOctNode* node);
    void outputMshPts(cOctNode* node);
    void outputMshFaces(cOctNode* node);

    void output_geoFFacesList(vector<cFeatureFace*>& geoFFacesList);
    void output_geoFPtsList(vector<cFeaturePt*>& geoFPtsList, vector<vector<double> >& geoPts3DList);
    void output_octree(const char* _fileName, cOctree* tree);
    void output_nodes(const char* _fileName, vector<cOctNode*> nodesList);

    void checkNodes_ErrorType1(vector<cOctNode*> nodesList);
    bool isErrorType1(cOctNode* node);
    vector<cOctNode*>  fix_ErrorType1(cOctNode* node);

 
};


#endif