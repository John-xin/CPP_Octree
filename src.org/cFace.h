/*
 * cFace.h
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#ifndef CFACE_H_
#define CFACE_H_

#include <iostream>
#include <vector>
//#include <set>
//#include <stdlib.h>
//#include <stdio.h>
//#include <math.h>
//#include <string>
//#include <sstream>
//#include <algorithm> // find, sort
//#include <utility>   // pair
#include "cTri.h"

using namespace std;

class cFace{
public:
	cFace();
	~cFace();
	string nid;
	int nbr; //neighbour
	int own; //owner
	int label;
	int nPts;
	vector<vector<double> > ptsList;
	vector<int> ptIndxList;
	vector<double> N;
	vector<double> centroid;


	bool bFlag; //boundary flag;
	string phyName;
	int phyNameIndx;
    vector<double> low;
    vector<double> upp;
    int mshVolIndx;
    int state; //0-boundary face ; 1 - internal face
    cOctNode* node;

    void getN();
    double getAngle(cFace *surface);
    vector<double> getCentroid();
    void findPhyName(vector<cTri*> _geoFFacesList);
    int findGeoFaceWithMinDist(cLine _ray,vector<cTri*> _geoFFacesList);
    void changeOrder();
    int findFaceRelationship(cFace* f1, cFace* f2);
    bool isPtInFace(vector<double> pt);

};



#endif /* CFACE_H_ */
