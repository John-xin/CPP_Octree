/*
 * cFeature.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#include "cFeature.h"
#include "cOctNode.h"
#include "cGeomData.h"

cFeaturePt::cFeaturePt()
{
    // Default cFeaturePt constructor
    indx = 0;
}

cFeaturePt::cFeaturePt(int _indx)
{
    // cFeaturePt constructor with label and vert
    indx    = _indx;
}

// cFeaturePt destructor
cFeaturePt::~cFeaturePt() {
    //cout << "Destroying cFeaturePt" << endl;
}

bool cFeaturePt::isInNode(cOctNode *node)
{
	vector<double> vert;
	vert=cGeomData::getInstance()->pts3DList[indx];
    if (vert[0] > node->upp[0]) return false;
    if (vert[1] > node->upp[1]) return false;
    if (vert[2] > node->upp[2]) return false;
    if (vert[0] < node->low[0]) return false;
    if (vert[1] < node->low[1]) return false;
    if (vert[2] < node->low[2]) return false;

    //if vert in or on the node then return true
    return true;
}

bool cFeaturePt::isNewInList(vector<cFeaturePt*> _geoFPtsList){
	for(unsigned i=0;i<_geoFPtsList.size();i++){
		if(indx==_geoFPtsList[i]->indx){
			return false;
		}
	}
	return true;
}

// ------------------------------------------------------

cFeatureEdge::cFeatureEdge()
{
    // Default cFeautureEdge constructor
    indx = 0;
    pts3D_db2.resize(2);
    for (vector<vector<double> >::iterator it=pts3D_db2.begin(); it!=pts3D_db2.end(); ++it)
        (*it).resize(2,0.0);

    getV();
    getLowerVert();
    getUpperVert();
}

cFeatureEdge::cFeatureEdge(int _indx, vector<vector<double> > _pts3D_db2)
{
    // cFeautureEdge constructor with label and vertices
    indx    = _indx;
    pts3D_db2 = _pts3D_db2;
    getLowerVert();
    getUpperVert();
}

cFeatureEdge::cFeatureEdge(int _indx, vector<int> _ptIndx_int2)
{
    // cFeautureEdge constructor with label and vertLabels
    indx    = _indx;
    ptIndx_int2 = _ptIndx_int2;
    getV();
    getLowerVert();
    getUpperVert();

}

// cFeautureEdge destructor
cFeatureEdge::~cFeatureEdge() {
    //cout << "Destroying cFeautureEdge" << endl;
}

void cFeatureEdge::getV()
{
    //get (vertices 2 x 3) by (vertlabels 2 x 1)
    int indx;
    for(int i=0; i<2; i++) {
         indx=ptIndx_int2[i];
         pts3D_db2.push_back(cGeomData::getInstance()->pts3DList[indx]);
    }
}


void cFeatureEdge::getLowerVert()
{
    // Lower vertices of cTri bounding box
    low.resize(3,1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<2; i++) {
            if (pts3D_db2[i][j] < low[j])
            {
                low[j] = pts3D_db2[i][j];
            }
        }
    }
}

void cFeatureEdge::getUpperVert()
{
    // Upper vertices of cTri bounding box
    upp.resize(3,-1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<2; i++) {
            if (pts3D_db2[i][j] > upp[j])
            {
                upp[j] = pts3D_db2[i][j];
            }
        }
    }
}

bool cFeatureEdge::isInNode(cOctNode *node)
{
    // Tests if bounding box of cTri is inside of or overlapping the given cOctNode
    // This is a simple test and even if bounding box is found to be inside the
    // cOctNode, the cTri itself may not be
    if (low[0] > node->upp[0]) return false;
    if (low[1] > node->upp[1]) return false;
    if (low[2] > node->upp[2]) return false;
    if (upp[0] < node->low[0]) return false;
    if (upp[1] < node->low[1]) return false;
    if (upp[2] < node->low[2]) return false;
    return true;
}
