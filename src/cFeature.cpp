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

//******************************************************************
//******************************************************************

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

//******************************************************************
//******************************************************************
cLine::cLine()
{
    // Default constructor for cLine
    // Default line is unit vector along x-axis
    p0.resize(3,0.0); p1.resize(3,0.0); dir.resize(3,0.0);
    p1[0]=1.0; dir[0]=1.0;
}

cLine::cLine(vector<double> &_p0, vector<double> &p1_dir, int isP1orDir)
{
    // cLine constructor with p0 and p1 or dir
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = _p0;
    if (isP1orDir==0) {
        p1 = p1_dir;
        getDir(); }
    else if (isP1orDir==1) {
        dir = p1_dir;
        getP1(); }
}

// cLine destructor
cLine::~cLine() {}

void cLine::getDir()
{
    // Get unit vector defining direction of cLine
    vector<double> p0p1(3); double dmag=0.0;
    for (unsigned int i=0; i<3; i++) {
        p0p1[i] = p1[i]-p0[i];
        dmag   += pow(p0p1[i],2.0);
    }
    dmag = sqrt(dmag);
    dir  = p0p1;
    for (vector<double>::iterator it=dir.begin(); it!=dir.end(); ++it)
        *it /= dmag;
}

void cLine::getP1()
{
    // Get a point on the cLine, p1, located 1.0 units away from the origin, p0
    vector<double> p1(3);
    for (unsigned int i=0; i<3; i++)
        p1[i] = p0[i]+dir[i];
}
// ------------------------------------------------------
cFeatureFace::cFeatureFace()
{
    // Default cTri constructor
    indx = 0;
    pts3D_db3.resize(3);
    phyName = "unAssigned";
    phyNameIndx = -100;
    for (vector<vector<double> >::iterator it=pts3D_db3.begin(); it!=pts3D_db3.end(); ++it)
        (*it).resize(3,0.0);
    pts3D_db3[1][0]=1.0; pts3D_db3[2][1]=1.0;
    getN();
    getD();
    getLowerVert();
    getUpperVert();
}

cFeatureFace::cFeatureFace(int _indx, vector<vector<double> > _pts3D_db3)
{
    // cTri constructor with label and pts3D
    indx    = _indx;
    pts3D_db3 = _pts3D_db3;
    phyName = "unAssigned";
    phyNameIndx = -100;
    getN();
    getD();
    getLowerVert();
    getUpperVert();
}

cFeatureFace::cFeatureFace(int _indx, int _phyNameIndx, vector<int> _ptIndx_int3)
{
    // cTri constructor with label and vertLabels
    indx    = _indx;
    phyNameIndx=_phyNameIndx;
    phyName=cGeomData::getInstance()->phyNamesList[phyNameIndx];
    ptIndx_int3 = _ptIndx_int3;
    getV();
    getN();
    getD();
    getLowerVert();
    getUpperVert();
}

// cTri destructor
cFeatureFace::~cFeatureFace() {
    //cout << "Destroying cTri" << endl;
}

void cFeatureFace::getV()
{
    //get (vertices 3 x 3) by (vertlabels 3 x 1)
    int indx;
    for(int i=0; i<3; i++) {
        indx=ptIndx_int3[i];
        pts3D_db3.push_back(cGeomData::getInstance()->pts3DList[indx]);
    }
}

void cFeatureFace::getN()
{
    // Get cTri face normal

    vector<vector<double> > v=pts3D_db3; //3 x 3
    vector<double> v1(3),v2(3),v3;

    // vector 1 = vert 0 - vert 1 (dir 1->0)  ---- vector 2 = vert 0 - vert 2
    for (unsigned int i=0; i<3; i++) {
        v1[i] = v[0][i] - v[1][i];
        v2[i] = v[0][i] - v[2][i]; }
    v3 = commonFunc::getInstance()->crossProduct(v1,v2);
    double v3mag = sqrt(commonFunc::getInstance()->dotProduct(v3,v3));
    N = v3;
    for (vector<double>::iterator it=N.begin(); it!=N.end(); ++it)
        *it /= v3mag;
}

void cFeatureFace::getD()
{
    // Perp distance of cTri face from origin which, along with the face normal,
    // defines the plane of the cTri face

    D = commonFunc::getInstance()->dotProduct(pts3D_db3[0],N);
}

void cFeatureFace::getLowerVert()
{
    // Lower vertices of cTri bounding box
    low.resize(3,1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            if (pts3D_db3[i][j] < low[j])
            {
                low[j] = pts3D_db3[i][j];
            }
        }
    }
}

void cFeatureFace::getUpperVert()
{
    // Upper vertices of cTri bounding box
    upp.resize(3,-1.0e+30);
    for (int j=0; j<3; j++) {
        for (int i=0; i<3; i++) {
            if (pts3D_db3[i][j] > upp[j])
            {
                upp[j] = pts3D_db3[i][j];
            }
        }
    }
}

bool cFeatureFace::isInNode(cOctNode *node)
{
    // Tests if bounding box of cTri is inside of or overlapping the given cOctNode
    // This is a simple test and even if bounding box is found to be inside the
    // cOctNode, the cTri itself may not be
	double tol=1e-6;
    if (low[0] - node->upp[0]>tol) return false;
    if (low[1] - node->upp[1]>tol) return false;
    if (low[2] - node->upp[2]>tol) return false;
    if (upp[0] - node->low[0]<-tol) return false;
    if (upp[1] - node->low[1]<-tol) return false;
    if (upp[2] - node->low[2]<-tol) return false;
    return true;
}

bool cFeatureFace::isPointInTri(vector<double> &p)
{
    // Determines if point p is within the cTri by computing and
    // testing the barycentric coordinates (u, v, w) of p

    // Find Barycentric coordinates of point (u,v,w)
    vector<double> v0 = commonFunc::getInstance()->vectSubtract(pts3D_db3[1],pts3D_db3[0]);
    vector<double> v1 = commonFunc::getInstance()->vectSubtract(pts3D_db3[2],pts3D_db3[0]);
    vector<double> v2 = commonFunc::getInstance()->vectSubtract(p,pts3D_db3[0]);
    double d00 = commonFunc::getInstance()->dotProduct(v0, v0);
    double d01 = commonFunc::getInstance()->dotProduct(v0, v1);
    double d11 = commonFunc::getInstance()->dotProduct(v1, v1);
    double d20 = commonFunc::getInstance()->dotProduct(v2, v0);
    double d21 = commonFunc::getInstance()->dotProduct(v2, v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0 - v - w;

    // Use Barycentric coordinates to work out if point lies within the cTri element
    double tol = 0.0;
    return ((v>=tol) && (w>=tol) && (u>=tol));
}

bool cFeatureFace::rayPlaneIntersectPoint(cLine &ray, bool entryOnly=false)
{
    // Tests if ray intersects with the cTri face
    // NOTE: Provide option to just check for entry intersections, not both
    //       entries/exits. This should cut down checking somewhat.
    double tol  = 1.0e-06;
    double sDen = commonFunc::getInstance()->dotProduct(ray.dir,N);
    if ((entryOnly && sDen>tol) || (!entryOnly && fabs(sDen)>tol))
    {
        double sNum = D - commonFunc::getInstance()->dotProduct(ray.p0,N);
        double s = sNum / sDen;
        vector<double> p = commonFunc::getInstance()->vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    }
    return false;
}

bool cFeatureFace::rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s)
{
    // Tests if ray intersects with the cTri face
    // Returns the coordinates of the intersection point and the distance, s, from
    // the origin of the ray
    double tol  = 1.0e-06;
    double sDen = commonFunc::getInstance()->dotProduct(ray.dir,N);
    if (fabs(sDen)> tol) // Normals cannot be perpendicular such that dot product equals 0
    {
        double sNum = D - commonFunc::getInstance()->dotProduct(ray.p0,N);
        s = sNum / sDen;
        p = commonFunc::getInstance()->vectAdd(ray.p0,ray.dir,s);
        return isPointInTri(p);
    }
    return false;
}

double cFeatureFace::getAngle(cFeatureFace *surface)
{
    double angle;
    angle=acos(commonFunc::getInstance()->dotProduct(N,surface->N))*180/3.141592653;
    if (angle==0){ angle=180;}
    return angle;

}

bool cFeatureFace::isPointSame(int& vertLabel1, int& vertLabel2)
{
    if(vertLabel1!=vertLabel2) return false;
    return true;
}

bool cFeatureFace::isAdjacent(cFeatureFace* surface)
{
	int flag=0;
	for(int i=0; i<3; i++){
	    for(int j=0; j<3; j++){
	        if (isPointSame(ptIndx_int3[i], surface->ptIndx_int3[j])){
		    flag=flag+1;
		}
	    }
	}
	if (flag==2){
	    return true;
	}else{
	    return false;
	}

}

vector<int> cFeatureFace::getAdjacentEdge(cFeatureFace* surface)
{
	//make sure isAjacent is true then call this function
   vector<int> _ptIndx_int2;
   for(int i=0; i<3; i++){
	    for(int j=0; j<3; j++){
	        if (isPointSame(ptIndx_int3[i], surface->ptIndx_int3[j])){
		              _ptIndx_int2.push_back(ptIndx_int3[i]);
		          }
		  }
	 }
	 return _ptIndx_int2;
}
