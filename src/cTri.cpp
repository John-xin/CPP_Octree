/*
 * cTri.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#include "cTri.h"
#include "cOctNode.h"
#include "cGeomData.h"
#include "commonFunc.h"

cTri::cTri()
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

cTri::cTri(int _indx, vector<vector<double> > _pts3D_db3)
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

cTri::cTri(int _indx, int _phyNameIndx, vector<int> _ptIndx_int3)
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
cTri::~cTri() {
    //cout << "Destroying cTri" << endl;
}

void cTri::getV()
{
    //get (vertices 3 x 3) by (vertlabels 3 x 1)
    int indx;
    for(int i=0; i<3; i++) {
        indx=ptIndx_int3[i];
        pts3D_db3.push_back(cGeomData::getInstance()->pts3DList[indx]);
    }
}

void cTri::getN()
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

void cTri::getD()
{
    // Perp distance of cTri face from origin which, along with the face normal,
    // defines the plane of the cTri face

    D = commonFunc::getInstance()->dotProduct(pts3D_db3[0],N);
}

void cTri::getLowerVert()
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

void cTri::getUpperVert()
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

bool cTri::isInNode(cOctNode *node)
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

bool cTri::isPointInTri(vector<double> &p)
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

bool cTri::rayPlaneIntersectPoint(cLine &ray, bool entryOnly=false)
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

bool cTri::rayPlaneIntersectPoint(cLine &ray, vector<double> &p, double &s)
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

double cTri::getAngle(cTri *surface)
{
    double angle;
    angle=acos(commonFunc::getInstance()->dotProduct(N,surface->N))*180/3.141592653;
    if (angle==0){ angle=180;}
    return angle;

}

bool cTri::isPointSame(int& vertLabel1, int& vertLabel2)
{
    if(vertLabel1!=vertLabel2) return false;
    return true;
}

bool cTri::isAdjacent(cTri* surface)
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

vector<int> cTri::getAdjacentEdge(cTri* surface)
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



