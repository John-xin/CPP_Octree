// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include "cOctNode.h"
#include "cFeature.h"
#include "cOctree.h"
#include "commonFunc.h"
#include "cOctreeOFMesh.h"

//**********************************************************************
//**********************************************************************
cBoundary::cBoundary()
{
    // Default cTri constructor
    phyName="unAssigned";
    phyNameIndx=-100;
    patchType="unAssigned";
    startFace=-100;
    nFaces=0;
}

cBoundary::~cBoundary() {
    //cout << "Destroying cBoundary" << endl;
}

//*************************************************************************
//***********************************************************************
cFace::cFace()
{
	nid="unAssigned";
    nbr=-100;
    own=-100;
    label=-100;
    nPts=0;
    bFlag=false;
    phyName = "unAssigned";
    phyNameIndx=-100;
    state = FaceState::unassigned;
    exportState = false;
    node=NULL;

}

cFace::~cFace(){}

void cFace::updatePtsList()
{
    ptsList.resize(0);
    for (int indx : ptIndxList) {
        vector<double> pt = cOctreeOFMesh::getInstance()->mshPtsList[indx];
        ptsList.push_back(pt);
    }
}

void cFace::getN()
{
    // Get cTri face normal

    vector<vector<double> > v=ptsList; //3 x 3
    vector<double> v1(3),v2(3),v3;
    vector<double> center;

    // vector 1 = vert 0 - vert 1 (dir 1->0)  ---- vector 2 = vert 0 - vert 2
    center = getCentroid(v);
    v1 = commonFunc::getInstance()->vectSubtract(v[0],center);
    v2 = commonFunc::getInstance()->vectSubtract(v[1], center);
    v3 = commonFunc::getInstance()->crossProduct(v1,v2);
    double v3mag = sqrt(commonFunc::getInstance()->dotProduct(v3,v3));
    N = v3; 
    for (vector<double>::iterator it = N.begin(); it != N.end(); ++it){ *it /= v3mag; }
        
    
}

double cFace::getAngle(cFace* surface)
{
	double tol=0.000001;
    getN();
    surface->getN();
	double angle;
    angle=acos(commonFunc::getInstance()->dotProduct(N,surface->N))*180/3.141592653;
    if (fabs(angle)< tol){ angle=180.0;}
    return angle;

}

vector<double> cFace::getCentroid(vector<vector<double>> _ptsList) {
	vector<double> pt;
	pt.resize(3);
	for(unsigned int i=0; i< _ptsList.size();i++){
		pt[0]+= _ptsList[i][0];
		pt[1]+= _ptsList[i][1];
		pt[2]+= _ptsList[i][2];
	}
	pt[0]/= _ptsList.size();
	pt[1]/= _ptsList.size();
	pt[2]/= _ptsList.size();
	return pt;
}

void cFace::findPhyName(vector<cFeatureFace*> _geoFFacesList) {
	int indx1=-1;
	int indx2=-1;
	getN();
	vector<double> centroid, pt2;
	centroid=getCentroid(ptsList);
	cLine ray1(centroid,N,1);
	indx1=findGeoFaceWithMinDist(ray1,_geoFFacesList);

//	pt2.resize(3);
//	pt2[0]=centroid[0]-0.01;
//	pt2[1]=centroid[1]-0.01;
//	pt2[2]=centroid[2]-0.01;
//	cLine ray2(low,N,1);
//	indx2=findGeoFaceWithMinDist(ray2,_geoFFacesList);

	ostringstream oss;

	if( indx1==-1){//no intersections
		phyNameIndx=-100;
		phyName="unAssigned";
	}else{
		phyNameIndx=_geoFFacesList[indx1]->phyNameIndx;
		phyName=_geoFFacesList[indx1]->phyName;

//		if(_geoFFacesList[indx1].phyName !=_geoFFacesList[indx2].phyName){
//			cout<<"Warning:" << nid << "-> phyName needs checking" <<"\n";
//			oss<<"Warning:" << nid << "-> phyName maybe geoFace-";
//			oss<<oss.str() << _geoFFacesList[indx1].indx << "-" <<_geoFFacesList[indx1].phyName.c_str() <<" or geoFace-" ;
//			oss<<oss.str() << _geoFFacesList[indx2].indx << "-" <<_geoFFacesList[indx2].phyName.c_str() <<"\n\n";
//			cout << oss.str();
//		}
	}

	//	if(nid=="node-0-3-6->face-0"){
	//		cout<<"selected index - " << indx <<"\n";
	//		for(unsigned int i=0; i<_geoFFacesList.size();i++){
	//			cout<< i <<" geoFace - " <<_geoFFacesList[i].indx << "-" << _geoFFacesList[i].phyName.c_str() << "\n" ;
	//		}
	//
	//	}

}

int cFace::findGeoFaceWithMinDist(cLine _ray, vector<cFeatureFace*> _geoFFacesList) {
	cFeatureFace* geoFFace;
	vector<double> pt;
	double s=0;
	double s_min=1e10;
	int indx=-1;
	bool flag;

	for(unsigned int i=0; i<_geoFFacesList.size();i++){
		geoFFace=_geoFFacesList[i];
		flag=geoFFace->rayPlaneIntersectPoint(_ray, pt, s);

		if(flag==true){
//			if(nid=="node-0-3-6->face-0"){
//				cout<< i <<" geoFace - " <<geoFFace.indx << " - " << geoFFace.phyName.c_str() << " - " << s <<"\n" ;
//			}
			if(fabs(s)<s_min){
				s_min=fabs(s);
				indx=i;
			}
		}
	}
	//if return -1 -> no intersections


//	if(nid=="node-0-3-6->face-0"){
//	cout<<"selected index - " << indx <<"\n";}

	return indx;
}

int numOfOrderChanged=0;

void cFace::changeOrder() {
	vector<int> _ptLabelList=ptIndxList;
	int numOfPts=ptIndxList.size();
	for(int i=0; i<numOfPts;i++){
		ptIndxList[i]=_ptLabelList[numOfPts-1-i];
	}
	numOfOrderChanged++;
	//cout<<"order changed -> "<< nid <<" ("<< ptIndxList[0] <<"," << ptIndxList[1]  <<"," << ptIndxList[2]  <<"," << ptIndxList[3]<<") \n";
}

int cFace::findFaceRelationship(cFace* f1, cFace* f2) {
	//assume f1 and f2 are rectangles
	double tol=1e-6;
    double angle = fabs(f1->getAngle(f2) - 180);
	if(angle <tol){
		if(fabs(f1->low[0]-f2->low[0])<tol && fabs(f1->low[1]-f2->low[1])<tol && fabs(f1->low[2]-f2->low[2])<tol &&
		   fabs(f1->upp[0]-f2->upp[0])<tol && fabs(f1->upp[1]-f2->upp[1])<tol && fabs(f1->upp[2]-f2->upp[2])<tol)
		{
			return 1;//f1 = f2

		}else if(f1->isPtInFace(f2->low) && f1->isPtInFace(f2->upp))	{
			return 3;//f1>= f2

		}else if(f2->isPtInFace(f1->low) && f2->isPtInFace(f1->upp))	{
		    return 2;//f2>= f1
		}else{
			return 0;
		}
	}
	return 0;

}


bool cFace::isPtInFace(vector<double> pt) {
	//assume face is a rectangle
	double tol=1e-6;
	if(pt.size()==3){
	    if (low[0] - pt[0]> tol) return false;
	    if (low[1] - pt[1]> tol) return false;
	    if (low[2] - pt[2]> tol) return false;
	    if (upp[0] - pt[0]< -tol) return false;
	    if (upp[1] - pt[1]< -tol) return false;
	    if (upp[2] - pt[2]< -tol) return false;
	    return true;
	}else{
		cout<<"error: pt size in <isPtInFace Func> is not 3" << "\n";
		return false;
	}
}

vector<int> cFace::getInFacePts(vector<int> _ptsIndxList)
{
    vector<int> inFacePts;
    vector<double> pt;
    for (size_t i = 0; i < _ptsIndxList.size();i++) {
        pt = cOctreeOFMesh::getInstance()->mshPtsList[_ptsIndxList[i]];
        if (isPtInFace(pt)) {
            inFacePts.push_back(_ptsIndxList[i]);
        }
    }
    return inFacePts;
}

vector<int> cFace::getInOrderPts(vector<int> _ptsIndxList)
{
    vector<int> inOrderPts;

    //get 3d points
    vector<vector<double>> givenPts;
    for (int indx: _ptsIndxList){
        vector<double> pt = cOctreeOFMesh::getInstance()->mshPtsList[indx];
        givenPts.push_back(pt);
    }

    vector<double> center = getCentroid(givenPts);
    
    //get angle -> center to point
    // angle = atan2(det, dot) -> to output angle range (-pi to +pi)
    // dot = v1 * v2
    // det = N * (v1 x v2) for 3d points in a plane
    vector<double> angles;
    vector<double> startRay= commonFunc::getInstance()->vectSubtract(givenPts[0], center);
    getN();
    for (vector<double> pt : givenPts) {
        vector<double> nextRay=commonFunc::getInstance()->vectSubtract(pt, center);
        double dot = commonFunc::getInstance()->dotProduct(startRay, nextRay);
        vector<double> vect = commonFunc::getInstance()->crossProduct(startRay, nextRay);
        double det = commonFunc::getInstance()->dotProduct(N, vect);
        double angle = atan2(det, dot);
        angles.push_back(angle);
    }

    //sort pairs by angle
    vector< pair<double, int> > pairs;
    for (size_t i = 0; i < _ptsIndxList.size();i++) {
        pairs.push_back(make_pair(angles[i], _ptsIndxList[i]));
    }
    
    sort(pairs.begin(),pairs.end());


    //get back ordered indx list
    for (pair<double, int> p : pairs) {
        inOrderPts.push_back(p.second);
    }

    return inOrderPts;
}

//***********************************************************************
//**********************************************************************

cOctNode::cOctNode()
{
    // Default octNode constructor
    level = 0;
    depth = 1;
    nid   = "";
    size  = 1.0;
    state = NodeState::unassigned;
    parent=NULL;
    position.resize(3,0.0);
    getLowUppVerts();

    //geoFFacesList.reserve(MAX_OCTNODE_OBJECTS);
}

cOctNode::cOctNode(int _level, string _nid, vector<double> _position, double _size,
		vector<cFeaturePt*> _geoFPtsList,
		vector<cFeatureEdge*> _geoFEdgesList,
		vector<cFeatureFace*> _geoFFacesList,
		cOctNode* _parent)
{
    // octNode constructor with level, node id (nid), position and size
    level    = _level;
    depth = 1;
    nid      = _nid;
    position = _position;
    size     = _size;
    mshVolIndx = -100;
    state = NodeState::unassigned;
    parent=_parent;
    geoFPtsList=_geoFPtsList;
    geoFEdgesList=_geoFEdgesList;
    geoFFacesList=_geoFFacesList;
    getLowUppVerts();
    mshPtsRepeatedList.resize(8);
    mshPtsIndxList.resize(8);
    mshFacesList.resize(6);

    calMshPts3D();//based on getLowUppVerts();
    //calMshFaces();//must after calMshPts3D();
    removeExtraFeats();
}

cOctNode::~cOctNode() {
    //cout << "Calling destructor for cOctnode " << nid << endl;
}

bool cOctNode::isLeafNode()
{
    // Checks if cOctNode is a leaf node by counting the number of branches. A
    // leaf node has no branches
    return children.size()==0;
}

void cOctNode::getLowUppVerts()
{
    // Get coordinates of the lower and upper vertices of the cOctNode
    low.resize(3);
    upp.resize(3);
    double halfSize = size/2.0;
    for (int i=0; i<3; i++) {
        low[i] = position[i] - halfSize;
        upp[i] = position[i] + halfSize;
    }
}

void cOctNode::calMshPts3D()//init mshPts3DList
{
	vector<double> pt;
    pt.resize(3);

	mshPts3DList.push_back(low);//left_down z=0     0
	pt[0]=upp[0]; pt[1]=low[1];	pt[2]=low[2];
	mshPts3DList.push_back(pt);//right_down z=0     1
	pt[0]=upp[0]; pt[1]=upp[1];	pt[2]=low[2];
	mshPts3DList.push_back(pt);//right_up z=0       2
	pt[0]=low[0]; pt[1]=upp[1];	pt[2]=low[2];
	mshPts3DList.push_back(pt);//left_up z=0        3

	pt[0]=low[0]; pt[1]=low[1];	pt[2]=upp[2];
	mshPts3DList.push_back(pt);//left_down z=1      4
	pt[0]=upp[0]; pt[1]=low[1];	pt[2]=upp[2];
	mshPts3DList.push_back(pt);//right_down z=1     5
	mshPts3DList.push_back(upp);//right_up z=1      6
	pt[0]=low[0]; pt[1]=upp[1];	pt[2]=upp[2];
	mshPts3DList.push_back(pt);//left_up z=1        7

	mshPtsRepeatedList[0]=false;
	mshPtsRepeatedList[1]=false;
	mshPtsRepeatedList[2]=false;
	mshPtsRepeatedList[3]=false;
	mshPtsRepeatedList[4]=false;
	mshPtsRepeatedList[5]=false;
}

//void cOctNode::addGeoFFacesIndx(int _indx) { geoFFacesIndxList.push_back(_indx); }
//
//void cOctNode::addGeoFPtsIndx(int _indx) { geoFPtsIndxList.push_back(_indx); }

int cOctNode::numOfGeoFFaces() { return (int)(geoFFacesList.size()); }

int cOctNode::numOfGeoFPts(){return (int)(geoFPtsList.size());}

int cOctNode::numOfGeoFEdges(){return (int)(geoFEdgesList.size());}

void cOctNode::addNode(int _level, string _nid, vector<double> _position, double _size,
						vector<cFeaturePt*> _geoFPtsList,
						vector<cFeatureEdge*> _geoFEdgesList,
						vector<cFeatureFace*> _geoFFacesList,
						cOctNode* _parent)
{
    children.push_back(
    		new cOctNode(_level, _nid, _position, _size, _geoFPtsList, _geoFEdgesList, _geoFFacesList, _parent)
    		);
}

bool cOctNode::sphereRayIntersect(cLine &ray)
{
    // Quick test for determining if a ray is *likely* to intersect a given node

    // Radius of sphere that contains node
    double radius = commonFunc::getInstance()->distBetweenPoints(low,position);

    // Project centre of sphere (node.position) onto ray
    vector<double> oc = commonFunc::getInstance()->vectSubtract(position, ray.p0); //vector - ray Origin -> node's Center 
    double s = commonFunc::getInstance()->dotProduct(oc,ray.dir); // projection length of oc onto ray
    vector<double> projpnt = commonFunc::getInstance()->vectAdd(ray.p0, ray.dir, s); //move ray origin along ray dir to the projection point
    double dist = commonFunc::getInstance()->distBetweenPoints(projpnt,position);

    // If distance between spherical centre and projected point is
    // less than the radius of the sphere, then an intersection is
    // *possible*
    return (dist<=radius);
}

bool cOctNode::boxRayIntersect(cLine &ray)
{
    // An accurate test for determining if a ray will intersect a given node.
    // Tests for intersections between the ray and all 6 faces of the node.

    vector<double> p; double D, sDen, sNum, s, tol = 1.0e-06; int i, j;

    for (int faceId=1; faceId<=6; faceId++)
    {
        // Get D (distance of plane to origin) and N (face normal) of node face
        vector<double> N(3,0.0);
        switch(faceId) {
            case 1: {
                D = -low[0]; N[0] = -1.0; break; } // -x face
            case 2: {
                D = -low[1]; N[1] = -1.0; break; } // -y face
            case 3: {
                D = -low[2]; N[2] = -1.0; break; } // -z face
            case 4: {
                D =  upp[0]; N[0] =  1.0; break; } // +x face
            case 5: {
                D =  upp[1]; N[1] =  1.0; break; } // +y face
            case 6: {
                D =  upp[2]; N[2] =  1.0; }        // +z face
        }

        // Get intersection point between face plane and ray. If no intersection is
        // possible (i.e. the normal of the face is perp. to the line) then skip face
        sDen = commonFunc::getInstance()->dotProduct(ray.dir,N);
        if (fabs(sDen)>tol) {

            // Find intersection point p
            sNum = D - commonFunc::getInstance()->dotProduct(ray.p0,N);
            s    = sNum / sDen;
            p    = commonFunc::getInstance()->vectAdd(ray.p0,ray.dir,s);

            // Check if intersection point is within bounds of face. If so, then
            // return true. If not, then skip face
            if      (faceId==1 || faceId==4) { i=1; j=2; } // -x,+x
            else if (faceId==2 || faceId==5) { i=0; j=2; } // -y,+y
            else if (faceId==3 || faceId==6) { i=0; j=1; } // -z,+z
            if ((p[i]>=low[i] && p[i]<=upp[i]) && (p[j]>=low[j] && p[j]<=upp[j])) {
                return true; }
        }
    }
    return false;
}

void cOctNode::calMshFaces()
{
    for (unsigned int i = 0; i < mshFacesList.size(); i++) {
        mshFacesList[i].low.resize(3, -100.0);
        mshFacesList[i].upp.resize(3, -100.0);
        //mshFacesList[i].ptLabelList.resize(4,-100);
        mshFacesList[i].nid="node-"+nid+"->face-"+commonFunc::getInstance()->NumberToString(i);
        mshFacesList[i].own=mshVolIndx;
        //mshFacesList[i].mshVolIndx=mshVolIndx;
        mshFacesList[i].nPts=4;
        mshFacesList[i].node=this;
    }
    mshFacesList[0].low = mshPts3DList[0]; //face0 - s;
    mshFacesList[1].low = mshPts3DList[1]; //face1 - e
    mshFacesList[2].low = mshPts3DList[3]; //face2 - n
    mshFacesList[3].low = mshPts3DList[0]; //face3 - w
    mshFacesList[4].low = mshPts3DList[4]; //face4 - t
    mshFacesList[5].low = mshPts3DList[0]; //face5 - b

    mshFacesList[0].upp = mshPts3DList[5]; //face0 - s;
    mshFacesList[1].upp = mshPts3DList[6]; //face1 - e
    mshFacesList[2].upp = mshPts3DList[6]; //face2 - n
    mshFacesList[3].upp = mshPts3DList[7]; //face3 - w
    mshFacesList[4].upp = mshPts3DList[6]; //face4 - t
    mshFacesList[5].upp = mshPts3DList[2]; //face5 - b

    mshFacesList[0].ptsList.push_back(mshPts3DList[0]);
    mshFacesList[0].ptsList.push_back(mshPts3DList[1]);
    mshFacesList[0].ptsList.push_back(mshPts3DList[5]);
    mshFacesList[0].ptsList.push_back(mshPts3DList[4]);

    mshFacesList[1].ptsList.push_back(mshPts3DList[1]);
    mshFacesList[1].ptsList.push_back(mshPts3DList[2]);
    mshFacesList[1].ptsList.push_back(mshPts3DList[6]);
    mshFacesList[1].ptsList.push_back(mshPts3DList[5]);

    mshFacesList[2].ptsList.push_back(mshPts3DList[2]);
    mshFacesList[2].ptsList.push_back(mshPts3DList[3]);
    mshFacesList[2].ptsList.push_back(mshPts3DList[7]);
    mshFacesList[2].ptsList.push_back(mshPts3DList[6]);

    mshFacesList[3].ptsList.push_back(mshPts3DList[0]);
    mshFacesList[3].ptsList.push_back(mshPts3DList[4]);
    mshFacesList[3].ptsList.push_back(mshPts3DList[7]);
    mshFacesList[3].ptsList.push_back(mshPts3DList[3]);

    mshFacesList[4].ptsList.push_back(mshPts3DList[5]);
    mshFacesList[4].ptsList.push_back(mshPts3DList[6]);
    mshFacesList[4].ptsList.push_back(mshPts3DList[7]);
    mshFacesList[4].ptsList.push_back(mshPts3DList[4]);

    mshFacesList[5].ptsList.push_back(mshPts3DList[1]);
    mshFacesList[5].ptsList.push_back(mshPts3DList[0]);
    mshFacesList[5].ptsList.push_back(mshPts3DList[3]);
    mshFacesList[5].ptsList.push_back(mshPts3DList[2]);
//----------------------------------------------------------
    update_MshFaces_PtIndxList();
}

void cOctNode::update_MshFaces_PtIndxList() //init eleFaces's ptLabelList
{
    mshFacesList[0].ptIndxList.push_back(mshPtsIndxList[0]);
    mshFacesList[0].ptIndxList.push_back(mshPtsIndxList[1]);
    mshFacesList[0].ptIndxList.push_back(mshPtsIndxList[5]);
    mshFacesList[0].ptIndxList.push_back(mshPtsIndxList[4]);

    mshFacesList[1].ptIndxList.push_back(mshPtsIndxList[1]);
    mshFacesList[1].ptIndxList.push_back(mshPtsIndxList[2]);
    mshFacesList[1].ptIndxList.push_back(mshPtsIndxList[6]);
    mshFacesList[1].ptIndxList.push_back(mshPtsIndxList[5]);

    if(this->nid=="0-1-7"){
    	cout<<"";
    }
    mshFacesList[2].ptIndxList.push_back(mshPtsIndxList[2]);
    mshFacesList[2].ptIndxList.push_back(mshPtsIndxList[3]);
    mshFacesList[2].ptIndxList.push_back(mshPtsIndxList[7]);
    mshFacesList[2].ptIndxList.push_back(mshPtsIndxList[6]);

    mshFacesList[3].ptIndxList.push_back(mshPtsIndxList[0]);
    mshFacesList[3].ptIndxList.push_back(mshPtsIndxList[4]);
    mshFacesList[3].ptIndxList.push_back(mshPtsIndxList[7]);
    mshFacesList[3].ptIndxList.push_back(mshPtsIndxList[3]);

    mshFacesList[4].ptIndxList.push_back(mshPtsIndxList[5]);
    mshFacesList[4].ptIndxList.push_back(mshPtsIndxList[6]);
    mshFacesList[4].ptIndxList.push_back(mshPtsIndxList[7]);
    mshFacesList[4].ptIndxList.push_back(mshPtsIndxList[4]);

    mshFacesList[5].ptIndxList.push_back(mshPtsIndxList[1]);
    mshFacesList[5].ptIndxList.push_back(mshPtsIndxList[0]);
    mshFacesList[5].ptIndxList.push_back(mshPtsIndxList[3]);
    mshFacesList[5].ptIndxList.push_back(mshPtsIndxList[2]);
}


//void cOctNode::addGeoFFacesIndx(int _indx) {
//	geoFFacesIndxList.push_back(_indx);
//}

void cOctNode::removeExtraFeats() {
//	cFeaturePt *fPt;
//	cFeatureEdge *fEdge;
//	cTri *fTri;

   //use a tmpList to keep the elems in the node and replace geoFPtsList by tmpList
    //remove feature point
    vector <cFeaturePt*> tmpGeoFPtsList;
    for (vector <cFeaturePt*>::iterator it = geoFPtsList.begin(); it != geoFPtsList.end(); it++) {
        if ((*it)->isInNode(this)) {
            tmpGeoFPtsList.push_back(*it);
        }
    }
    geoFPtsList = tmpGeoFPtsList;

    //remove feature edge
    vector <cFeatureEdge*> tmpGeoFEdgesList;
    for (vector <cFeatureEdge*>::iterator it = geoFEdgesList.begin(); it != geoFEdgesList.end(); it++) {
        if ((*it)->isInNode(this)) {
            tmpGeoFEdgesList.push_back(*it);
        }
    }
    geoFEdgesList = tmpGeoFEdgesList;

    //remove feature face
    vector <cFeatureFace*> tmpGeoFFacesList;
    for (vector <cFeatureFace*>::iterator it = geoFFacesList.begin(); it != geoFFacesList.end(); it++) {
        if ((*it)->isInNode(this)) {
            tmpGeoFFacesList.push_back(*it);
        }
    }
    geoFFacesList=tmpGeoFFacesList;
}

bool cOctNode::isPtInNode(vector<double> pt) {
	double tol=1e-6;
    if (low[0] - pt[0] >tol) return false;
    if (low[1] - pt[1] >tol) return false;
    if (low[2] - pt[2] >tol) return false;
    if (upp[0] - pt[0] <-tol) return false;
    if (upp[1] - pt[1] <-tol) return false;
    if (upp[2] - pt[2] <-tol) return false;
    return true;
}
