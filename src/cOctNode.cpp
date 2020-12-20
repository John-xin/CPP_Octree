#include "cOctNode.h"
#include "cFeature.h"
#include "cTri.h"
#include "cOctree.h"
#include "commonFunc.h"

cOctNode::cOctNode()
{
    // Default octNode constructor
    level = 0;
    depth = 1;
    nid   = "";
    size  = 1.0;
    isBoundaryNode=-100;
    isInteriorNode=-100;
    parent=NULL;
    position.resize(3,0.0);
    getLowUppVerts();

    //geoFFacesList.reserve(MAX_OCTNODE_OBJECTS);
}

cOctNode::cOctNode(int _level, string _nid, vector<double> _position, double _size,
		vector<cFeaturePt*> _geoFPtsList,
		vector<cFeatureEdge*> _geoFEdgesList,
		vector<cTri*> _geoFFacesList,
		cOctNode* _parent)
{
    // octNode constructor with level, node id (nid), position and size
    level    = _level;
    depth = 1;
    nid      = _nid;
    position = _position;
    size     = _size;
    mshVolIndx = -100;
    isBoundaryNode=-100;
    isInteriorNode=-100;
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
    setGeoFFacesIndxList();
    //geoFFacesIndxList.reserve(MAX_OCTNODE_OBJECTS);
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
						vector<cTri*> _geoFFacesList,
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
    vector<double> oc = commonFunc::getInstance()->vectSubtract(position, ray.p0);
    double s = commonFunc::getInstance()->dotProduct(oc,ray.dir);
    vector<double> projpnt = commonFunc::getInstance()->vectAdd(ray.p0, ray.dir, s);
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


void cOctNode::updateMshFaces_mshVolIndx(int _mshVolIndx)
{
    for (unsigned i = 0; i < mshFacesList.size(); i++) {
        mshFacesList[i].mshVolIndx = _mshVolIndx;
    }
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

void cOctNode::setGeoFFacesIndxList() {
	for(unsigned int i=0; i< geoFFacesList.size();i++){
		geoFFacesIndxList.push_back(geoFFacesList[i]->indx);
		}
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
    vector <cTri*> tmpGeoFFacesList;
    for (vector <cTri*>::iterator it = geoFFacesList.begin(); it != geoFFacesList.end(); it++) {
        if ((*it)->isInNode(this)) {
            tmpGeoFFacesList.push_back(*it);
        }
    }
    tmpGeoFFacesList = geoFFacesList;
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
