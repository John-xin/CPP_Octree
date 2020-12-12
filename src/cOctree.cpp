
// Copyright (C) 2017 Michael Hogg

// This file is part of pyoctree - See LICENSE.txt for information on usage and redistribution

#include "cOctree.h"
#include "cOctNode.h"
#include "commonFunc.h"
#include <algorithm> // find, sort
#include <iomanip> // set decimals with 3 digits

//#include "cTri.h"
//#include "cFeature.h"
//#include "cGeomData.h"
//#include "cOctNode.h"

cOctree::cOctree()
{
	//octree setup depends on geomData - run geomData read STL first
    numOfGeoPts=cGeomData::getInstance()->numOfPts;
    numOfGeoTris=cGeomData::getInstance()->numOfTris;
    numOfGeoPhyNames=cGeomData::getInstance()->numOfPhyNames;
    ptInGeom.resize(3);

    //extract features from geom data
    setup_geoFFacesList(); //first step
    setup_geoFEdgesList(); //find feature edge from setup_geoFFacesList
    setup_geoFPtsList();  //find feature point from setup_geoFEdgesList
    cout<<"\n";
    //++++++++++++++++++++++++++++++++++++++
    int _offsets[][3] = {{-1,-1,-1},{+1,-1,-1},{-1,+1,-1},{+1,+1,-1},
                         {-1,-1,+1},{+1,-1,+1},{-1,+1,+1},{+1,+1,+1}};

    for (int i=0; i<8; i++) {
        for (int j=0; j<3; j++) {
            branchOffsets[i][j] = _offsets[i][j];
        }
    }

    //buildOctree();
    //saveAsOFMesh();
}
cOctree::~cOctree()
{
    cout << "Destroying the cOctree" << endl;
}
int cOctree::countAllNodes(cOctNode& node)
{
    int numOfNodes = 0;

    if (node.isLeafNode()) {
        numOfNodes++;
        return numOfNodes;
    }
    else {
        numOfNodes++;
        for (unsigned int i = 0; i < node.branches.size();i++) {
            numOfNodes+=countAllNodes(node.branches[i]);

        }
        return numOfNodes;

    }
}
int cOctree::countIntlNodes(cOctNode& node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<nonBNodesList.size();i++){
    	if(nonBNodesList[i]->state==1){
    		numOfNodes++;
    	}
    }
    return numOfNodes;
}
int cOctree::countBNodes(cOctNode& node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<bNodesList.size();i++){
    		numOfNodes++;
    }
    return numOfNodes;
}
int cOctree::countExtNodes(cOctNode& node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<nonBNodesList.size();i++){
    	if(nonBNodesList[i]->state==-1){
    		numOfNodes++;
    	}
    }
    return numOfNodes;
}
void cOctree::outputNodeName(cOctNode& node) {
	    if (node.isLeafNode()) {
	    	cout<< node.nid.c_str() <<"\n";
	    }
	    else {
	        for (unsigned int i = 0; i < node.branches.size();i++) {
	        	outputNodeName(node.branches[i]);
	        }
	    }
}
void cOctree::outputMshPts(cOctNode& node) {
    vector<int> ptIndxList;
	if (node.isLeafNode()) {
		ptIndxList=node.mshPtsIndxList;
		cout<< node.nid.c_str() <<" state:"<<node.state<<" ptIndxList (" << ptIndxList[0]<<"," << ptIndxList[1]<<","<< ptIndxList[2]<<","<< ptIndxList[3]<<","<< ptIndxList[4]<<","<< ptIndxList[5]<<","<< ptIndxList[6]<<","<< ptIndxList[7]<<")"<<"\n";

    }
    else {
        for (unsigned int i = 0; i < node.branches.size();i++) {
        	outputMshPts(node.branches[i]);
        }
    }
}
void cOctree::outputMshFaces(cOctNode& node) {
    vector<cFace> faceList;
    cFace face;
	if (node.isLeafNode()) {
		faceList=node.mshFacesList;
		for(unsigned j=0; j<faceList.size();j++){
			face=faceList[j];
			cout<< face.nid.c_str() <<" state:"<<face.state<<" ptIndxList (" << face.ptIndxList[0]<<"," << face.ptIndxList[1]<<","<< face.ptIndxList[2]<<","<< face.ptIndxList[3]<<")"<<"\n";
		}


    }
    else {
        for (unsigned int i = 0; i < node.branches.size();i++) {
        	outputMshFaces(node.branches[i]);
        }
    }
}

// ++++++ setup feature list
void cOctree::setup_geoFFacesList()
{
    vector<int> geoPtIndx_int3(3);
    int phyNameIndx;
    int geoFFaceIndx;
    geoFFacesList.reserve(numOfGeoTris);
    for (int i=0; i<numOfGeoTris; i++) {
        geoPtIndx_int3[0]=geoTrisList[i][0];
        geoPtIndx_int3[1]=geoTrisList[i][1];
        geoPtIndx_int3[2]=geoTrisList[i][2];
        phyNameIndx=geoTrisList[i][3];
        geoFFaceIndx=i;
        cTri* fFace=new cTri(geoFFaceIndx,phyNameIndx,geoPtIndx_int3);
        geoFFacesList.push_back(fFace);
    }
    std::cout<<"number of geoFFaces is " <<geoFFacesList.size()<<"\n";
}
void cOctree::setup_geoFEdgesList()
{
    vector<int> geoPtIndx_int2;
    int geoFEdgeIndx=0;
    cTri* surface1;
	cTri* surface2;

    unsigned int numOfSurfaces=geoFFacesList.size();
    int k=0;
    //this needs to be optimized - cost too much time
    for (unsigned int i=0; i<numOfSurfaces; i++) {
        for (unsigned int j=k; j<numOfSurfaces; j++) {
            surface1=geoFFacesList[i];
            surface2=geoFFacesList[j];
            if (i!=j) {
                if (surface1->isAdjacent(surface2)) {
                    if (surface1->getAngle(surface2) < 170) {
                        geoPtIndx_int2 = surface1->getAdjacentEdge(surface2);
                        cFeatureEdge* fEdge=new cFeatureEdge(geoFEdgeIndx, geoPtIndx_int2);
                        geoFEdgesList.push_back(fEdge);
                        geoFEdgeIndx = geoFEdgeIndx + 1;
                    }
                }
            }

        }
        k=k+1;
    }
    std::cout<<"number of feature edges is " <<geoFEdgesList.size()<<"\n";
}
void cOctree::setup_geoFPtsList()
{
    cFeaturePt* fPt=new cFeaturePt();
    for (unsigned int i=0; i<geoFEdgesList.size(); i++) {
    	for (int j=0; j<2; j++) {
    		fPt->indx=geoFEdgesList[i]->ptIndx_int2[j];
    		if(fPt->isNewInList(geoFPtsList)){
    			geoFPtsList.push_back(fPt);
    		}
        }
    }
    std::cout<<"number of feature points is " <<geoFPtsList.size()<<"\n";
}
int cOctree::numFFaces() { return (int)(geoFFacesList.size()); }
// ++++++ setup feature list

void cOctree::defineBody(vector<double> _ptInGeom) {
	if (_ptInGeom.size()!=3){
		cout<<"Must give a pt3D to define body"<<"\n";
	}else{
		ptInGeom=_ptInGeom;
	}

}

void cOctree::buildOctree() {
    //1. init root node
    setup_root();

    //2. +++++++++++++++++++++++++++
	std::cout << "Min Octree Level is " << MIN_OCTREE_LEVELS << "\n";
    splitOctreeByMinLevel(root);
    //std::cout << "After splitOctreeByMinLevel " << "\n";
    //outputNodeName(root);

	std::cout<<"Max FeatPt in node is "<< MAX_OCTNODE_FEATS<<"\n\n";
    splitOctreeByFeaturePt(root);
    //splitNodeById("0-0");
    setup_leafNodesList(&root);
    //std::cout << "After splitOctreeByFeaturePt " << "\n";
    //outputNodeName(root);
    //cout<<"\n";

    //3. identify boundary / interior / exterior node +++++++++++++++++++++++++++
    setup_bNodesList(&root); //leaf nodes -> identify boundary / non-boundary node
    std::cout<<"bNodesList length is "<< bNodesList.size() <<"\n";
    std::cout<<"nonBNodesList length is "<< nonBNodesList.size() <<"\n";
    setup_leafNodesNbr();//n*logN
    setup_nodesState(); //non-boundary node -> identify interior / exterior node
    //delExtNodes(root);
    cOctNode* aTMPNode= getNodeFromId("0-0-3-3-3-2");
    cout<<aTMPNode->state;


    //4. +++++++++++++++++++++++++++
    //from leaf nodes -> identify non-repeated mshPt

    setup_mshPtList();
    cout<<"num of mshPts is "<< mshPtsList.size() <<"\n";
    //outputMshPts(root);
    //from leaf nodes -> define node -> 6 faces by mshPt index
    //init mshFace -> own, nid, low, upp, ptIndxList
    //and mark mshVolIndx
    setup_nodeMshFaces();


    //5. from boundary / interior node -> identify internal / boundary mshFace +++++++++++++++++++++++++++
    //update own, nei, ptIndxList
    setup_mshFaceList();
    std::cout<<"mshIntlFacesList length is "<< mshIntlFacesList.size() <<"\n";
    std::cout<<"mshBFacesList length is "<< mshBFacesList.size() <<"\n\n";
    //outputMshFaces(root);
    //6. from boundary mshFace -> identify phyName of boundary mshFace +++++++++++++++++++++++++++
    setup_nodePhyName();
    setup_mshBdsList();


    //
    num_allNodes=countAllNodes(root);
    num_intlNodes=countIntlNodes(root);
    num_extNodes=countExtNodes(root);
    num_bNodes=countBNodes(root);
    std::cout<<"num of Nodes(root->all leaves) is "<< num_allNodes <<"\n";
    std::cout<<"num of leaf Nodes is "<< leafNodesList.size() <<"\n";
    std::cout<<"num of Interior Nodes is "<< num_intlNodes <<"\n";
    std::cout<<"num of Exterior Nodes is "<< num_extNodes <<"\n";
    std::cout<<"num of Boundary Nodes is "<< num_bNodes <<"\n\n";
}

// +++++ setup root
void cOctree::setup_root() {
    vector<double> position = getPositionRoot();
    double size = getSizeRoot();
//    vector<int> geoFPtsIndxList;
//    vector<int> geoEdgesIndxList;
//    vector<int> geoFFacesIndxList;
//    for(unsigned i=0;i<geoFPtsList.size();i++){
//    	geoFPtsIndxList.push_back(geoFPtsList[i].indx);
//    }
//    for(unsigned i=0;i<geoFEdgesList.size();i++){
//    	geoFPtsIndxList.push_back(geoFEdgesList[i].indx);
//    }
//    for(unsigned i=0;i<geoFFacesList.size();i++){
//    	geoFPtsIndxList.push_back(geoFEdgesList[i].indx);
//    }
    root = cOctNode(0,"0", position, size,geoFPtsList,geoFEdgesList,geoFFacesList, NULL);
}
double cOctree::getSizeRoot() {

    // Get low and upp
    vector<double> low, upp, range;
    low = cGeomData::getInstance()->low;
    upp = cGeomData::getInstance()->upp;

    // Range is the size of the node in each coord direction
    range = commonFunc::getInstance()-> vectSubtract(upp,low);
    double size = range[0];
    for (int i=1; i<3; i++) {
        if (range[i] > size) { size = range[i]; }
    }
    // Scale up size of node by 5%
    size *= 1.01;
    return size;
}
vector<double> cOctree::getPositionRoot() {

    // Get low and upp
    vector<double> low, upp, position(3);
    low = cGeomData::getInstance()->low;
    upp = cGeomData::getInstance()->upp;
    // Center of node is average of low and upp
    for (int i=0; i<3; i++) {
        position[i] = 0.5 * (low[i]+upp[i]);
    }
    return position;
}
// +++++ setup root

// +++++++++ split node
void cOctree::splitOctreeByFeaturePt(cOctNode& node) {

	if(node.isLeafNode()){
		splitNodeByFeaturePt(node);
	} else{
		for(int i=0; i<node.NUM_BRANCHES_OCTNODE; i++){
			splitOctreeByFeaturePt(node.branches[i]);
		}
	}
}
void cOctree::splitNodeByFeaturePt(cOctNode &node)
{
	if (node.numOfGeoFPts()> MAX_OCTNODE_FEATS && node.level<MAX_OCTREE_LEVELS){
	        // Split node into 8 branches
	        vector<double> newPos;
	        newPos.resize(3);
	        for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
	            for(int j=0; j<3; j++) {
	                newPos[j] = node.position[j] + 0.25*node.size*branchOffsets[i][j];
	            }
	            string nid = node.nid + "-" + commonFunc::getInstance()->NumberToString(i);
	            node.addNode(node.level+1,nid,newPos,0.5*node.size,
	            			 node.geoFPtsList,
	            			 node.geoFEdgesList,
	            			 node.geoFFacesList,
	            			 &node);
	        }
	        // Reallocate date from node to branches
	        for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
	            if (node.branches[i].numOfGeoFPts()> MAX_OCTNODE_FEATS  && node.branches[i].level<MAX_OCTREE_LEVELS)
	            {
	            	splitNodeByFeaturePt(node.branches[i]);
	            }
	        }

	        node.geoFPtsList.resize(0);
	        node.geoFEdgesList.resize(0);
	        node.geoFFacesList.resize(0);
	        node.geoFFacesIndxList.resize(0);

	    }
}
void cOctree::splitOctreeByMinLevel(cOctNode& node) {

	if(node.isLeafNode()){
		splitNodeByMinLevel(node);
	} else{
		for(int i=0; i<node.NUM_BRANCHES_OCTNODE; i++){
			splitOctreeByMinLevel(node.branches[i]);
		}
	}
}
void cOctree::splitNodeByMinLevel(cOctNode &node)
{
    if (node.level<MIN_OCTREE_LEVELS){
        // Split node into 8 branches
        vector<double> newPos;
        newPos.resize(3);
        for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
            for(int j=0; j<3; j++) {
                newPos[j] = node.position[j] + 0.25*node.size*branchOffsets[i][j];
            }
            string nid = node.nid + "-" + commonFunc::getInstance()->NumberToString(i);
            node.addNode(node.level+1,nid,newPos,0.5*node.size,
            			 node.geoFPtsList,
            			 node.geoFEdgesList,
            			 node.geoFFacesList,
            			 &node);
        }
        // split child node until min level
        for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
            if (node.branches[i].level < MIN_OCTREE_LEVELS) {
                    splitNodeByMinLevel(node.branches[i]);
                }
        }

        node.geoFPtsList.resize(0);
        node.geoFEdgesList.resize(0);
        node.geoFFacesList.resize(0);
        node.geoFFacesIndxList.resize(0);
    }

}
void cOctree::splitNodeById(string node_id) {
	cOctNode* node=getNodeFromId(node_id);

        // Split node into 8 branches
        vector<double> newPos;
        newPos.resize(3);
        for (int i=0; i<node->NUM_BRANCHES_OCTNODE; i++) {
            for(int j=0; j<3; j++) {
                newPos[j] = node->position[j] + 0.25*node->size*branchOffsets[i][j];
            }
            string nid = node->nid + "-" + commonFunc::getInstance()->NumberToString(i);
            node->addNode(node->level+1,nid,newPos,0.5*node->size,
            			 node->geoFPtsList,
            			 node->geoFEdgesList,
            			 node->geoFFacesList,
            			 node);
        }

        node->geoFPtsList.resize(0);
        node->geoFEdgesList.resize(0);
        node->geoFFacesList.resize(0);
        node->geoFFacesIndxList.resize(0);

}

void cOctree::setup_leafNodesList(cOctNode* node) {
    if (node->isLeafNode()) {
    	leafNodesList.push_back(node);
    }
    else {
        for (unsigned int i = 0; i < node->branches.size(); i++) {
        	setup_leafNodesList(&(node->branches[i]));
        }
    }
}
void cOctree::setup_leafNodesNbr() {
    //assume node level in the tree is less than one

	//need to return list of nodes for each nbr, e.g. low-level node nbr has 4  high-level nodes
	//need to modify getLeafNodeByPt() to return list of nodes

	for(unsigned i=0;i<leafNodesList.size();i++){
		cOctNode* node=leafNodesList[i];
		node->nbr.resize(6);
		for(unsigned j=0;j<6;j++){
			node->nbr[j].resize(0);
		}

	    	vector<double> pt;
	    	vector<vector <double> > vect(6);
	    	vect[0].resize(3);//nbr s
	    	vect[0][0]=0;
	    	vect[0][1]=-0.75*node->size;
	    	vect[0][2]=0;

	    	vect[1].resize(3);//nbr e
	    	vect[1][0]=0.75*node->size;
	    	vect[1][1]=0;
	    	vect[1][2]=0;

	    	vect[2].resize(3);//nbr n
	    	vect[2][0]=0;
	    	vect[2][1]=0.75*node->size;
	    	vect[2][2]=0;

	    	vect[3].resize(3);//nbr w
	    	vect[3][0]=-0.75*node->size;
	    	vect[3][1]=0;
	    	vect[3][2]=0;

	    	vect[4].resize(3);//nbr t
	    	vect[4][0]=0;
	    	vect[4][1]=0;
	    	vect[4][2]=0.75*node->size;

	    	vect[5].resize(3);//nbr b
	    	vect[5][0]=0;
	    	vect[5][1]=0;
	    	vect[5][2]=-0.75*node->size;

	    	for(int j=0;j<6;j++){
	    		pt=node->position;
	        	pt=commonFunc::getInstance()->vectAdd(pt,vect[j]);
	        	node->nbr[j]=getLeafNodeByPt(pt,root);

	        	if(node->nbr[j].size()==0){
	        			//cout<<"leaf node" << node->nid.c_str() <<"->nbr-" << j << " Null" <<"\n";

	        	}else{
	        		vector <cOctNode*>::iterator it = node->nbr[j].begin();
	        		while (it != node->nbr[j].end())
	        		{
	        		      //cout<<"leaf node" << node->nid.c_str() <<"->nbr-" << j << " "<<(*it)->nid.c_str() <<"\n";
	        		      ++it;
	        		}
	        	}
	        }
    }
}
// +++++++++ split node

//vector<cOctNode*> cOctree::getNodesFromLabel(int polyLabel)
//{
//    // Function for finding all the nodes that contains tri with given label
//    vector<cOctNode*> nodeList;
//    findBranchesByLabel(polyLabel,root,nodeList);
//    return nodeList;
//}

//void cOctree::findBranchesByLabel(int polyLabel, cOctNode &node, vector<cOctNode*> &nodeList)
//{
//    // Recursive function used by getNodesFromLabel
//    if (node.isLeafNode()) {
//        vector<int>::iterator it;
//        it = find(node.geoFFacesList.begin(),node.geoFFacesList.end(),polyLabel);
//        if (it != node.geoFFacesList.end()) { nodeList.push_back(&node); }
//    } else {
//        for (unsigned int i=0; i<node.branches.size(); i++) {
//            findBranchesByLabel(polyLabel, node.branches[i], nodeList);
//        }
//    }
//}

vector<cOctNode*> cOctree::getLeafNodeByPt(vector<double> pt, cOctNode& node) {
	tmpNodesList.resize(0);
	addLeafNodeByPt(pt, node);
	return tmpNodesList;
}
void cOctree::addLeafNodeByPt(vector<double> pt, cOctNode& node) {
	if(node.isPtInNode(pt)){
	    if (!node.isLeafNode()) {
	        for (unsigned int i = 0; i < node.branches.size();i++) {
	        	addLeafNodeByPt(pt, node.branches[i]);
	        }
	    }
	    else {
	    	tmpNodesList.push_back(&node);
	    }
	}
}
cOctNode* cOctree::getNodeFromId(string nodeId)
{
    return findBranchById(nodeId,root);
}

cOctNode* cOctree::findBranchById(string nodeId, cOctNode &node)
{
    if (nodeId.compare(node.nid)==0) {
        return &node;
    } else {
        for (unsigned int i=0; i<node.branches.size(); i++) {
            cOctNode *branch = findBranchById(nodeId, node.branches[i]);
            if (branch != NULL) { return branch; }
        }
    }
    return NULL;
}

//+++++++++++++ Intersection
vector<Intersection> cOctree::findRayIntersect(cLine &ray)
{
    // Get polys to check
    set<int> polyListCheck = getListPolysToCheck(ray);

    // Loop through all polys in check list to find a possible intersection
    vector<Intersection> intersectList;
    set<int>::iterator it;
    vector<double> ip;
    double s;
    for (it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
        int polyLabel = *it;
        if (geoFFacesList[polyLabel]->rayPlaneIntersectPoint(ray,ip,s)) {
            intersectList.push_back(Intersection(polyLabel,ip,s)); }
    }

    // Sort list in terms of distance of the intersection from the ray origin
    sort(intersectList.begin(),intersectList.end());

    return intersectList;
}
set<int> cOctree::getListPolysToCheck(cLine &ray)
{
    // Returns a list of all polygons that are within OctNodes hit by a given ray
    set<int> intTestPolys;
    getPolysToCheck(root,ray,intTestPolys);
    return intTestPolys;
}
void cOctree::getPolysToCheck(cOctNode &node, cLine &ray, set<int> &intTestPolys)
{
    // Utility function for getListPolysToCheck. Finds all OctNodes hit by a given ray
    // and returns a list of the objects contained within
    if (node.sphereRayIntersect(ray)) {
        if (node.boxRayIntersect(ray)) {
            if (node.isLeafNode()) {
                for (int i=0; i<node.numOfGeoFFaces(); i++) {
                    intTestPolys.insert(node.geoFFacesIndxList[i]); }
            } else {
                for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
                    getPolysToCheck(node.branches[i],ray,intTestPolys);
                }
            }
        }
    }
}
//vector<cOctNode*> cOctree::getSortedNodesToCheck(cLine &ray)
//{
//    // Finds all the nodes that intersect with given ray. Uses the nodes "position"
//    // to sort the nodes by distance from the ray origin (in ascending order).
//    // Nodes that are closest to the ray origin will be checked first for poly
//    // intersections
//    vector<pair<cOctNode*,double> > nodeList;
//    getNodesToCheck(root,ray,nodeList);
//    sort(nodeList.begin(),nodeList.end(),commonFunc::getInstance()->sortNodes);
//    vector<cOctNode*> nodeListSorted;
//    vector<pair<cOctNode*,double> >::iterator it;
//    for (it=nodeList.begin(); it!=nodeList.end(); it++) {
//        nodeListSorted.push_back((*it).first); }
//    return nodeListSorted;
//}

//void cOctree::getNodesToCheck(cOctNode &node, cLine &ray, vector<pair<cOctNode*,double> > &nodeList)
//{
//    // Utility function for getSortedNodesToCheck
//    // Finds all the nodes that intersect with given ray. Projects the node "position" (node
//    // centre) onto the ray to facilitate sorting of the nodes by distance from ray origin
//    if (node.sphereRayIntersect(ray)) {
//        if (node.boxRayIntersect(ray)) {
//            if (node.isLeafNode()) {
//                // Project node "position" on to ray and find distance from ray origin
//                vector<double> oc = commonFunc::getInstance()->vectSubtract(node.position, ray.p0);
//                double s = commonFunc::getInstance()->dotProduct(oc,ray.dir);
//                // Add node and distance to list
//                //nodeList.push_back(std::make_pair(&node,s));
//            } else {
//                for (int i=0; i<node.NUM_BRANCHES_OCTNODE; i++) {
//                    getNodesToCheck(node.branches[i],ray,nodeList);
//                }
//            }
//        }
//    }
//}

//vector<int> cOctree::findRayIntersects(vector<cLine> &rayList)
//{
//    // For each ray provided, determines if ray hits a poly in the tree and
//    // returns a boolean integer. Uses openmp to speed up the calculation
//    // Function findRayIntersectsSorted is a similar function that sorts the
//    // triangles in order of closest octNodes. For ray casting, this alternative
//    // function should be faster in most cases
//
//    int numRays = (int)(rayList.size());
//    vector<int> foundIntsects(numRays,0);
//    #pragma omp parallel for
//    for (int i=0; i<numRays; i++)
//    {
//        cLine ray = rayList[i];
//        set<int> polyListCheck = getListPolysToCheck(ray);
//        for (set<int>::iterator it=polyListCheck.begin(); it!=polyListCheck.end(); ++it) {
//            int polyLabel = *it;
//            //if (geoFFacesList[polyLabel].rayPlaneIntersectPoint(ray)) {
//                foundIntsects[i] = 1; break;
//        //}
//        }
//    }
//    return foundIntsects;
//}

//vector<int> cOctree::findRayIntersectsSorted(vector<cLine> &rayList)
//{
//    // For each ray provided, determines if ray hits a poly in the tree and
//    // returns a boolean integer.
//    // Uses "getSortedNodesToCheck", which returns a list of nodes that intersect
//    // with the given ray sorted in ascending order of distance from the ray origin
//    // Uses openmp to speed up the calculation
//
//    int numRays = (int)(rayList.size());
//    vector<int> foundIntsects(numRays,0);
//
//    #pragma omp parallel for
//    for (int i=0; i<numRays; i++)
//    {
//        // Get branches to check. Branches are sorted in ascending distance
//        // from ray origin
//        cLine *ray = &rayList[i];
//        vector<cOctNode*> nodeList = getSortedNodesToCheck(*ray);
//
//        // Loop through sorted branches, checking the polys contained within each
//        vector<cOctNode*>::iterator it;
//        for (it=nodeList.begin(); it!=nodeList.end(); ++it) {
//            cOctNode *node = *it;
//            for (unsigned int j=0; j<node->geoFFacesIndxList.size(); j++) {
//                int polyLabel = node->geoFFacesIndxList[j];
//                //if (geoFFacesList[polyLabel].rayPlaneIntersectPoint(*ray)) {
//                    foundIntsects[i]=1; break;
//            	//}
//            }
//            // If any poly from current node is hit, proceed on to the next node
//            if (foundIntsects[i]==1) break;
//        }
//    }
//    return foundIntsects;
//}

//++++++++++++++ Intersection


void cOctree::delExtNodes(cOctNode &node)
{
	if(!node.isLeafNode()){
		for (unsigned int i = 0; i < node.branches.size();i++) {
			if(node.branches[i].state==-1){
				node.branches.erase(node.branches.begin()+i);

			}
			delExtNodes(node.branches[i]);
		}
	}
}

void cOctree::addGapNodes()
{
    //find boundary face of nodes -> single face

	//find normal of vertices of boundary face

	//find intersection point between normal and tri

	//add new node -> 4 vertices and 4 intersection points
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctree::setup_mshPtList()
{
	    for(unsigned i=0;i<leafNodesList.size();i++){
	    	addMshPtsOfNode(leafNodesList[i]);

	    }
}
void cOctree::addMshPtsOfNode(cOctNode* node) //init node.elePtLabels
{
	//need to optimize efficiency -> only loop adjacent nodes to unique

//	int mshPtIndx;
//    for(unsigned int j=0; j<node->mshPts3DList.size();j++){
//    	mshPtIndx=-1;
//
//        for (unsigned int i = 0; i < mshPtsList.size(); i++) {
//            if (isPtSame(node->mshPts3DList[j],mshPtsList[i])) {
//            	mshPtIndx = i;//old pt
//            }
//        }
//
//        if (mshPtIndx == -1) {//new pt -> elePtList
//        	mshPtsList.push_back(node->mshPts3DList[j]);
//            node->mshPtsIndxList[j]=mshPtIndxCount;
//            mshPtIndxCount++;
//            //cout<<"pt-"<< j << " indx is " << node->mshPtsIndxList[j] << "\n";
//        }else
//        //not new pt -> find index in elePtList
//        //add pt label -> node.elePtLabels
//        {
//        	node->mshPtsIndxList[j]=mshPtIndx;
//        	//cout<<"pt-"<< j << " indx is " << node->mshPtsIndxList[j] << "\n";
//        }
//    }

	vector<double> pt;
	vector<cOctNode*> nodesList;
    for(unsigned int i=0; i<node->mshPts3DList.size();i++){
    	pt=node->mshPts3DList[i];
    	if(node->mshPtsRepeatedList[i]!=true){
    		nodesList=getLeafNodeByPt(pt, root);
    		if(nodesList.size()==0){
    			cout<<"this point (" << pt[0]<< "," <<pt[1]<< "," << pt[2]<< "," << ") is not in the octree\n";

    		}else if(nodesList.size()==1){
    			// this is a non-repeated pt
    			mshPtsList.push_back(pt);
    			node->mshPtsIndxList[i]=mshPtIndxCount;
    			mshPtIndxCount++;

    		}else{
    			mshPtsList.push_back(pt);
    			node->mshPtsIndxList[i]=mshPtIndxCount;
    			for(unsigned j=0;j<nodesList.size();j++){
        			if(nodesList[j]->nid!=node->nid){
        				for(unsigned k=0; k< nodesList[j]->mshPts3DList.size();k++){
        		            if (isPtSame(pt,nodesList[j]->mshPts3DList[k])) {
        		            	nodesList[j]->mshPtsRepeatedList[k]=true;//old pt
        		            	nodesList[j]->mshPtsIndxList[k]=mshPtIndxCount;
        		            }
        				}
        			}
    			}
    			mshPtIndxCount++;
    		}
    	}
    }

    //update node's mshFaces->ptIndxList
    //node.update_MshFaces_PtIndxList(node.mshPtsIndxList);

}
bool cOctree::isPtSame(vector<double>& pt1, vector<double>& pt2)
{
    double tol = 0.000001;
    if ((abs(pt1[0] - pt2[0])) > tol) return false;
    if ((abs(pt1[1] - pt2[1])) > tol) return false;
    if ((abs(pt1[2] - pt2[2])) > tol) return false;
    return true;
}
void cOctree::setup_nodeMshFaces() {
    for(unsigned i=0;i<leafNodesList.size();i++){
    	//mark mshVolIndx to all leaf node
    	if(leafNodesList[i]->state!=-1){
        	leafNodesList[i]->mshVolIndx=mshVolIndxCount;
        	leafNodesList[i]->updateMshFaces_mshVolIndx(mshVolIndxCount);//mark all node's mshFaces
            leafNodesList[i]->calMshFaces();
            mshVolIndxCount++;
    	}else{
    		leafNodesList[i]->mshVolIndx=-1;
    		leafNodesList[i]->updateMshFaces_mshVolIndx(-1);//mark all node's mshFaces
    	}

    }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctree::setup_bNodesList(cOctNode* node) {
    for(unsigned i=0;i<leafNodesList.size();i++){
    	findBNode(leafNodesList[i]);
    }
}
void cOctree::findBNode(cOctNode* node) {
	if(node->geoFFacesList.size()!=0){
		node->state=0; //boundary node
		bNodesList.push_back(node);
	}else{
		nonBNodesList.push_back(node);
	}
}
void cOctree::findNodeState(cOctNode* node) {
    //retrieve a node from nonBNodeslist
	//make a ray -> node position to defined point
	//find num of intersections
	//even -> ext node
	//odd -> intl node

	vector<Intersection> ints;
	int numOfInts;


		cLine ray(node->position,ptInGeom,0);
		ints=findRayIntersect(ray);
		numOfInts=0;
		//锟斤拷偶锟斤拷锟斤拷锟叫憋拷
		for(unsigned i=0;i<ints.size();i++){
			if(ints[i].s>0) {numOfInts++;}
		}

		if((numOfInts%2)==0){
			node->state=-1; //1 - exterior node
		}else{
			node->state=1; //1 - internal node
			for(unsigned j=0; j<node->mshFacesList.size();j++){
				node->mshFacesList[j].state=1;
			}
		}
}
void cOctree::setup_nodesState() {
	cOctNode *node;
	node=nonBNodesList[0];
	findNodeState(node);
	int numOfNodes=0;

	if (node->state==1){
		setup_nbrNodesState(node);
		for(unsigned i=0; i<nonBNodesList.size();i++){
			if(nonBNodesList[i]->state!=1){
				nonBNodesList[i]->state=-1;
				numOfNodes++;
			}
		}
	}else if (node->state==-1){
		setup_nbrNodesState(node);
		for(unsigned i=0; i<nonBNodesList.size();i++){
			if(nonBNodesList[i]->state!=-1){
				nonBNodesList[i]->state=1;
				numOfNodes++;
			}
		}
	}

	if(node->state==-1){
		cout<<"number of Interior nodes is "<<numOfNodes<<"\n";
	}else{
		cout<<"number of Exterior nodes is "<<numOfNodes<<"\n\n";
	}
}

void cOctree::setup_nbrNodesState(cOctNode* node) {
	cOctNode* nbr;
	for(int i=0; i<6; i++){
		if(node->nbr[i].size()!=0){
			for(unsigned j=0;j<node->nbr[i].size();j++){
				nbr=node->nbr[i][j];
				if(nbr->state!=node->state && nbr->state!=0){
					nbr->state=node->state;
					setup_nbrNodesState(nbr);
				}
			}
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void cOctree::setup_mshFaceList()
{
	cFace* currentMshFace;
	cFace* nbrMshFace;
	cOctNode* node;
	vector<cOctNode*> nbr;
    string curFaceName;
    string nbrMshFaceName;
    int indxArr[6]={2, 3, 0, 1, 5, 4};

	for(unsigned int i = 0; i < bNodesList.size(); i++){
		node=bNodesList[i];
		for(unsigned j=0; j< node->mshFacesList.size();j++ ){
			mshAllFacesList.push_back(&(node->mshFacesList[j]));

			currentMshFace=&node->mshFacesList[j];
			curFaceName=currentMshFace->nid;
			if(currentMshFace->state==-100){
				nbr=node->nbr[j];
				if(nbr.size()!=0){
					unsigned int count=0;
					for(unsigned k=0;k<nbr.size();k++){
						if(nbr[k]->state!=-1){
							nbrMshFace=(&nbr[k]->mshFacesList[indxArr[j]]);
							nbrMshFaceName=nbrMshFace->nid;
							int flag=currentMshFace->findFaceRelationship(currentMshFace, nbrMshFace);
							put2List(flag,currentMshFace,nbrMshFace);
						}else{
							count++;
						}
					}
					if(count==nbr.size()){
						currentMshFace->state=0;
						mshBFacesList.push_back(currentMshFace);
					}
				}else{
					currentMshFace->state=0;
					mshBFacesList.push_back(currentMshFace);
				}
			}
		}
	}

	for(unsigned int i = 0; i < nonBNodesList.size(); i++){
		node=nonBNodesList[i];
		if(node->state==1){
			for(unsigned j=0; j< node->mshFacesList.size();j++ ){
				mshAllFacesList.push_back(&(node->mshFacesList[j]));

				currentMshFace=&node->mshFacesList[j];
				curFaceName=currentMshFace->nid;

				if(currentMshFace->state==-100){

				nbr=node->nbr[j];
				if(nbr.size()!=0){
					for(unsigned k=0;k<nbr.size();k++){
						nbrMshFace=(&nbr[k]->mshFacesList[indxArr[j]]);
						nbrMshFaceName=nbrMshFace->nid;
						int flag=currentMshFace->findFaceRelationship(currentMshFace, nbrMshFace);
						put2List(flag,currentMshFace,nbrMshFace);
					}
				}else{
		    		currentMshFace->state=0;
		    		mshBFacesList.push_back(currentMshFace);
				}
				}
			}
		}
	}
	//鈥斺��----------------------------------------------------------------

	//findMshIntlFaces();

//////////////////////////////////////////////////////////////////

    num_mshBFaces=0;
    num_mshIntlFaces=mshIntlFacesList.size();

//	for(unsigned int i = 0; i < bNodesList.size(); i++){
//			for(unsigned j=0; j< bNodesList[i]->mshFacesList.size();j++ ){
//				if(bNodesList[i]->mshFacesList[j].state==0){
//					num_mshBFaces++;
//					//cout<< bNodesList[i]->mshFacesList[j].nid << " " << bNodesList[i]->mshFacesList[j].state <<"\n";
//				}else{
//					num_mshIntlFaces++;
//				}
//			}
//	}
//
//
//	for(unsigned int i = 0; i < nonBNodesList.size(); i++){
//		if(nonBNodesList[i]->state==1){
//			for(unsigned j=0; j< nonBNodesList[i]->mshFacesList.size();j++ ){
//				if(nonBNodesList[i]->mshFacesList[j].state==0){
//					num_mshBFaces++;
//					//cout<< nonBNodesList[i]->mshFacesList[j].nid << " " << nonBNodesList[i]->mshFacesList[j].state <<"\n";
//				}else{
//					num_mshIntlFaces++;
//				}
//			}
//		}
//	}

    int num_IntlFaces_inMshBFacesList=0;
    for (unsigned int i = 0; i < mshBFacesList.size(); i++) {
    	if(mshBFacesList[i]->state==1){
    		num_IntlFaces_inMshBFacesList++;
    	}else{
    		num_mshBFaces++;
    	}
    }


    cout<<"Number of Internal Faces in mshBFacesList is "<<num_IntlFaces_inMshBFacesList<<"\n";
    cout<<"Number of Boundary Faces in mshBFacesList is "<<num_mshBFaces<<"\n\n";

}
void cOctree::findMshIntlFaces()
{
    //add face to internalEleFacesList
   //add face to boundaryEleFacesList - phycial name 1
   //add face to boundaryEleFacesList - phycial name 2
   //add face to boundaryEleFacesList - phycial name 3

    int flag=-1;
    int bFlag=0;
    cFace* currentMshFace;
    cFace* listMshFace;

    string curFaceName;
    string listFaceName;

    int numOfFaces;
    numOfFaces=mshAllFacesList.size();
    int k=0;

    for ( int i = 0; i < numOfFaces; i++) {
    	k++;
    	currentMshFace=mshAllFacesList[i];
    	curFaceName=currentMshFace->nid;
    	if(currentMshFace->state!=1){
        	for (int j=k;j<numOfFaces;j++){
        		listMshFace=mshAllFacesList[j];
        		listFaceName=listMshFace->nid;
        		//if(listMshFace->state!=1){
        			flag=currentMshFace->findFaceRelationship(currentMshFace, listMshFace);
        			put2List(flag,currentMshFace,listMshFace);
        			if(flag!=0){
        				bFlag++;
        			}
        		//}
        	}

        	if(bFlag==0 && currentMshFace->state!=10 && currentMshFace->state!=1){
        		currentMshFace->state=0;
        		mshBFacesList.push_back(currentMshFace);
        		//cout<< currentMshFace->nid << " (" << currentMshFace->ptIndxList[0] <<"," << currentMshFace->ptIndxList[1]  <<"," << currentMshFace->ptIndxList[2]  <<"," << currentMshFace->ptIndxList[3]<< ") " << currentMshFace->state <<"\n";
        	}
        	bFlag=0;
    	}
    }
}
void cOctree::put2List(int flag,cFace* currentMshFace,cFace* listMshFace)
{
    // 0 - not in
	// 1 - same face
	// 2 - f2 contains f1
	// 3 - f1 contains f2
 	if(flag==1){//same faces
        //keep only one internal face dir based on "small volLabel point to large volLabel"
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, then swap with the one in the list
        	currentMshFace->state=1;
        	listMshFace->state=1;
            currentMshFace->nbr = listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
        else {//current is larger, then keep it in the list
        	currentMshFace->state=1;
        	listMshFace->state=1;
        	listMshFace->nbr = currentMshFace->mshVolIndx; //current is large, then keep the list
        	mshIntlFacesList.push_back(listMshFace);
        }
	}else if(flag==2){// current face belongs to face in list-> then swap with the list
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, then keep anti-clock-wise order of the face
            listMshFace->state = 10;//this is big face
            currentMshFace->state=1;
            currentMshFace->nbr = listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
        else {//current is larger, then change to clock-wise order of the current face
        	listMshFace->state=10;//this is big face
        	currentMshFace->state=1;
        	currentMshFace->changeOrder();
            currentMshFace->nbr=currentMshFace->mshVolIndx;
            currentMshFace->own=listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
	}else if(flag==3){//current face contains face in list -> keep the one in list
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, change order of the face in list
        	listMshFace->state=1;
        	currentMshFace->state=10;//this is big face

        	listMshFace->changeOrder();
            listMshFace->nbr = listMshFace->mshVolIndx;
            listMshFace->own = currentMshFace->mshVolIndx;
            mshIntlFacesList.push_back(listMshFace);

            //listMshFace=currentMshFace;
        }
        else {//current is larger, then keep the order of the face in list
        	listMshFace->state=1;
        	currentMshFace->state=10;//this is big face

        	listMshFace->nbr=currentMshFace->mshVolIndx;
            mshIntlFacesList.push_back(listMshFace);

            //listMshFace=currentMshFace;
        }
	}else if(flag==0){

	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctree::setup_nodePhyName() {

	cOctNode *node=NULL;
	//vector<cTri>& geoFFacesList;
	cFace* mshFace=NULL;
    int numOfBFaces=0;

   //cout<<"BFaces in bNodesList "<<"\n";

	for(unsigned i=0; i<bNodesList.size();i++){
		node=bNodesList[i];
		vector<cTri*> geoFFacesList=node->geoFFacesList;
		for(unsigned j=0; j<node->mshFacesList.size();j++){
			mshFace=&(node->mshFacesList[j]);
			//judge mshFace -> phyName
			if(mshFace->state==0){
//				if(mshFace->nid=="node-0-3-6->face-0"){
//					cout<<"\n";
//				}
				mshFace->findPhyName(geoFFacesList);
				numOfBFaces++;

				//cout<<mshFace->nid<<" "<< mshFace->state<<"\n";
			}
		}
	}
	cout<<"Number of PhyName assigned is "<<numOfBFaces<<"\n";
}
void cOctree::setup_mshBdsList() {
	//init mshBdsList based on numOfPhynames
	mshBdsList.resize(numOfGeoPhyNames+1);
	for(int i=0; i<numOfGeoPhyNames;i++){
		mshBdsList[i].phyName=geoPhyNamesList[i];
		mshBdsList[i].phyNameIndx=i;
	}
	mshBdsList[numOfGeoPhyNames].phyName="unAssigned";
	mshBdsList[numOfGeoPhyNames].phyNameIndx=-100;

	//put mshFace into mshBdsList based on it's phyNameIndx
	cOctNode *node=NULL;
	cFace* mshFace=NULL;
	int indx;
	int numOfBFaces=0;
	for(unsigned int i=0; i<bNodesList.size();i++){
		node=bNodesList[i];
		for(unsigned int j=0; j<node->mshFacesList.size();j++){
			mshFace=&(node->mshFacesList[j]);
			//judge mshFace -> phyName
			if(mshFace->state==0){
				numOfBFaces++;
				if(mshFace->phyName!="unAssigned"){
					indx=mshFace->phyNameIndx;
					mshBdsList[indx].mshFacesList.push_back(mshFace);
				}else{
					mshBdsList[numOfGeoPhyNames].mshFacesList.push_back(mshFace);
				}

			}else{
				continue;
			}
		}
	}

	//complete other properties of the List
	mshBdsList[0].startFace=mshIntlFacesList.size();
	mshBdsList[0].nFaces=mshBdsList[0].mshFacesList.size();
	for(int i=1; i<numOfGeoPhyNames+1;i++){
		mshBdsList[i].nFaces=mshBdsList[i].mshFacesList.size();
		mshBdsList[i].startFace=mshBdsList[i-1].nFaces+mshBdsList[i-1].startFace;
	}

	cout<<"Number of BFaces in mshBdsList is "<<numOfBFaces<<"\n\n";

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// +++++ save mesh
void cOctree::saveAsGmsh(const char* _fileName)
{
	ofstream myFile;
	myFile.open(_fileName, ios::app);
	if (!myFile.is_open()){
		cout << "Open file failure" << endl;
	}

	//general mesh format
	myFile << "$MeshFormat" << endl;
	myFile << "2.2 0 8" << endl; //format_version Ascii_type data_size
	myFile << "$EndMeshFormat" << endl;

	//physical entity -> boundary name
	myFile << "$PhysicalNames" << endl;
	myFile << cGeomData::getInstance()->numOfPhyNames << endl;
	for(int i=0;i<cGeomData::getInstance()->numOfPhyNames; i++)
	{
		myFile << cGeomData::getInstance()->phyNamesList[i] <<endl;
	}
	myFile << "$EndPhysicalNames" << endl;

	//vertices list
	myFile << "$Nodes" << endl;

	myFile << "$EndNodes" << endl;

	//element list
	myFile << "$Elements" << endl;
	//boundary 2D element -> associate with geom entity and physical entity

	//3D element -> associate with geom entity and physical entity

	myFile << "$End$Elements" << endl;

	myFile.close();
}

void cOctree::saveAsOFMesh() {
	saveAsOFMeshPts();
	saveAsOFMeshFaces();
	saveAsOFMeshNeis();
	saveAsOFMeshOwns();
	saveAsOFMeshBds();
}

void cOctree::saveAsOFMeshPts() {
	string fileName="points";
	ofstream myFile;
		myFile.open(fileName, ios::ate|ios::out);
		if (!myFile.is_open()){
			cout << "Open file failure" << endl;
		}
		//general mesh format
		myFile << "FoamFile" << endl;
		myFile << "{" << endl; //format_version Ascii_type data_size
		myFile << "    version     2.0;" << endl;
		myFile << "    format      ascii;" << endl;
		myFile << "    class       vectorField;" << endl;
		myFile << "    location    \"constant/polyMesh\";" << endl;
		myFile << "    object      points;" << endl;
		myFile << "}" << endl;
		myFile << "" << endl;
		//points
		myFile << mshPtsList.size() << endl;
		myFile << "(" << endl;
		for(unsigned i=0;i<mshPtsList.size();i++){
			vector<double> pt;
			pt=mshPtsList[i];
			myFile << "(" << pt[0] << " " << pt[1] << " " << pt[2] << ")\n";
		}
		myFile << ")" << endl;
		myFile.close();
}

void cOctree::saveAsOFMeshFaces() {
	string fileName="faces";
	ofstream myFile;
		myFile.open(fileName, ios::ate|ios::out);//open file for writing; if file exists, overwrite
		if (!myFile.is_open()){
			cout << "Open file failure" << endl;
		}
		//general mesh format
		myFile << "FoamFile" << endl;
		myFile << "{" << endl; //format_version Ascii_type data_size
		myFile << "    version     2.0;" << endl;
		myFile << "    format      ascii;" << endl;
		myFile << "    class       faceList;" << endl;
		myFile << "    location    \"constant/polyMesh\";" << endl;
		myFile << "    object      faces;" << endl;
		myFile << "}" << endl;
		myFile << "" << endl;
		//points
		myFile << num_mshBFaces + num_mshIntlFaces << endl;
		myFile << "(" << endl;
		for(unsigned i=0;i<mshIntlFacesList.size();i++){
			cFace* face;
			face=mshIntlFacesList[i];
			int nPts=face->nPts ;
			ostringstream oss;
			oss<< nPts << "(";
			for(int j=0;j<nPts;j++){
				oss<< face->ptIndxList[j]<< " ";
			}
			myFile << oss.str() << ")\n";

//			cout<<setiosflags(ios::fixed)<<setprecision(2);
//			cout<< face->nid << " <" << face->ptIndxList[0]<< "("<<face->ptsList[0][0] <<","<<face->ptsList[0][1] <<","<<face->ptsList[0][2] <<"), "
//									 << face->ptIndxList[1]<< "("<<face->ptsList[1][0] <<","<<face->ptsList[1][1] <<","<<face->ptsList[1][2] <<"), "
//					                 << face->ptIndxList[2]<< "("<<face->ptsList[2][0] <<","<<face->ptsList[2][1] <<","<<face->ptsList[2][2] <<"), "
//					                 << face->ptIndxList[3]<< "("<<face->ptsList[3][0] <<","<<face->ptsList[3][1] <<","<<face->ptsList[3][2] <<")> "
//			    << face->state <<" "<< face->phyName <<"\n";

			oss.clear();
		}

		for(unsigned i=0;i<mshBdsList.size();i++){
			vector<cFace*> mshFacesList=mshBdsList[i].mshFacesList;
			for(unsigned j=0;j<mshFacesList.size();j++){
				cFace* face;
				face=mshFacesList[j];
				int nPts=face->nPts ;
				ostringstream oss;
				oss<< nPts << "(";
				for(int k=0;k<nPts;k++){
					oss<< face->ptIndxList[k]<< " ";
				}
				myFile << oss.str() << ")\n";
				//cout<< face->nid << " (" << face->ptIndxList[0] <<"," << face->ptIndxList[1]  <<"," << face->ptIndxList[2]  <<"," << face->ptIndxList[3]<< ") " << face->state <<" "<< face->phyName <<"\n";

				oss.clear();
			}
		}
		myFile << ")" << endl;
		myFile.close();
}

void cOctree::saveAsOFMeshNeis() {
	string fileName="neighbour";
	ofstream myFile;
		myFile.open(fileName, ios::ate|ios::out);
		if (!myFile.is_open()){
			cout << "Open file failure" << endl;
		}
		//general mesh format
		myFile << "FoamFile" << endl;
		myFile << "{" << endl; //format_version Ascii_type data_size
		myFile << "    version     2.0;" << endl;
		myFile << "    format      ascii;" << endl;
		myFile << "    class       labelList;" << endl;
		myFile << "    note       \"nPoints:" << mshPtsList.size()<< "\t nCells:"<< countAllNodes(root) << "\t nFaces:"<<num_mshBFaces + num_mshIntlFaces<<"\t nInternalFaces:" <<num_mshIntlFaces<<"\";\n";
		myFile << "    location    \"constant/polyMesh\";" << endl;
		myFile << "    object      neighbour;" << endl;
		myFile << "}" << endl;
		myFile << "" << endl;
		//points
		myFile << num_mshIntlFaces << endl;
		myFile << "(" << endl;
		for(int i=0;i<num_mshIntlFaces;i++){
			cFace* face;
			face=mshIntlFacesList[i];
			myFile << face->nbr << "\n";
		}
		myFile << ")" << endl;
		myFile.close();
}

void cOctree::saveAsOFMeshOwns() {
	string fileName="owner";
	ofstream myFile;
		myFile.open(fileName, ios::ate|ios::out);
		if (!myFile.is_open()){
			cout << "Open file failure" << endl;
		}
		//general mesh format
		myFile << "FoamFile" << endl;
		myFile << "{" << endl; //format_version Ascii_type data_size
		myFile << "    version     2.0;" << endl;
		myFile << "    format      ascii;" << endl;
		myFile << "    class       labelList;" << endl;
		myFile << "    note       \"nPoints:" << mshPtsList.size()<< "\t nCells:"<< countAllNodes(root) << "\t nFaces:"<<num_mshBFaces + num_mshIntlFaces<<"\t nInternalFaces:" <<num_mshIntlFaces<<"\";\n";
		myFile << "    location    \"constant/polyMesh\";" << endl;
		myFile << "    object      owner;" << endl;
		myFile << "}" << endl;
		myFile << "" << endl;
		//points
		myFile << num_mshBFaces+num_mshIntlFaces << endl;
		myFile << "(" << endl;
		for(unsigned i=0;i<mshIntlFacesList.size();i++){
			cFace* face;
			face=mshIntlFacesList[i];
			myFile << face->own << "\n";
		}

		for(unsigned i=0;i<mshBdsList.size();i++){
			vector<cFace*> mshFacesList=mshBdsList[i].mshFacesList;
			for(unsigned j=0;j<mshFacesList.size();j++){
				cFace* face;
				face=mshFacesList[j];
				myFile << face->own << "\n";
			}
		}

		myFile << ")" << endl;
		myFile.close();
}

void cOctree::saveAsOFMeshBds() {
	string fileName="boundary";
	ofstream myFile;
		myFile.open(fileName, ios::ate|ios::out);
		if (!myFile.is_open()){
			cout << "Open file failure" << endl;
		}
		//general mesh format
		myFile << "FoamFile" << endl;
		myFile << "{" << endl; //format_version Ascii_type data_size
		myFile << "    version     2.0;" << endl;
		myFile << "    format      ascii;" << endl;
		myFile << "    class       polyBoundaryMesh;" << endl;
		myFile << "    location    \"constant/polyMesh\";" << endl;
		myFile << "    object      boundary;" << endl;
		myFile << "}" << endl;
		myFile << "" << endl;
		//points
		myFile << mshBdsList.size() << endl;
		myFile << "(" << endl;

		for(unsigned i=0;i<mshBdsList.size();i++){
			string phyName=mshBdsList[i].phyName;
			string patchTpye=mshBdsList[i].patchType;
			int nFaces=mshBdsList[i].nFaces;
			int startFace=mshBdsList[i].startFace;
			ostringstream oss;
			oss<<"\t"<< phyName << "\n\t{\n";
			oss<<"\t\t"<<"type    "<<patchTpye<<";\n";
			oss<<"\t\t"<<"nFaces    "<<nFaces<<";\n";
			oss<<"\t\t"<<"startFace    "<<startFace<<";\n";
			oss<<"\t}\n";
			myFile << oss.str() << endl;
			oss.clear();
		}

		myFile << ")" << endl;
		myFile.close();
}


// +++++ save mesh

// ------------------------------------------------------

