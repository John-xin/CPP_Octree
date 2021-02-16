
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
    //++++++++++++++++++++++++++++++++++++++
    int _offsets[][3] = {{-1,-1,-1},{+1,-1,-1},{-1,+1,-1},{+1,+1,-1},
                         {-1,-1,+1},{+1,-1,+1},{-1,+1,+1},{+1,+1,+1}};

    for (int i=0; i<8; i++) {
        for (int j=0; j<3; j++) {
            branchOffsets[i][j] = _offsets[i][j];
        }
    }

    depth=0;

}
cOctree::~cOctree()
{
    cout << "Destroying the cOctree" << endl;
}

int cOctree::countAllNodes(cOctNode* node)
{
    int numOfNodes = 0;

    if (node->isLeafNode()) {
        numOfNodes++;
        return numOfNodes;
    }
    else {
        numOfNodes++;
        for (unsigned int i = 0; i < node->children.size();i++) {
            numOfNodes+=countAllNodes(node->children[i]);

        }
        return numOfNodes;

    }
}

int cOctree::countIntlNodes(cOctNode* node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<leafNodesList.size();i++){
    	if(leafNodesList[i]!=NULL && leafNodesList[i]->state==NodeState::interior){
    		numOfNodes++;
    	}
    }
    return numOfNodes;
}

int cOctree::countBNodes(cOctNode* node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<leafNodesList.size();i++){
    	if(leafNodesList[i]!=NULL && leafNodesList[i]->state==NodeState::boundary){
    		numOfNodes++;
    	}
    }
    return numOfNodes;
}

void cOctree::showOctreeNodes()
{
	num_allNodes = countAllNodes(&root);
	num_intlNodes = countIntlNodes(&root);
	num_extNodes = countExtNodes(&root);
	num_bNodes = countBNodes(&root);
	std::cout << "num of Nodes(root->all leaves) is " << num_allNodes << "\n";
	std::cout << "num of leaf Nodes is " << leafNodesList.size() << "\n";
	std::cout << "num of Interior Nodes is " << num_intlNodes << "\n";
	std::cout << "num of Exterior Nodes is " << num_extNodes << "\n";
	std::cout << "num of Boundary Nodes is " << num_bNodes << "\n\n";
}

int cOctree::countExtNodes(cOctNode* node) {
    int numOfNodes = 0;
    for(unsigned i=0; i<leafNodesList.size();i++){
    	if(leafNodesList[i]==NULL){
    		numOfNodes++;
    	}
    }
    return numOfNodes;
}


// +++++ setup root
void cOctree::setup_root() {
    vector<double> position = getPositionRoot();
	//vector<double> position = setCustomPosition(3,3,63);
    double size = getSizeRoot();
    root = cOctNode(0,"0", position, size, cGeomData::getInstance()->geoFPtsList, cGeomData::getInstance()->geoFEdgesList, cGeomData::getInstance()->geoFFacesList, NULL);
}
double cOctree::getSizeRoot() {

    // Get low and upp
    vector<double> low, upp, range;
    low = cGeomData::getInstance()->low;
    upp = cGeomData::getInstance()->upp;

    // Range is the size of the node in each coord direction
    //range[0]=len_x; range[1]=len_y; range[2]=len_z
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
vector<double> cOctree::setCustomPosition(double x,double y,double z) {
	vector<double> pos(3);
	pos[0]=x;
	pos[1]=y;
	pos[2]=z;
	return pos;
}
// +++++ setup root

void cOctree::cutNode(cOctNode* node)
{
	vector<double> newPos;
	newPos.resize(3);
	for (int i = 0; i < node->NUM_BRANCHES_OCTNODE; i++) {
		for (int j = 0; j < 3; j++) {
			newPos[j] = node->position[j] + 0.25 * node->size * branchOffsets[i][j];
		}
		string nid = node->nid + "-" + commonFunc::getInstance()->NumberToString(i);
		node->addNode(node->level + 1, nid, newPos, 0.5 * node->size,
			node->geoFPtsList,
			node->geoFEdgesList,
			node->geoFFacesList,
			node);
	}

	for (cOctNode* child : node->children) {
		setup_boundaryNode(child);
		setup_interiorNode(child);
	}

	node->geoFPtsList.resize(0);
	node->geoFEdgesList.resize(0);
	node->geoFFacesList.resize(0);
}

// +++++++++ split node
void cOctree::splitNode(cOctNode* node) {
	vector<cOctNode*> nodesList;
    // Split node into 8 branches
	cutNode(node);
	for (cOctNode* child : node->children) {
		nodesList.push_back(child);
	}
	
	// if the node level is larger than the min-octree level, then make the nbr nodes balance to 2:1
	// e.g. cut node is level 2, min level is 1. After cutting, node's children is level 3, the octree becomes 3:1 
	if (node->level - MIN_OCTREE_LEVELS >= 1) {
		vector<cOctNode*> nbrNodes = getNbr18Nodes(node);
		for (cOctNode* nbrNode : nbrNodes) {
			if (node->level - nbrNode->level >= 1) {
				cutNode(nbrNode);

				for (cOctNode* child : nbrNode->children) {
					nodesList.push_back(child);
				}

				vector<cOctNode*> nbr6Nodes = getNbr6Nodes(nbrNode);
				for (cOctNode* myNode : nbr6Nodes) {
					nodesList.push_back(myNode);
				}
			}
		}
	}
	//deal with nbr setting
	//need to optimize. there are many repeated nodes in the list
	vector<cOctNode*> nbr6Nodes = getNbr6Nodes(node);
	for (cOctNode* myNode : nbr6Nodes) {
		nodesList.push_back(myNode);
	}


	for (cOctNode* myNode : nodesList) {
		if (myNode->nid=="0-1-6") {
			cout << "";
		}
		setLeafNodeNbr(myNode);
	}

}
void cOctree::splitOctreeByFeaturePt(cOctNode* node) {

	if(node->isLeafNode()){
		splitNodeByFeaturePt(node);
	} else{
		for(int i=0; i<node->NUM_BRANCHES_OCTNODE; i++){
			splitOctreeByFeaturePt(node->children[i]);
		}
	}
}
void cOctree::splitNodeByFeaturePt(cOctNode *node)
{
	if (node->numOfGeoFPts()> MAX_OCTNODE_FEATS && node->level<MAX_OCTREE_LEVELS){
		// Split node into 8 branches
		splitNode(node);

		// Reallocate date from node to branches
		for (int i = 0; i < node->NUM_BRANCHES_OCTNODE; i++) {
			splitNodeByFeaturePt(node->children[i]);
		}
	}
}
void cOctree::splitOctreeByMinLevel(cOctNode* node) {
	//this sets global nodes level
	if(node->isLeafNode()){
		splitNodeByLevel(MIN_OCTREE_LEVELS,node);
	} else{
		for(int i=0; i<node->NUM_BRANCHES_OCTNODE; i++){
			splitOctreeByMinLevel(node->children[i]);
		}
	}
}
void cOctree::splitNodeByLevel(int level, cOctNode *node)
{
    if (node->level< level){
        // Split node into 8 branches
		splitNode(node);
        // split child node until min level
        for (int i=0; i<node->NUM_BRANCHES_OCTNODE; i++) {
            splitNodeByLevel(level, node->children[i]);
        }
    }
}
void cOctree::getOctreeDepth(cOctNode *node)
{
	if(node->isLeafNode()){
		if(node->level > depth){
			depth = node->level;
		}
	} else{
		for(int i=0; i<node->NUM_BRANCHES_OCTNODE; i++){
			getOctreeDepth(node->children[i]);
		}
	}
}
void cOctree::balanceOctree(cOctNode* node) {
	if(node->isLeafNode()){
		splitNodeByLevelDiff(1, node);
	} else{
		for(int i=0; i<node->NUM_BRANCHES_OCTNODE; i++){
			balanceOctree(node->children[i]);
		}
	}
}

void cOctree::balOctree2to1(vector<cOctNode*> nodesList)
{
	int checkCount = 0;
	queue<cOctNode*> toCheckNodes;
	for (cOctNode* node : nodesList) {
		toCheckNodes.push(node);
	}

	while (toCheckNodes.size() != 0) {
		cOctNode* myNode = toCheckNodes.front();
		if (myNode->isLeafNode()) {
			if (!is2to1Bal(myNode)) {
				splitNodeBy2to1Bal(myNode);
				vector<cOctNode*> introducedNodes = getNbr6Nodes(myNode);
				checkCount++;
				for (cOctNode* node : introducedNodes) {
					toCheckNodes.push(node);
				}
			}
		}
		toCheckNodes.pop();
	}

	cout << "Check 2to1Bal: Number of nodes made 2to1 is " << checkCount << "\n";
	update_leafNodesList(&root);
}
bool cOctree::is2to1Bal(cOctNode* node)
{
	//if any nbr level - node level > 1 then 2to1 Bal test is false
	for (vector<cOctNode*> nbr : node->nbrsList) {
		for (cOctNode* nbrNode : nbr) {
			if (nbrNode->level - node->level > 1) {
				return false;
			}
		}
	}
	return true;
}
void cOctree::splitNodeBy2to1Bal(cOctNode* node)
{
	if (node->isLeafNode()) {
		splitNode(node);
		//2. Setup the new split nodes
		
		for (cOctNode* child : node->children) {
				setup_boundaryNode(child);
				setup_interiorNode(child);
				setLeafNodeNbr(child);
		}

		//3. Repair the nbr relationship of the split node's nbrs and introduce new possible error nodes to check
		vector<cOctNode*> nbr6Nodes = getNbr6Nodes(node);
		for (cOctNode* myNode : nbr6Nodes) {
			setLeafNodeNbr(myNode);
		}
		//4. check children
		for (cOctNode* child : node->children) {
			if (!is2to1Bal(node->children[0])) {
				splitNodeBy2to1Bal(child);
			}
		}
	}


}
void cOctree::splitNodeByLevelDiff(int lvlDiff, cOctNode *node){
	if (depth- node->level > lvlDiff){
		splitNode(node);
	    // split child node until level diff less than 1
	    
		//for (int i=0; i<node->NUM_BRANCHES_OCTNODE; i++) {
	    //   if (depth-node->children[i]->level > lvlDiff) {
	    //       splitNodeByLevelDiff(lvlDiff, node->children[i]);
	    //    }
	    //}
	}
}
void cOctree::splitNodeById(string node_id) {
	cOctNode* node=getNodeFromId(node_id);
	splitNode(node);
}
void cOctree::splitNodeByPhyName(string phyName, int level, cOctNode* node) {
	if (node->isLeafNode() && node->level < level) {
		if (node->geoFFacesList.size() != 0) {
			for (unsigned i = 0; i < node->geoFFacesList.size(); i++) {
				if (node->geoFFacesList[i]->phyName == phyName) {
					break;
				}
			}
			splitNode(node);
			for (unsigned i = 0; i < node->children.size(); i++) {
				splitNodeByPhyName(phyName, level, node->children[i]);
			}
		}
	}
	else {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			splitNodeByPhyName(phyName, level, node->children[i]);
		}
	}
}
// +++++++++ split node

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctree::setup_boundaryNode(cOctNode* node) {
	if (node->isLeafNode()) {
		if (node->geoFFacesList.size() != 0) {
			node->state = NodeState::boundary; //boundary node
		}
		else {
			node->state = NodeState::nonBoundary;; //non-boundary node
		}
	}
	else {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			setup_boundaryNode(node->children[i]);
		}
	}
}
void cOctree::setup_interiorNode(cOctNode* node) {
//	cOctNode *node;
//	node=nonBNodesList[5];
//	int flag;
//	flag=isInteriorNode(node);
//	node->isInteriorNode=flag;
//	int numOfNodes=0;
//
//	if (flag==1){
//		setup_nbrNodesState(node);
//		for(unsigned i=0; i<nonBNodesList.size();i++){
//			if(nonBNodesList[i]->isInteriorNode!=1){
//				nonBNodesList[i]->isInteriorNode=0;
//				numOfNodes++;
//			}
//		}
//	}else if (flag==0){
//		setup_nbrNodesState(node);
//		for(unsigned i=0; i<nonBNodesList.size();i++){
//			if(nonBNodesList[i]->isInteriorNode!=0){
//				nonBNodesList[i]->isInteriorNode=1;
//				numOfNodes++;
//			}
//		}
//	}

//	if(node->isInteriorNode==1){
//		cout<<"number of Interior nodes is "<<numOfNodes<<"\n";
//	}else{
//		cout<<"number of Exterior nodes is "<<numOfNodes<<"\n\n";
//	}

//----brutal force method---------------------------------------------------
	if (node->isLeafNode()) {
		if (node->state== NodeState::nonBoundary) {
			if (isInteriorNode(node)) {
				node->state = NodeState::interior;
			}else {
				node->state = NodeState::exterior;
			}
		}
	}
	else {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			setup_interiorNode(node->children[i]);
		}
	}
}
int cOctree::isInteriorNode(cOctNode* node) {
	//assume input STL is closed
    //retrieve a node from nonBNodeslist
	//make a ray -> node position to defined point
	//find num of intersections
	//even -> ext node
	//odd -> intl node

	vector<Intersection> intersections;
	cLine ray(node->position, cGeomData::getInstance()->ptInGeom,0); //depends on defindBody() to initialize ptInGeom;
	intersections=findRayIntersect(ray);

	int numOfInts=0;
	for(unsigned i=0;i<intersections.size();i++){
		if(intersections[i].s>0) {numOfInts++;}
	}

	if((numOfInts%2)==0){
		return 0;//0 - exterior node
		//node->isInteriorNode=0;
	}else{
		return 1; //1 - interior node
		//node->isInteriorNode=1; //1 - internal node
		//for(unsigned j=0; j<node->mshFacesList.size();j++){
		//	node->mshFacesList[j].isBoundaryFace=0;
		//}
	}
}

void cOctree::setup_nbrNodesState(cOctNode* node) {
	cOctNode* nbr;
	for(int i=0; i<6; i++){
		if(node->nbrsList[i].size()!=0){
			for(unsigned j=0;j<node->nbrsList[i].size();j++){
				nbr=node->nbrsList[i][j];
				if(nbr->state== NodeState::nonBoundary){
					nbr->state=node->state;
					setup_nbrNodesState(nbr);
				}
			}
		}
	}
}

vector<cOctNode*> cOctree::getNbr6Nodes(cOctNode* node)
{
	//get nbrs by moving center pt to a circle edge with a radius of 0.75 * size
	//prerequisite:
	//1. octree must be 2to1 balance
	//2. node must be a cube

	vector<vector <double> > vect(6);
	for (int i = 0; i < 6; i++) {
		vect[i].resize(3);
	}
	double r = 0.75 * node->size;
	vect[0] = { 0, -r, 0 };//nbr s
	vect[1] = { r, 0, 0 };//nbr e
	vect[2] = { 0,  r, 0 };//nbr n
	vect[3] = { -r, 0, 0 };//nbr w
	vect[4] = { 0, 0,  r };//nbr t
	vect[5] = { 0, 0, -r };//nbr b

	vector<cOctNode*> nbr6Nodes;
	for (int j = 0; j < 6; j++) {
		vector<double> pt = node->position;
		pt = commonFunc::getInstance()->vectAdd(pt, vect[j]);

		//return list of nodes for each nbr, e.g. low-level node nbr has 4  high-level nodes
		//0-s 1-e 2-n 3-w 4-b 5-t
		vector<cOctNode*> nbrNodes = getLeafNodeByPt(pt, &root);
		for (cOctNode* nbrNode : nbrNodes) {
			bool flag = true;
			for (cOctNode* aNode : nbr6Nodes) {
				if (nbrNode->nid == aNode->nid) {
					flag = false;
				}
			}
			if (flag) {
				nbr6Nodes.push_back(nbrNode);
			}
		}
	}
	return nbr6Nodes;
}

vector<cOctNode*> cOctree::getNbr18Nodes(cOctNode* node)
{
	//get nbrs by moving center pt to a circle edge with a radius of 0.75 * size
	//prerequisite:
	//1. octree must be 2to1 balance
	//2. node must be a cube

	vector<vector <double> > vect(18);
	for (int i = 0; i < 18; i++) {
		vect[i].resize(3);
	}
	double r = 0.75*node->size;
	double pi = 3.1415;
	vect[0] = { 0, -r, 0 };//nbr s
	vect[1] = {  r, 0, 0 };//nbr e
	vect[2] = { 0,  r, 0 };//nbr n
	vect[3] = { -r, 0, 0 };//nbr w
	vect[4] = { 0, 0,  r };//nbr t
	vect[5] = { 0, 0, -r };//nbr b

	vect[6] = { r * sin(pi * 0.25), -r * cos(pi * 0.25), 0 };//nbr z-se
	vect[7] = { r * sin(pi * 0.75), -r * cos(pi * 0.75), 0 };//nbr z-ne
	vect[8] = { r * sin(pi * 1.25), -r * cos(pi * 1.25), 0 };//nbr z-nw
	vect[9] = { r * sin(pi * 1.75), -r * cos(pi * 1.75), 0 };//nbr z-sw
	vect[10] = { 0, r * sin(pi * 0.25), -r * cos(pi * 0.25) };//nbr x-bn
	vect[11] = { 0, r * sin(pi * 0.75), -r * cos(pi * 0.75) };//nbr x-tn
	vect[12] = { 0, r * sin(pi * 1.25), -r * cos(pi * 1.25) };//nbr x-ts
	vect[13] = { 0, r * sin(pi * 1.75), -r * cos(pi * 1.75) };//nbr x-bs
	vect[14] = { r * sin(pi * 0.25), 0, -r * cos(pi * 0.25) };//nbr y-be
	vect[15] = { r * sin(pi * 0.75), 0, -r * cos(pi * 0.75) };//nbr y-te
	vect[16] = { r * sin(pi * 1.25), 0, -r * cos(pi * 1.25) };//nbr y-tw
	vect[17] = { r * sin(pi * 1.75), 0, -r * cos(pi * 1.75) };//nbr y-bw

	vector<cOctNode*> nbr18Nodes;
	for (int j = 0; j < 18; j++) {
		vector<double> pt = node->position;
		pt = commonFunc::getInstance()->vectAdd(pt, vect[j]);

		//return list of nodes for each nbr, e.g. low-level node nbr has 4  high-level nodes
		//0-s 1-e 2-n 3-w 4-b 5-t
		vector<cOctNode*> nbrNodes = getLeafNodeByPt(pt, &root);
		for (cOctNode* nbrNode : nbrNodes) {
			bool flag = true;
			for (cOctNode* aNode: nbr18Nodes) {
				if (nbrNode->nid==aNode->nid) {
					flag = false;
				}
			}
			if (flag) {
				nbr18Nodes.push_back(nbrNode);
			}	
		}
	}
	return nbr18Nodes;
}

void cOctree::setup_leafNodesList(cOctNode* node) {

	if (node->isLeafNode()) {
		leafNodesList.push_back(node);
	}
	else {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			setup_leafNodesList(node->children[i]);
		}
	}
}

void cOctree::update_leafNodesList(cOctNode* node)
{
	leafNodesList.resize(0);
	setup_leafNodesList(node);
}

void cOctree::setup_leafNodesNbr() {
	//assume node level in the tree is less than one
	for (unsigned i = 0; i < leafNodesList.size(); i++) {
		cOctNode* node = leafNodesList[i];
		setLeafNodeNbr(node);
	}
}
void cOctree::setLeafNodeNbr(cOctNode* node) {
	node->nbrsList.resize(6);
	for (unsigned j = 0; j < 6; j++) {
		node->nbrsList[j].resize(0);
	}

	vector<vector <double> > vect(6);
	for (int i = 0; i < 6; i++) {
		vect[i].resize(3);
	}
	double r = 0.75 * node->size;
	vect[0] = { 0, -r, 0 };//nbr s
	vect[1] = { r, 0, 0 };//nbr e
	vect[2] = { 0,  r, 0 };//nbr n
	vect[3] = { -r, 0, 0 };//nbr w
	vect[4] = { 0, 0,  r };//nbr t
	vect[5] = { 0, 0, -r };//nbr b

	for (int j = 0; j < 6; j++) {
		vector<double> pt = node->position;
		pt = commonFunc::getInstance()->vectAdd(pt, vect[j]);

		//return list of nodes for each nbr, e.g. low-level node nbr has 4  high-level nodes
		//0-s 1-e 2-n 3-w 4-t 5-b
		node->nbrsList[j] = getLeafNodeByPt(pt, &root);
	}
}

void cOctree::splitExtNodeNbr(cOctNode* node) {
	if (node->isLeafNode()) {
		if (node->state == NodeState::exterior) {
			setLeafNodeNbr(node);
			for (int i = 0; i < 6; i++) {
				for (unsigned j = 0; j < node->nbrsList[i].size(); j++) {
					if (node->level > node->nbrsList[i][j]->level) {
						splitNode(node);
					}
				}
			}
		}
	}
	else {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			splitExtNodeNbr(node->children[i]);
		}
	}
}
void cOctree::delExtNodes()
{
	int indxArr[6]={2,3,0,1,5,4};
	for(unsigned i=0;i<leafNodesList.size();i++){
		cOctNode* &node=leafNodesList[i];

    	if(node->state==NodeState::exterior){

    		for(unsigned j=0; j<6;j++){
    			if(node->nbrsList[j].size()!=0){
					for(unsigned k=0;k<node->nbrsList[j].size();k++){
						vector<cOctNode*> &nbr=node->nbrsList[j][k]->nbrsList[indxArr[j]];
						nbr.resize(0);

						//this method becomes invalid if a small node is ext, then the big nbr node
						//would set the big nbr node's 4 small nbr node all 0, as a result big node cannot find small nodes but
						//small node can find big node

						//to avoid this, big node should be split to small node level

//						if(nbr.size()==1){
//							nbr.resize(0);
//						}else if(nbr.size()!=1 && nbr.size()!=0){
//							for(unsigned l=0;l<nbr.size();l++){
//								if(nbr[l]->nid==node->nid) nbr.erase(nbr.begin()+l);
//							}
//						}

					}
				}
			}

    		string nid=node->nid;
    		cout<<nid.c_str()<<"\n";
    		node=NULL;
    	}
	}
}
void cOctree::delExtNodes2(cOctNode* &node)
{
	if (node->isLeafNode()) {
		if(node->state==NodeState::exterior){
			node=NULL;
		}
	}else{
	  for (unsigned int i = 0; i < node->children.size();i++) {
		  delExtNodes2(node->children[i]);
	  }
	}
}

void cOctree::update_allNodesList(cOctNode* node)
{
	allNodesList.push_back(node);
	if (node->children.size() != 0) {
		for (unsigned int i = 0; i < node->children.size(); i++) {
			update_allNodesList(node->children[i]);
		}
	}
}

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

vector<cOctNode*> cOctree::getLeafNodeByPt(vector<double> pt, cOctNode* node) {
	tmpNodesList.resize(0);
	addLeafNodeByPt(pt, node);
	return tmpNodesList;
}
void cOctree::addLeafNodeByPt(vector<double> pt, cOctNode* node) {
	if(node->isPtInNode(pt)){
	    if (!node->isLeafNode()) {
	        for (unsigned int i = 0; i < node->children.size();i++) {
	        	addLeafNodeByPt(pt, node->children[i]);
	        }
	    }
	    else {
	    	tmpNodesList.push_back(node);
	    }
	}
}
cOctNode* cOctree::getNodeFromId(string nodeId)
{
    return findBranchById(nodeId,&root);
}

cOctNode* cOctree::findBranchById(string nodeId, cOctNode *node)
{
    if (nodeId.compare(node->nid)==0) {
        return node;
    } else {
        for (unsigned int i=0; i<node->children.size(); i++) {
            cOctNode *branch = findBranchById(nodeId, node->children[i]);
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
        if (cGeomData::getInstance()->geoFFacesList[polyLabel]->rayPlaneIntersectPoint(ray,ip,s)) {
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
    getPolysToCheck(&root,ray,intTestPolys);
    return intTestPolys;
}
void cOctree::getPolysToCheck(cOctNode *node, cLine &ray, set<int> &intTestPolys)
{
    // Utility function for getListPolysToCheck. Finds all OctNodes hit by a given ray
    // and returns a list of the objects contained within
    if (node->sphereRayIntersect(ray)) {
        if (node->boxRayIntersect(ray)) {
            if (node->isLeafNode()) {
                for (int i=0; i<node->numOfGeoFFaces(); i++) {
                    intTestPolys.insert(node->geoFFacesList[i]->indx); }
            } else {
                for (int i=0; i<node->NUM_BRANCHES_OCTNODE; i++) {
                    getPolysToCheck(node->children[i],ray,intTestPolys);
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

void cOctree::addGapNodes()
{
    //find boundary face of nodes -> single face

	//find normal of vertices of boundary face

	//find intersection point between normal and tri

	//add new node -> 4 vertices and 4 intersection points
}



// ------------------------------------------------------

