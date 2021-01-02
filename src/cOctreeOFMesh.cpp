#include "cOctreeOFMesh.h"
#include "cOctree.h"

cOctreeOFMesh::cOctreeOFMesh()
{
	octree = NULL;
}

cOctreeOFMesh::~cOctreeOFMesh()
{
}

cOctreeOFMesh* cOctreeOFMesh::getInstance()
{
	static cOctreeOFMesh theInstance;
	//	if (theCommonFunc==NULL)
	//    {
	//        theCommonFunc = new commonFunc();
	//    }
	return &theInstance;
}

void cOctreeOFMesh::setup_octree(cOctree* _octree)
{
	octree = _octree;
}

void cOctreeOFMesh::setup_mshNodesList()
{
	for (unsigned int i = 0; i < octree->leafNodesList.size(); i++) {
		if (octree->leafNodesList[i]->state != NodeState::exterior) {
			mshNodesList.push_back(octree->leafNodesList[i]);
		}
	}
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeOFMesh::setup_mshPtList()
{
	for (unsigned i = 0; i < octree->leafNodesList.size(); i++) {
		if (octree->leafNodesList[i] == NULL) {
			cout << "setup_mshPtList failed at " << octree->leafNodesList[i]->nid << "\n";
		}
		else {
			addMshPtsOfNode(octree->leafNodesList[i]);
		}
	}
}
void cOctreeOFMesh::addMshPtsOfNode(cOctNode* node) //init node.elePtLabels
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
	for (unsigned int i = 0; i < node->mshPts3DList.size(); i++) {
		pt = node->mshPts3DList[i];
		if (node->mshPtsRepeatedList[i] != true) {
			nodesList = octree->getLeafNodeByPt(pt, &octree->root);
			if (nodesList.size() == 0) {
				cout << "this point (" << pt[0] << "," << pt[1] << "," << pt[2] << "," << ") is not in the octree\n";

			}
			else if (nodesList.size() == 1) {
				// this is a non-repeated pt
				mshPtsList.push_back(pt);
				node->mshPtsIndxList[i] = mshPtIndxCount;
				mshPtIndxCount++;

			}
			else {
				mshPtsList.push_back(pt);
				node->mshPtsIndxList[i] = mshPtIndxCount;
				for (unsigned j = 0; j < nodesList.size(); j++) {
					if (nodesList[j]->nid != node->nid) {
						for (unsigned k = 0; k < nodesList[j]->mshPts3DList.size(); k++) {
							if (isPtSame(pt, nodesList[j]->mshPts3DList[k])) {
								nodesList[j]->mshPtsRepeatedList[k] = true;//old pt
								nodesList[j]->mshPtsIndxList[k] = mshPtIndxCount;
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
bool cOctreeOFMesh::isPtSame(vector<double>& pt1, vector<double>& pt2)
{
	double tol = 0.000001;
	if ((abs(pt1[0] - pt2[0])) > tol) return false;
	if ((abs(pt1[1] - pt2[1])) > tol) return false;
	if ((abs(pt1[2] - pt2[2])) > tol) return false;
	return true;
}
void cOctreeOFMesh::setup_nodeMshFaces() {
	cOctNode* node;
	cFace face;
	int count = 0;
	for (unsigned int i = 0; i < mshNodesList.size(); i++) {
		node = mshNodesList[i];
		node->mshVolIndx = count;
		count++;
		node->calMshFaces();
		for (size_t j = 0; j < node->mshFacesList.size(); j++) {
			face = node->mshFacesList[j];
			face.own = node->mshVolIndx;
		}
	}
}
void cOctreeOFMesh::update_nodeMshFaces()
{
	cOctNode* node;
	cFace* face;

	for (unsigned i = 0; i < mshNodesList.size(); i++) {
		node = mshNodesList[i];
		vector<int> surrPts = getSurrPts(node);

		for (size_t j = 0; j < node->mshFacesList.size(); j++) {
			face = &node->mshFacesList[j];
			vector<int> inFacePts = face->getInFacePts(surrPts);
			if (inFacePts.size() != 9 && inFacePts.size() != 4) {
				vector<int> inOrderPts = face->getInOrderPts(inFacePts);
				face->nPts = inOrderPts.size();
				face->ptIndxList = inOrderPts;
				face->updatePtsList();
			}
		}
	}
}
vector<int> cOctreeOFMesh::getSurrPts(cOctNode* node)
{
	vector<cOctNode*> surrNodes;
	vector<cOctNode*> tmpNodes;
	vector<cOctNode*> nbr;
	vector<int> surrPts;
	int ptIndx;

	//get surr nodes -> node's nbr and nbr's nbr -> 18+6 nos
	//check if surr node is split
	//if split -> addNode's pts

	surrNodes = getSurrNodes(node);
	for (unsigned i = 0; i < 6; i++) {
		nbr = node->nbrsList[i];
		for (unsigned j = 0; j < nbr.size(); j++) {
			tmpNodes=getSurrNodes(nbr[j]); //warning -------contains repeated nodes
			for (cOctNode* node : tmpNodes) {
				surrNodes.push_back(node);
			}
		}
	}

	for (size_t i = 0; i < surrNodes.size(); i++){
		if (surrNodes[i]->level > node->level) {
			for (size_t j = 0; j < surrNodes[i]->mshPtsIndxList.size();j++) {
				ptIndx = surrNodes[i]->mshPtsIndxList[j];
				surrPts.push_back(ptIndx);
			}
		}
	}

	for (size_t j = 0; j < node->mshPtsIndxList.size(); j++) {
		ptIndx = node->mshPtsIndxList[j];
		surrPts.push_back(ptIndx);
	}

    //unique surrPts
	sort(surrPts.begin(), surrPts.end());
	surrPts.erase(unique(surrPts.begin(), surrPts.end()), surrPts.end());
	return surrPts;
}

vector<cOctNode*> cOctreeOFMesh::getSurrNodes(cOctNode* node)
{
	vector<cOctNode*> surrNodes;
	vector<cOctNode*> nbr;
	for (unsigned i = 0; i < 6; i++) {
		nbr = node->nbrsList[i];
		for (unsigned j = 0; j < nbr.size(); j++) {
		    surrNodes.push_back(nbr[j]);
		}
	}
	return surrNodes;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeOFMesh::setup_nodeMshFaceState()
{
	//boundary node -> set boundary face
	//exterior node -> set boundary face
	//boundary node next to exterior node -> not set boundary face

	cFace* currentMshFace;
	cFace* nbrMshFace;
	vector<cOctNode*> nbr;
	string curFaceName;
	string nbrMshFaceName;
	int indxArr[6] = { 2, 3, 0, 1, 5, 4 };

	for (unsigned int i = 0; i < mshNodesList.size(); i++) {
		cOctNode*& node = mshNodesList[i]; //this is wrong when null is assign to reference
		for (unsigned j = 0; j < node->mshFacesList.size(); j++) {
			nbr = node->nbrsList[j];
			currentMshFace = &node->mshFacesList[j];
			if (nbr.size() == 0) {
				node->mshFacesList[j].state = FaceState::boundary;
				node->mshFacesList[j].exportState = true;
				mshBFacesList.push_back(&node->mshFacesList[j]);
			}
			else {
				for (size_t k = 0; k < nbr.size(); k++) {
					if (nbr[k]->state == NodeState::exterior) {
						node->mshFacesList[j].state = FaceState::boundary;
						node->mshFacesList[j].exportState = true;
						mshBFacesList.push_back(&node->mshFacesList[j]);

					}
					else {
						nbrMshFace = (&nbr[k]->mshFacesList[indxArr[j]]);
						int flag = currentMshFace->findFaceRelationship(currentMshFace, nbrMshFace);
						saveInterlMshFaces(flag, currentMshFace, nbrMshFace);
					}
				}
			}
		}
	}

	// ================================================================
//    int num_IntlFaces_inMshBFacesList=0;
//    for (unsigned int i = 0; i < mshBFacesList.size(); i++) {
//    	if(mshBFacesList[i]->isBoundaryFace==0){
//    		num_IntlFaces_inMshBFacesList++;
//    	}else{
//    		num_mshBFaces++;
//    	}
//    }
//    cout<<"Number of Internal Faces in mshBFacesList is "<<num_IntlFaces_inMshBFacesList<<"\n";
//    cout<<"Number of Boundary Faces in mshBFacesList is "<<num_mshBFaces<<"\n\n";
}
void cOctreeOFMesh::setup_intlMshFace()
{
	cFace* currentMshFace;
	cFace* nbrMshFace;
	vector<cOctNode*> nbr;
	string curFaceName;
	string nbrMshFaceName;
	int indxArr[6] = { 2, 3, 0, 1, 5, 4 };

	for (unsigned int i = 0; i < mshNodesList.size(); i++) {
		cOctNode*& node = mshNodesList[i];
		if (node != NULL) {

			for (unsigned j = 0; j < node->mshFacesList.size(); j++) {
				currentMshFace = &node->mshFacesList[j];
				curFaceName = currentMshFace->nid;
				if (currentMshFace->state == FaceState::unassigned) {
					nbr = node->nbrsList[j];
					if (nbr.size() != 0) {

						for (unsigned k = 0; k < nbr.size(); k++) {
							nbrMshFace = (&nbr[k]->mshFacesList[indxArr[j]]);
							nbrMshFaceName = nbrMshFace->nid;
							if (nbrMshFaceName == "node-0-0-3-3-3-3->face-4") {
								cout << "";
							}
							int flag = currentMshFace->findFaceRelationship(currentMshFace, nbrMshFace);
							put2List(flag, currentMshFace, nbrMshFace);
						}

					}
				}
			}

		}
	}

	//////////////////////////////////////////////////////////////////

	num_mshBFaces = 0;
	num_mshIntlFaces = mshIntlFacesList.size();

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

	int num_IntlFaces_inMshBFacesList = 0;
	for (unsigned int i = 0; i < mshBFacesList.size(); i++) {
		if (mshBFacesList[i]->state == FaceState::internal) {
			num_IntlFaces_inMshBFacesList++;
			cout << mshBFacesList[i]->nid.c_str() << "\n";
		}
		else {
			num_mshBFaces++;
		}
	}
	cout << "Number of Internal Faces in mshBFacesList is " << num_IntlFaces_inMshBFacesList << "\n";
	cout << "Number of Boundary Faces in mshBFacesList is " << num_mshBFaces << "\n\n";

}
void cOctreeOFMesh::findMshIntlFaces()
{
	//add face to internalEleFacesList
   //add face to boundaryEleFacesList - phycial name 1
   //add face to boundaryEleFacesList - phycial name 2
   //add face to boundaryEleFacesList - phycial name 3

	int flag = -1;
	int bFlag = 0;
	cFace* currentMshFace;
	cFace* listMshFace;

	string curFaceName;
	string listFaceName;

	int numOfFaces;
	numOfFaces = mshAllFacesList.size();
	int k = 0;

	for (int i = 0; i < numOfFaces; i++) {
		k++;
		currentMshFace = mshAllFacesList[i];
		curFaceName = currentMshFace->nid;
		if (currentMshFace->state == FaceState::internal) {
			for (int j = k; j < numOfFaces; j++) {
				listMshFace = mshAllFacesList[j];
				listFaceName = listMshFace->nid;
				//if(listMshFace->isBoundaryNode!=1){
				flag = currentMshFace->findFaceRelationship(currentMshFace, listMshFace);
				put2List(flag, currentMshFace, listMshFace);
				if (flag != 0) {
					bFlag++;
				}
				//}
			}

			/*
			if(bFlag==0 && currentMshFace->isBoundaryFace!=10 && currentMshFace->isBoundaryFace!=1){
				currentMshFace->isBoundaryFace=0;
				mshBFacesList.push_back(currentMshFace);
				//cout<< currentMshFace->nid << " (" << currentMshFace->ptIndxList[0] <<"," << currentMshFace->ptIndxList[1]  <<"," << currentMshFace->ptIndxList[2]  <<"," << currentMshFace->ptIndxList[3]<< ") " << currentMshFace->isBoundaryNode <<"\n";
			}
			*/
			bFlag = 0;
		}
	}
}
void cOctreeOFMesh::put2List(int flag, cFace* currentMshFace, cFace* listMshFace)
{
	// 0 - not in
	// 1 - same face
	// 2 - f2 contains f1
	// 3 - f1 contains f2
	if (currentMshFace->exportState == false && listMshFace->exportState == false) {
		if (flag == 1) {//same faces
			//keep only one internal face dir based on "small volLabel point to large volLabel"
			if (currentMshFace->own < listMshFace->own) { //current is smaller, then swap with the one in the list
				currentMshFace->state = FaceState::internal;
				listMshFace->state = FaceState::internal;
				currentMshFace->nbr = listMshFace->own;
				currentMshFace->exportState = true;
				mshIntlFacesList.push_back(currentMshFace);
			}
			else {//current is larger, then keep it in the list
				currentMshFace->state = FaceState::internal;
				listMshFace->state = FaceState::internal;
				listMshFace->nbr = currentMshFace->own; //current is large, then keep the list
				listMshFace->exportState = true;
				mshIntlFacesList.push_back(listMshFace);
			}
		}
		else if (flag == 2) {// current face belongs to face in list-> then swap with the list
			if (currentMshFace->own < listMshFace->own) { //current is smaller, then keep anti-clock-wise order of the face
				listMshFace->state = FaceState::internal;//this is big face
				currentMshFace->state = FaceState::internal;
				currentMshFace->nbr = listMshFace->own;
				currentMshFace->exportState = true;
				mshIntlFacesList.push_back(currentMshFace);
			}
			else {//current is larger, then change to clock-wise order of the current face
				listMshFace->state = FaceState::internal;//this is big face
				currentMshFace->state = FaceState::internal;
				//currentMshFace->changeOrder();
				currentMshFace->nbr = listMshFace->own;
				currentMshFace->exportState = true;
				////currentMshFace->nbr = currentMshFace->own;
				////currentMshFace->own = listMshFace->own;
				mshIntlFacesList.push_back(currentMshFace);
			}
		}
		else if (flag == 3) {//current face contains face in list -> keep the one in list
			if (currentMshFace->own < listMshFace->own) { //current is smaller, change order of the face in list
				listMshFace->state = FaceState::internal;
				currentMshFace->state = FaceState::internal;//this is big face
				//listMshFace->changeOrder();
				listMshFace->nbr = currentMshFace->own;
				listMshFace->exportState = true;
				////listMshFace->nbr = listMshFace->own;
				////listMshFace->own = currentMshFace->own;
				mshIntlFacesList.push_back(listMshFace);

				//listMshFace=currentMshFace;
			}
			else {//current is larger, then keep the order of the face in list
				listMshFace->state = FaceState::internal;
				currentMshFace->state = FaceState::internal;//this is big face

				listMshFace->nbr = currentMshFace->own;
				listMshFace->exportState = true;
				mshIntlFacesList.push_back(listMshFace);

				//listMshFace=currentMshFace;
			}
		}
		else if (flag == 0) {

		}
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeOFMesh::saveInterlMshFaces(int flag, cFace* currentMshFace, cFace* listMshFace)
{
	// 0 - not in
	// 1 - same face
	// 2 - f2 contains f1
	// 3 - f1 contains f2
	if (currentMshFace->exportState == false && listMshFace->exportState == false) {
		if (flag == 1) {//same faces
			//keep only one internal face dir based on "small volLabel point to large volLabel"
			currentMshFace->state = FaceState::internal;
			listMshFace->state = FaceState::internal;
				
			currentMshFace->nbr = listMshFace->own;
			listMshFace->nbr = currentMshFace->own; 

			if (currentMshFace->nPts >= listMshFace->nPts) {
				currentMshFace->exportState = true;
				mshIntlFacesList.push_back(currentMshFace);
			}
			else {
				listMshFace->exportState = true;
				mshIntlFacesList.push_back(listMshFace);
			}
		}
		else if (flag == 2) {// current face belongs to face in list-> then swap with the list
				listMshFace->state = FaceState::internal;//this is big face
				currentMshFace->state = FaceState::internal;
				currentMshFace->nbr = listMshFace->own;
				currentMshFace->exportState = true;
				mshIntlFacesList.push_back(currentMshFace);
		}
		else if (flag == 3) {//current face contains face in list -> keep the one in list
				listMshFace->state = FaceState::internal;
				currentMshFace->state = FaceState::internal;//this is big face
				listMshFace->nbr = currentMshFace->own;
				listMshFace->exportState = true;
				mshIntlFacesList.push_back(listMshFace);
		}
		else if (flag == 0) {

		}
	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeOFMesh::setup_nodePhyName() {

	cOctNode* node = NULL;
	//vector<cTri>& geoFFacesList;
	cFace* mshFace = NULL;
	int numOfBFaces = 0;

	//cout<<"BFaces in bNodesList "<<"\n";

	for (unsigned i = 0; i < mshNodesList.size(); i++) {
		node = mshNodesList[i];
		if (node != NULL) {
			vector<cFeatureFace*> geoFFacesList = node->geoFFacesList;
			for (unsigned j = 0; j < node->mshFacesList.size(); j++) {
				mshFace = &(node->mshFacesList[j]);
				//judge mshFace -> phyName
				if (mshFace->state == FaceState::boundary) {
					//				if(mshFace->nid=="node-0-3-6->face-0"){
					//					cout<<"\n";
					//				}
					mshFace->findPhyName(geoFFacesList);
					numOfBFaces++;
					//cout<<mshFace->nid<<" "<< mshFace->state<<"\n";
				}
			}
		}
	}
	cout << "Number of PhyName assigned is " << numOfBFaces << "\n";
}
void cOctreeOFMesh::setup_mshBdsList() {
	//init mshBdsList based on numOfPhynames
	cGeomData* geom=cGeomData::getInstance();
	mshBdsList.resize(geom->numOfPhyNames + 1);
	for (int i = 0; i < geom->numOfPhyNames; i++) {
		mshBdsList[i].phyName = geom->phyNamesList[i];
		mshBdsList[i].phyNameIndx = i;
	}
	mshBdsList[geom->numOfPhyNames].phyName = "unAssigned";
	mshBdsList[geom->numOfPhyNames].phyNameIndx = -100;

	//put mshFace into mshBdsList based on it's phyNameIndx
	cOctNode* node = NULL;
	cFace* mshFace = NULL;
	int indx;
	int numOfBFaces = 0;
	for (unsigned int i = 0; i < mshNodesList.size(); i++) {
		node = mshNodesList[i];
		if (node != NULL) {
			for (unsigned int j = 0; j < node->mshFacesList.size(); j++) {
				mshFace = &(node->mshFacesList[j]);
				//judge mshFace -> phyName
				if (mshFace->state == FaceState::boundary) {
					numOfBFaces++;
					if (mshFace->phyName != "unAssigned") {
						indx = mshFace->phyNameIndx;
						mshBdsList[indx].mshFacesList.push_back(mshFace);
					}
					else {
						mshBdsList[geom->numOfPhyNames].mshFacesList.push_back(mshFace);
					}
				}
				else {
					continue;
				}
			}
		}
	}

	//complete other properties of the List
	mshBdsList[0].startFace = mshIntlFacesList.size();
	mshBdsList[0].nFaces = mshBdsList[0].mshFacesList.size();
	for (int i = 1; i < geom->numOfPhyNames + 1; i++) {
		mshBdsList[i].nFaces = mshBdsList[i].mshFacesList.size();
		mshBdsList[i].startFace = mshBdsList[i - 1].nFaces + mshBdsList[i - 1].startFace;
	}

	cout << "Number of BFaces in mshBdsList is " << numOfBFaces << "\n\n";

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeOFMesh::saveAsOFMesh() {
	saveAsOFMeshPts();
	saveAsOFMeshFaces();
	saveAsOFMeshNeis();
	saveAsOFMeshOwns();
	saveAsOFMeshBds();
}

void cOctreeOFMesh::saveAsOFMeshPts() {
	string fileName = "./output/constant/polyMesh/points";
	ofstream myFile;
	myFile.open(fileName, ios::ate | ios::out);
	if (!myFile.is_open()) {
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
	for (unsigned i = 0; i < mshPtsList.size(); i++) {
		vector<double> pt;
		pt = mshPtsList[i];
		myFile << "( " << pt[0] << " " << pt[1] << " " << pt[2] << " )\n";
	}
	myFile << ")" << endl;
	myFile.close();
}

void cOctreeOFMesh::saveAsOFMeshFaces() {
	string fileName = "./output/constant/polyMesh/faces";
	int num_mshIntlFaces = mshIntlFacesList.size();
	int num_mshBFaces = mshBFacesList.size();
	ofstream myFile;
	myFile.open(fileName, ios::ate | ios::out);//open file for writing; if file exists, overwrite
	if (!myFile.is_open()) {
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
	for (unsigned i = 0; i < mshIntlFacesList.size(); i++) {
		cFace* face;
		face = mshIntlFacesList[i];
		int nPts = face->nPts;
		ostringstream oss;
		oss << nPts << " ( ";
		for (int j = 0; j < nPts; j++) {
			oss << face->ptIndxList[j] << " ";
		}
		myFile << oss.str() << ")\n";

		//			cout<<setiosflags(ios::fixed)<<setprecision(2);
		//			cout<< face->nid << " <" << face->ptIndxList[0]<< "("<<face->ptsList[0][0] <<","<<face->ptsList[0][1] <<","<<face->ptsList[0][2] <<"), "
		//									 << face->ptIndxList[1]<< "("<<face->ptsList[1][0] <<","<<face->ptsList[1][1] <<","<<face->ptsList[1][2] <<"), "
		//					                 << face->ptIndxList[2]<< "("<<face->ptsList[2][0] <<","<<face->ptsList[2][1] <<","<<face->ptsList[2][2] <<"), "
		//					                 << face->ptIndxList[3]<< "("<<face->ptsList[3][0] <<","<<face->ptsList[3][1] <<","<<face->ptsList[3][2] <<")> "
		//			    << face->isBoundaryNode <<" "<< face->phyName <<"\n";

		oss.clear();
	}

	for (unsigned i = 0; i < mshBdsList.size(); i++) {
		vector<cFace*> mshFacesList = mshBdsList[i].mshFacesList;
		for (unsigned j = 0; j < mshFacesList.size(); j++) {
			cFace* face;
			face = mshFacesList[j];
			int nPts = face->nPts;
			ostringstream oss;
			oss << nPts << " ( ";
			for (int k = 0; k < nPts; k++) {
				oss << face->ptIndxList[k] << " ";
			}
			myFile << oss.str() << ")\n";
			//cout<< face->nid << " (" << face->ptIndxList[0] <<"," << face->ptIndxList[1]  <<"," << face->ptIndxList[2]  <<"," << face->ptIndxList[3]<< ") " << face->isBoundaryNode <<" "<< face->phyName <<"\n";

			oss.clear();
		}
	}
	myFile << ")" << endl;
	myFile.close();
}

void cOctreeOFMesh::saveAsOFMeshNeis() {
	string fileName = "./output/constant/polyMesh/neighbour";
	int num_mshIntlFaces = mshIntlFacesList.size();
	ofstream myFile;
	myFile.open(fileName, ios::ate | ios::out);
	if (!myFile.is_open()) {
		cout << "Open file failure" << endl;
	}
	//general mesh format
	myFile << "FoamFile" << endl;
	myFile << "{" << endl; //format_version Ascii_type data_size
	myFile << "    version     2.0;" << endl;
	myFile << "    format      ascii;" << endl;
	myFile << "    class       labelList;" << endl;
	myFile << "    note       \"nPoints:" << mshPtsList.size() << "\t nCells:" << octree->countAllNodes(&octree->root) << "\t nFaces:" << num_mshBFaces + num_mshIntlFaces << "\t nInternalFaces:" << num_mshIntlFaces << "\";\n";
	myFile << "    location    \"constant/polyMesh\";" << endl;
	myFile << "    object      neighbour;" << endl;
	myFile << "}" << endl;
	myFile << "" << endl;
	//points
	myFile << num_mshIntlFaces << endl;
	myFile << "(" << endl;
	for (int i = 0; i < num_mshIntlFaces; i++) {
		cFace* face;
		face = mshIntlFacesList[i];
		myFile << face->nbr << "\n";
	}
	myFile << ")" << endl;
	myFile.close();
}

void cOctreeOFMesh::saveAsOFMeshOwns() {
	string fileName = "./output/constant/polyMesh/owner";
	int num_mshIntlFaces = mshIntlFacesList.size();
	int num_mshBFaces = mshBFacesList.size();
	ofstream myFile;
	myFile.open(fileName, ios::ate | ios::out);
	if (!myFile.is_open()) {
		cout << "Open file failure" << endl;
	}
	//general mesh format
	myFile << "FoamFile" << endl;
	myFile << "{" << endl; //format_version Ascii_type data_size
	myFile << "    version     2.0;" << endl;
	myFile << "    format      ascii;" << endl;
	myFile << "    class       labelList;" << endl;
	myFile << "    note       \"nPoints:" << mshPtsList.size() << "\t nCells:" << octree->countAllNodes(&octree->root) << "\t nFaces:" << num_mshBFaces + num_mshIntlFaces << "\t nInternalFaces:" << num_mshIntlFaces << "\";\n";
	myFile << "    location    \"constant/polyMesh\";" << endl;
	myFile << "    object      owner;" << endl;
	myFile << "}" << endl;
	myFile << "" << endl;
	//points
	myFile << num_mshBFaces + num_mshIntlFaces << endl;
	myFile << "(" << endl;
	for (unsigned i = 0; i < mshIntlFacesList.size(); i++) {
		cFace* face;
		face = mshIntlFacesList[i];
		myFile << face->own << "\n";
	}

	for (unsigned i = 0; i < mshBdsList.size(); i++) {
		vector<cFace*> mshFacesList = mshBdsList[i].mshFacesList;
		for (unsigned j = 0; j < mshFacesList.size(); j++) {
			cFace* face;
			face = mshFacesList[j];
			myFile << face->own << "\n";
		}
	}

	myFile << ")" << endl;
	myFile.close();
}

void cOctreeOFMesh::saveAsOFMeshBds() {
	string fileName = "./output/constant/polyMesh/boundary";
	ofstream myFile;
	myFile.open(fileName, ios::ate | ios::out);
	if (!myFile.is_open()) {
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

	for (unsigned i = 0; i < mshBdsList.size(); i++) {
		string phyName = mshBdsList[i].phyName;
		string patchTpye = mshBdsList[i].patchType;
		int nFaces = mshBdsList[i].nFaces;
		int startFace = mshBdsList[i].startFace;
		ostringstream oss;
		oss << "\t" << phyName << "\n\t{\n";
		oss << "\t\t" << "type    " << patchTpye << ";\n";
		oss << "\t\t" << "nFaces    " << nFaces << ";\n";
		oss << "\t\t" << "startFace    " << startFace << ";\n";
		oss << "\t}\n";
		myFile << oss.str() << endl;
		oss.clear();
	}

	myFile << ")" << endl;
	myFile.close();
}


// +++++ save mesh