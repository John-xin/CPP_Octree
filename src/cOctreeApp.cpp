#include "cOctreeApp.h"
#include "cOctree.h"
#include "cOctreeUtil.h"

cOctreeApp::cOctreeApp()
{
	geom = cGeomData::getInstance();
	octree = new cOctree();
	util = new cOctreeUtil();
}

cOctreeApp::~cOctreeApp()
{
    cout << "Destroying the cOctreeApp" << endl;
}

void cOctreeApp::defineBody(vector<double> _ptInGeom) {
	if (_ptInGeom.size()!=3){
		cout<<"Must give a pt3D to define body"<<"\n";
	}else{
		geom->ptInGeom=_ptInGeom;
	}
}

void cOctreeApp::buildOctree() {

	
	//1. init root node
    octree->setup_root();

    //2. +++++++++++++++++++++++++++
	std::cout << "Min Octree Level is " << MIN_OCTREE_LEVELS << "\n";
	std::cout << "Max FeatPt in node is " << MAX_OCTNODE_FEATS << "\n\n";

	octree->splitOctreeByMinLevel(&(octree->root));
    //std::cout << "After splitOctreeByMinLevel " << "\n";
    //outputNodeName(&root);
	//util->output_octree("./output/byMinLevel_octree.txt", octree);
	
	octree->splitOctreeByFeaturePt(&(octree->root));
	//util->output_octree("./output/byFeaturePt_octree.txt", octree);
    //splitNodeByPhyName("bldg",6,&root);
    //splitNodeById("0-0");
    //std::cout << "After splitOctreeByFeaturePt " << "\n";
    //outputNodeName(&root);
	octree->getOctreeDepth(&(octree->root));
	octree->balanceOctree(&(octree->root));
	//util->output_octree("./output/balanced_octree.txt", octree);
	//util->outputNodeName(&(octree->root));
	octree->setup_leafNodesList(&(octree->root));
	octree->setup_leafNodesNbr();//n*logN - set nbrs for all leafNodes

	//3. identify boundary / interior / exterior node +++++++++++++++++++++++++++
	//find boundary node if it includes geoFFaces
	octree->setup_boundaryNode(&(octree->root)); //leaf nodes -> identify boundary / non-boundary node
	//closed STL -> node center - defined body pt -> ray -> intersections 
	octree->setup_interiorNode(&(octree->root)); //non-boundary node -> identify interior / exterior node

	//split nbrs fo extNode to same level as extNode

	//octree->splitExtNodeNbr(&(octree->root));
	//octree->setup_boundaryNode(&(octree->root));
	//octree->setup_interiorNode(&(octree->root));



	//octree->delExtNodes();
	//util->output_octree("./output/delExtNodes_octree.txt", octree);
    //std::cout<<"bNodesList length is "<< bNodesList.size() <<"\n";
    //std::cout<<"nonBNodesList length is "<< nonBNodesList.size() <<"\n";


    //4. +++++++++++++++++++++++++++
    //from leaf nodes -> identify non-repeated mshPt
	setup_mshPtList();
    cout<<"num of mshPts is "<< mshPtsList.size() <<"\n";
	util->output_leafNodes("./output/leafNodes.txt", octree->leafNodesList);
    //outputMshPts(root);
    //from leaf nodes -> define node -> 6 faces by mshPt index
    //init mshFace -> own, nid, low, upp, ptIndxList
    //and mark mshVolIndx
    setup_nodeMshFaces();
	//util->output_octree("./output/setup_mshFaces_octree.txt", octree);
   
	//5. from boundary / interior node -> identify internal / boundary mshFace +++++++++++++++++++++++++++
    //update own, nei, ptIndxList
    setup_boundaryMshFace();
    setup_intlMshFace();
    std::cout<<"mshIntlFacesList length is "<< mshIntlFacesList.size() <<"\n";
    std::cout<<"mshBFacesList length is "<< mshBFacesList.size() <<"\n\n";
    //outputMshFaces(root);
    //6. from boundary mshFace -> identify phyName of boundary mshFace +++++++++++++++++++++++++++
    setup_nodePhyName();
    setup_mshBdsList();

    //

}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeApp::setup_mshPtList()
{
	    for(unsigned i=0;i<octree->leafNodesList.size();i++){
	    	if(octree->leafNodesList[i]!=NULL){
	    		addMshPtsOfNode(octree->leafNodesList[i]);
	    	}
	    }
}
void cOctreeApp::addMshPtsOfNode(cOctNode* node) //init node.elePtLabels
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
    		nodesList= octree->getLeafNodeByPt(pt, &octree->root);
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
bool cOctreeApp::isPtSame(vector<double>& pt1, vector<double>& pt2)
{
    double tol = 0.000001;
    if ((abs(pt1[0] - pt2[0])) > tol) return false;
    if ((abs(pt1[1] - pt2[1])) > tol) return false;
    if ((abs(pt1[2] - pt2[2])) > tol) return false;
    return true;
}
void cOctreeApp::setup_nodeMshFaces() {
//    for(unsigned i=0;i<leafNodesList.size();i++){
//    	//mark mshVolIndx to all leaf node
//    	if(leafNodesList[i]->isInteriorNode!=0){
//        	leafNodesList[i]->mshVolIndx=mshVolIndxCount;
//        	leafNodesList[i]->updateMshFaces_mshVolIndx(mshVolIndxCount);//mark all node's mshFaces
//            leafNodesList[i]->calMshFaces();
//            mshVolIndxCount++;
//    	}else{
//    		leafNodesList[i]->mshVolIndx=-1;
//    		leafNodesList[i]->updateMshFaces_mshVolIndx(-1);//mark all node's mshFaces
//    	}
//
//    }
//

	for(unsigned i=0;i< octree->leafNodesList.size();i++){
    	//mark mshVolIndx to all leaf node

    	if(octree->leafNodesList[i]!= NULL){
			octree->leafNodesList[i]->mshVolIndx=mshVolIndxCount;
			octree->leafNodesList[i]->updateMshFaces_mshVolIndx(mshVolIndxCount);//mark all node's mshFaces
			octree->leafNodesList[i]->calMshFaces();
            mshVolIndxCount++;
    	}
    }
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeApp::setup_boundaryMshFace()
{
	//boundary node -> set boundary face
	//exterior node -> set boundary face
	//boundary node next to exterior node -> not set boundary face
	//cOctNode* node;
	for(unsigned int i = 0; i < octree->leafNodesList.size(); i++){
		cOctNode* &node= octree->leafNodesList[i]; //this is wrong when null is assign to reference

		//node=leafNodesList[i];

		if(node!=NULL){

		//if(node->nid=="0-0-3-3-3-3" ||node->nid=="0-0-3-3-3-7-0"){cout<<"";}
		for(unsigned j=0; j<6;j++ ){
			if(node->nbrsList[j].size()==0){
				node->mshFacesList[j].state=FaceState::boundary;
				mshBFacesList.push_back(&node->mshFacesList[j]);
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
void cOctreeApp::setup_intlMshFace()
{
	cFace* currentMshFace;
	cFace* nbrMshFace;
	vector<cOctNode*> nbr;
    string curFaceName;
    string nbrMshFaceName;
    int indxArr[6]={2, 3, 0, 1, 5, 4};

	for(unsigned int i = 0; i < octree->leafNodesList.size(); i++){
		cOctNode* &node= octree->leafNodesList[i];
		if(node!=NULL){

		for(unsigned j=0; j< node->mshFacesList.size();j++ ){
			currentMshFace=&node->mshFacesList[j];
			curFaceName=currentMshFace->nid;
			if(currentMshFace->state==FaceState::unassigned){
				nbr=node->nbrsList[j];
				if(nbr.size()!=0){

				for(unsigned k=0;k<nbr.size();k++){
					nbrMshFace=(&nbr[k]->mshFacesList[indxArr[j]]);
					nbrMshFaceName=nbrMshFace->nid;
					if(nbrMshFaceName=="node-0-0-3-3-3-3->face-4"){
						cout<<"";
					}
					int flag=currentMshFace->findFaceRelationship(currentMshFace, nbrMshFace);
					put2List(flag,currentMshFace,nbrMshFace);
				}

				}
			}
		}

		}
	}

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
    	if(mshBFacesList[i]->state==FaceState::internal){
    		num_IntlFaces_inMshBFacesList++;
    		cout<<mshBFacesList[i]->nid.c_str()<<"\n";
    	}else{
    		num_mshBFaces++;
    	}
    }
    cout<<"Number of Internal Faces in mshBFacesList is "<<num_IntlFaces_inMshBFacesList<<"\n";
    cout<<"Number of Boundary Faces in mshBFacesList is "<<num_mshBFaces<<"\n\n";

}
void cOctreeApp::findMshIntlFaces()
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
    	if(currentMshFace->state==FaceState::internal){
        	for (int j=k;j<numOfFaces;j++){
        		listMshFace=mshAllFacesList[j];
        		listFaceName=listMshFace->nid;
        		//if(listMshFace->isBoundaryNode!=1){
        			flag=currentMshFace->findFaceRelationship(currentMshFace, listMshFace);
        			put2List(flag,currentMshFace,listMshFace);
        			if(flag!=0){
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
        	bFlag=0;
    	}
    }
}
void cOctreeApp::put2List(int flag,cFace* currentMshFace,cFace* listMshFace)
{
    // 0 - not in
	// 1 - same face
	// 2 - f2 contains f1
	// 3 - f1 contains f2
 	if(flag==1){//same faces
        //keep only one internal face dir based on "small volLabel point to large volLabel"
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, then swap with the one in the list
        	currentMshFace->state=FaceState::internal;
        	listMshFace->state = FaceState::internal;
            currentMshFace->nbr = listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
        else {//current is larger, then keep it in the list
        	currentMshFace->state = FaceState::internal;
        	listMshFace->state = FaceState::internal;
        	listMshFace->nbr = currentMshFace->mshVolIndx; //current is large, then keep the list
        	mshIntlFacesList.push_back(listMshFace);
        }
	}else if(flag==2){// current face belongs to face in list-> then swap with the list
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, then keep anti-clock-wise order of the face
            listMshFace->state = FaceState::internal;//this is big face
            currentMshFace->state = FaceState::internal;
            currentMshFace->nbr = listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
        else {//current is larger, then change to clock-wise order of the current face
        	listMshFace->state = FaceState::internal;//this is big face
        	currentMshFace->state=FaceState::internal;
        	currentMshFace->changeOrder();
            currentMshFace->nbr=currentMshFace->mshVolIndx;
            currentMshFace->own=listMshFace->mshVolIndx;
            mshIntlFacesList.push_back(currentMshFace);
        }
	}else if(flag==3){//current face contains face in list -> keep the one in list
        if (currentMshFace->mshVolIndx < listMshFace->mshVolIndx) { //current is smaller, change order of the face in list
        	listMshFace->state = FaceState::internal;
        	currentMshFace->state = FaceState::internal;//this is big face

        	listMshFace->changeOrder();
            listMshFace->nbr = listMshFace->mshVolIndx;
            listMshFace->own = currentMshFace->mshVolIndx;
            mshIntlFacesList.push_back(listMshFace);

            //listMshFace=currentMshFace;
        }
        else {//current is larger, then keep the order of the face in list
        	listMshFace->state = FaceState::internal;
        	currentMshFace->state = FaceState::internal;//this is big face

        	listMshFace->nbr=currentMshFace->mshVolIndx;
            mshIntlFacesList.push_back(listMshFace);

            //listMshFace=currentMshFace;
        }
	}else if(flag==0){

	}
}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void cOctreeApp::setup_nodePhyName() {

	cOctNode *node=NULL;
	//vector<cTri>& geoFFacesList;
	cFace* mshFace=NULL;
    int numOfBFaces=0;

   //cout<<"BFaces in bNodesList "<<"\n";

	for(unsigned i=0; i< octree->leafNodesList.size();i++){
		node= octree->leafNodesList[i];
		if(node!=NULL){
			vector<cFeatureFace*> geoFFacesList=node->geoFFacesList;
			for(unsigned j=0; j<node->mshFacesList.size();j++){
				mshFace=&(node->mshFacesList[j]);
				//judge mshFace -> phyName
				if(mshFace->state==FaceState::boundary){
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
	cout<<"Number of PhyName assigned is "<<numOfBFaces<<"\n";
}
void cOctreeApp::setup_mshBdsList() {
	//init mshBdsList based on numOfPhynames
	mshBdsList.resize(geom->numOfPhyNames+1);
	for(int i=0; i< geom->numOfPhyNames;i++){
		mshBdsList[i].phyName= geom->phyNamesList[i];
		mshBdsList[i].phyNameIndx=i;
	}
	mshBdsList[geom->numOfPhyNames].phyName="unAssigned";
	mshBdsList[geom->numOfPhyNames].phyNameIndx=-100;

	//put mshFace into mshBdsList based on it's phyNameIndx
	cOctNode *node=NULL;
	cFace* mshFace=NULL;
	int indx;
	int numOfBFaces=0;
	for(unsigned int i=0; i< octree->leafNodesList.size();i++){
		node= octree->leafNodesList[i];
		if(node!=NULL){
			for(unsigned int j=0; j<node->mshFacesList.size();j++){
				mshFace=&(node->mshFacesList[j]);
				//judge mshFace -> phyName
				if(mshFace->state == FaceState::boundary){
					numOfBFaces++;
					if(mshFace->phyName!="unAssigned"){
						indx=mshFace->phyNameIndx;
						mshBdsList[indx].mshFacesList.push_back(mshFace);
					}else{
						mshBdsList[geom->numOfPhyNames].mshFacesList.push_back(mshFace);
						}
				}else{
					continue;
				}
			}
		}
	}

	//complete other properties of the List
	mshBdsList[0].startFace=mshIntlFacesList.size();
	mshBdsList[0].nFaces=mshBdsList[0].mshFacesList.size();
	for(int i=1; i< geom->numOfPhyNames +1;i++){
		mshBdsList[i].nFaces=mshBdsList[i].mshFacesList.size();
		mshBdsList[i].startFace=mshBdsList[i-1].nFaces+mshBdsList[i-1].startFace;
	}

	cout<<"Number of BFaces in mshBdsList is "<<numOfBFaces<<"\n\n";

}
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// +++++ save mesh
void cOctreeApp::saveAsGmsh(const char* _fileName)
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

void cOctreeApp::saveAsOFMesh() {
	saveAsOFMeshPts();
	saveAsOFMeshFaces();
	saveAsOFMeshNeis();
	saveAsOFMeshOwns();
	saveAsOFMeshBds();
}

void cOctreeApp::saveAsOFMeshPts() {
	string fileName="./output/constant/polyMesh/points";
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
			myFile << "( " << pt[0] << " " << pt[1] << " " << pt[2] << " )\n";
		}
		myFile << ")" << endl;
		myFile.close();
}

void cOctreeApp::saveAsOFMeshFaces() {
	string fileName="./output/constant/polyMesh/faces";
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
//			    << face->isBoundaryNode <<" "<< face->phyName <<"\n";

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
				//cout<< face->nid << " (" << face->ptIndxList[0] <<"," << face->ptIndxList[1]  <<"," << face->ptIndxList[2]  <<"," << face->ptIndxList[3]<< ") " << face->isBoundaryNode <<" "<< face->phyName <<"\n";

				oss.clear();
			}
		}
		myFile << ")" << endl;
		myFile.close();
}

void cOctreeApp::saveAsOFMeshNeis() {
	string fileName="./output/constant/polyMesh/neighbour";
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
		myFile << "    note       \"nPoints:" << mshPtsList.size()<< "\t nCells:"<< octree->countAllNodes(&octree->root) << "\t nFaces:"<<num_mshBFaces + num_mshIntlFaces<<"\t nInternalFaces:" <<num_mshIntlFaces<<"\";\n";
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

void cOctreeApp::saveAsOFMeshOwns() {
	string fileName="./output/constant/polyMesh/owner";
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
		myFile << "    note       \"nPoints:" << mshPtsList.size()<< "\t nCells:"<< octree->countAllNodes(&octree->root) << "\t nFaces:"<<num_mshBFaces + num_mshIntlFaces<<"\t nInternalFaces:" <<num_mshIntlFaces<<"\";\n";
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


void cOctreeApp::saveAsOFMeshBds() {
	string fileName="./output/constant/polyMesh/boundary";
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