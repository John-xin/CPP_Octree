#include "cOctreeApp.h"
#include "cOctree.h"
#include "cOctreeUtil.h"
#include "cOctreeOFMesh.h"

cOctreeApp::cOctreeApp()
{
	geom = cGeomData::getInstance();
	octree = new cOctree();
	util = new cOctreeUtil();
	util->setup_octree(octree);
	ofMesh = cOctreeOFMesh::getInstance();
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
	std::cout << "Min Octree Level is " << octree->MIN_OCTREE_LEVELS << "\n";
	std::cout << "Max Octree Level is " << octree->MAX_OCTREE_LEVELS << "\n";
	std::cout << "Max FeatPt in node is " << octree->MAX_OCTNODE_FEATS << "\n\n";
	//1. init root node
    octree->setup_root();

    //2. +++++++++++++++++++++++++++
	octree->splitOctreeByMinLevel(&(octree->root));
	util->output_octree("./output/byMinLevel_octree.txt", octree);
	//octree->getNbr18Nodes(&(octree->root));
	octree->splitOctreeByFeaturePt(&(octree->root));
	//util->output_octree("./output/byFeaturePt_octree.txt", octree);
    //splitNodeByPhyName("bldg",6,&root);
	octree->getOctreeDepth(&(octree->root));
	//octree->balanceOctree(&(octree->root));
	//3. identify boundary / interior / exterior node +++++++++++++++++++++++++++
    //find boundary node if it includes geoFFaces
	//octree->setup_boundaryNode(&(octree->root)); //leaf nodes -> identify boundary / non-boundary node
    //closed STL -> node center - defined body pt -> ray -> intersections 
	//octree->setup_interiorNode(&(octree->root)); //non-boundary node -> identify interior / exterior node
	//util->outputNodeName(&(octree->root));
	octree->setup_leafNodesList(&(octree->root));
	//octree->setup_leafNodesNbr();//n*logN - set nbrs for all leafNodes
    
    util->checkNodes_ErrorType1(octree->leafNodesList);

	//*******note this********
	//octree->setup_boundaryNode(&(octree->root));
	//octree->setup_interiorNode(&(octree->root));
	//octree->setup_leafNodesNbr();



	//split nbrs fo extNode to same level as extNode
	//octree->splitExtNodeNbr(&(octree->root));
	//octree->delExtNodes();
	//util->output_octree("./output/delExtNodes_octree.txt", octree);
}

void cOctreeApp::buildOFMesh()
{
	ofMesh->setup_octree(octree);
	ofMesh->setup_mshNodesList();
	cout << "num of mshNodes is " << ofMesh->mshNodesList.size() << "\n";
	//4. +++++++++++++++++++++++++++
    //from leaf nodes -> identify non-repeated mshPt
	ofMesh->setup_mshPtList(octree->leafNodesList);
	//ofMesh->setup_mshPtList(ofMesh->mshNodesList);

	//ofMesh->setup_mshPtList(ofMesh->mshNodesList);
	cout << "num of mshPts is " << ofMesh->mshPtsList.size() << "\n";
	//from leaf nodes -> define node -> 6 faces by mshPt index
    //init mshFace -> own, nid, low, upp, ptIndxList
    //and mark mshVolIndx
	ofMesh->setup_nodeMshFaces();//setup mshFaces with 4 pts
	ofMesh->update_nodeMshFaces();//update mshFaces with actual pts
	ofMesh->setup_nodeMshFaceState();

	//std::cout<<"bNodesList length is "<< bNodesList.size() <<"\n";
	//std::cout<<"nonBNodesList length is "<< nonBNodesList.size() <<"\n";

	util->output_octree("./output/setup_mshFaces_octree.txt", octree);

	//5. from boundary / interior node -> identify internal / boundary mshFace +++++++++++++++++++++++++++
	//update own, nei, ptIndxList
	//setup_boundaryMshFace();
	util->output_nodes("./output/leafNodes.txt", octree->leafNodesList);
	util->output_nodes("./output/mshNodes.txt", ofMesh->mshNodesList);
	std::cout << "mshIntlFacesList length is " << ofMesh->mshIntlFacesList.size() << "\n";
	std::cout << "mshBFacesList length is " << ofMesh->mshBFacesList.size() << "\n\n";
	//outputMshFaces(root);
	//6. from boundary mshFace -> identify phyName of boundary mshFace +++++++++++++++++++++++++++
	ofMesh->setup_nodePhyName();
	ofMesh->setup_mshBdsList();
}

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
	ofMesh->saveAsOFMeshPts();
	ofMesh->saveAsOFMeshFaces();
	ofMesh->saveAsOFMeshNeis();
	ofMesh->saveAsOFMeshOwns();
	ofMesh->saveAsOFMeshBds();
}