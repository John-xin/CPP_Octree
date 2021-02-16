#ifndef COCTREEOFMESH_H
#define COCTREEOFMESH_H

#include "cFeature.h"
#include "cGeomData.h"
#include "cOctNode.h"
class cOctree;
using namespace std;

class cOctreeOFMesh {
public:
	cOctreeOFMesh();
	~cOctreeOFMesh();
	static cOctreeOFMesh* getInstance();
    cOctree* octree;

    //msh data
    vector<cOctNode*> mshNodesList;
    vector<vector<double>> mshPtsList;
    int mshPtIndxCount = 0;
    int mshVolIndxCount = 0;
    //vector<cFace> mshExtFacesList;
    vector<cFace*> mshIntlFacesList;
    vector<cFace*> mshBFacesList;
    vector<cFace*> mshAllFacesList;
    vector<cBoundary> mshBdsList;
    int num_mshBFaces, num_mshIntlFaces;
    //vector<int> eleVolList;

    void setup_octree(cOctree* _octree);
    void setup_mshNodesList();
    void setup_mshPtList(vector<cOctNode*> nodesList);
    void addMshPtsOfNode(cOctNode* node);
    bool isPtSame(vector<double>& pt1, vector<double>& pt2);
    void setup_nodeMshFaces();
    void update_nodeMshFaces();
    vector<int> getSurrPts(cOctNode* node);
    vector<cOctNode*> getSurrNodes(cOctNode* node);


    void setup_nodeMshFaceState();
    void setup_intlMshFace();
    void findMshIntlFaces();
    void put2List(int flag, cFace* currentMshFace, cFace* listMshFace);
    void saveInterlMshFaces(int flag, cFace* currentMshFace, cFace* listMshFace);

    void setup_nodePhyName();
    void setup_mshBdsList();

    void saveAsOFMesh();
    void saveAsOFMeshPts();
    void saveAsOFMeshFaces();
    void saveAsOFMeshNeis();
    void saveAsOFMeshOwns();
    void saveAsOFMeshBds();
};


#endif
