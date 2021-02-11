#include "cOctree.h"
#include "cOctNode.h"
#include "cFeature.h"
#include "cOctreeUtil.h"



cOctreeUtil::cOctreeUtil(){

}

cOctreeUtil::~cOctreeUtil(){

}

void cOctreeUtil::setup_octree(cOctree* _octree)
{
    octree = _octree;
}


void cOctreeUtil::outputNodeName(cOctNode* node) {
    if (node->isLeafNode()) {
        cout << node->nid.c_str() << "\n";
    }
    else {
        for (unsigned int i = 0; i < node->children.size(); i++) {
            outputNodeName(node->children[i]);
        }
    }
}

void cOctreeUtil::outputMshPts(cOctNode* node) {
    vector<int> ptIndxList;
    if (node->isLeafNode()) {
        ptIndxList = node->mshPtsIndxList;
        cout << node->nid.c_str() << " state:" << (node->state) << " ptIndxList (" << ptIndxList[0] << "," << ptIndxList[1] << "," << ptIndxList[2] << "," << ptIndxList[3] << "," << ptIndxList[4] << "," << ptIndxList[5] << "," << ptIndxList[6] << "," << ptIndxList[7] << ")" << "\n";

    }
    else {
        for (unsigned int i = 0; i < node->children.size(); i++) {
            outputMshPts(node->children[i]);
        }
    }
}

void cOctreeUtil::outputMshFaces(cOctNode* node) {
    vector<cFace> faceList;
    cFace face;
    if (node->isLeafNode()) {
        faceList = node->mshFacesList;
        for (unsigned j = 0; j < faceList.size(); j++) {
            face = faceList[j];
            cout << face.nid.c_str() << " state:" << face.state << " ptIndxList (" << face.ptIndxList[0] << "," << face.ptIndxList[1] << "," << face.ptIndxList[2] << "," << face.ptIndxList[3] << ")" << "\n";

        }
    }
    else {
        for (unsigned int i = 0; i < node->children.size(); i++) {
            outputMshFaces(node->children[i]);
        }
    }
}

void cOctreeUtil::output_geoFPtsList(vector<cFeaturePt*>& geoFPtsList, vector<vector<double> >& geoPts3DList)
{
    string fileName = "./output/featruePts.txt";
    ofstream myFile;
    myFile.open(fileName, ios::ate | ios::out);
    if (!myFile.is_open()) {
        cout << "Open file failure" << endl;
    }
    //**********output pts************
    myFile << "number x y z index" << endl;
    //output pts
    int count = 1;
    for (unsigned i = 0; i < geoFPtsList.size(); i++) {
        cFeaturePt* fPt = geoFPtsList[i];
        vector<double> pt = geoPts3DList[fPt->indx];
        ostringstream oss;
        oss << count << " " << pt[0] << " " << pt[1] << " " << pt[2] << " " << fPt->indx;
        myFile << oss.str() << endl;
        oss.clear();
        count++;
    }
    //**********output pts************
    myFile.close();
}

void cOctreeUtil::output_octree(const char* _fileName, cOctree* tree)
{
    tree->allNodesList.resize(0);
    tree->update_allNodesList(&(tree->root));

    ofstream myFile;
    myFile.open(_fileName, ios::ate | ios::out);
    if (!myFile.is_open()) {
        cout << "Open file failure" << endl;
    }
    //**********output octree************
    myFile << "id numOfFeatures mshPtsIndxList mshFaces-nid mshFaces-PtsIndxList" << endl;
    //output 
    for (unsigned i = 0; i < tree->allNodesList.size(); i++) {
        cOctNode* node = tree->allNodesList[i];
        string nid = node->nid;
        int numOfGeoFPts = node->geoFPtsList.size();
        int numOfGeoFEdges = node->geoFEdgesList.size();
        int numOfGeoFFaces = node->geoFFacesList.size();

        vector<int> mshPtsIndxList;
        mshPtsIndxList.resize(8, -1);
        if (node->mshPtsIndxList.size() != 0) {
            for (unsigned int j = 0; j < node->mshPtsIndxList.size(); j++) {
                mshPtsIndxList[j]=node->mshPtsIndxList[j];
            }
        }

        vector<string> mshFaces_nidList;
        mshFaces_nidList.resize(6, "-1");
        for (int j = 0; j < 6; j++) {
            mshFaces_nidList[j] = node->mshFacesList[j].nid;
        }

        vector<vector<int>> mshFaces_ptIndxList;
        mshFaces_ptIndxList.resize(6);
  
        if (node->mshFacesList.size() != 0) {
            for (unsigned int j = 0; j < node->mshFacesList.size(); j++) {
                if (node->mshFacesList[j].ptsList.size() != 0) {
                    for (unsigned int k = 0; k < node->mshFacesList[j].ptIndxList.size(); k++) {
                        mshFaces_ptIndxList[j].push_back(node->mshFacesList[j].ptIndxList[k]);
                    }
                } else {
                    mshFaces_ptIndxList[j].resize(0);
                }
            }
        }

        ostringstream oss;
        oss << nid << " [ " << numOfGeoFFaces << "," << numOfGeoFEdges << "," <<  numOfGeoFPts << " ] [ ";
        for (int j = 0; j < 8; j++) {
            oss << mshPtsIndxList[j] << ",";
        }
        oss.seekp(-1, oss.cur);
        oss << " ] [ ";

        for (int j = 0; j < 6; j++) {
            oss << mshFaces_nidList[j] << ",";
        }
        oss.seekp(-1, oss.cur);
        oss << " ] [ ";

        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < mshFaces_ptIndxList[j].size(); k++) {
                oss << mshFaces_ptIndxList[j][k] << ",";
            }
            oss.seekp(-1, oss.cur);
            oss << " | ";
        }
        oss.seekp(-2, oss.cur);
        oss << "]";
        myFile << oss.str() << endl;
        oss.clear();
    }
    //**********output octree************
    myFile.close();
}

void cOctreeUtil::output_nodes(const char* _fileName, vector<cOctNode*> nodesList)
{
    ofstream myFile;
    myFile.open(_fileName, ios::ate | ios::out);
    if (!myFile.is_open()) {
        cout << "Open file failure" << endl;
    }
    //**********output octree************
    myFile << "id mshVolIndx mshPtIndxList node-state numOfFeatures mshFaces-state mshFaces-ptsIndx" << endl;
    //output 
    for (unsigned int i = 0; i < nodesList.size(); i++) {
        cOctNode* node = nodesList[i];

        string nid = node->nid;
        int numOfGeoFFaces = node->geoFFacesList.size();
        int numOfGeoFEdges = node->geoFEdgesList.size();
        int numOfGeoFPts = node->geoFPtsList.size();

        vector<int> mshPtsIndxList;
        mshPtsIndxList.resize(8, -1);
        if (node->mshPtsIndxList.size() != 0) {
            for (unsigned int j = 0; j < node->mshPtsIndxList.size(); j++) {
                mshPtsIndxList[j] = node->mshPtsIndxList[j];
            }
        }

        vector<string> mshFaces_stateList;
        mshFaces_stateList.resize(6, "-1");
        for (int j = 0; j < 6; j++) {
            mshFaces_stateList[j] = to_string(node->mshFacesList[j].state);
        }

        vector<vector<int>> mshFaces_ptIndxList;
        mshFaces_ptIndxList.resize(6);
        for (int j = 0; j < 6; j++) {
            mshFaces_ptIndxList[j].resize(0);
        }

        if (node->mshFacesList.size() != 0) {
            for (unsigned int j = 0; j < node->mshFacesList.size(); j++) {
                if (node->mshFacesList[j].ptsList.size() != 0) {
                    for (unsigned int k = 0; k < node->mshFacesList[j].ptIndxList.size(); k++) {
                        mshFaces_ptIndxList[j].push_back(node->mshFacesList[j].ptIndxList[k]);
                    }
                }
                else {
                    mshFaces_ptIndxList[j].resize(0);
                }
            }
        }

        ostringstream oss;
        oss << nid << " [ ";
        for (int j = 0; j < 8; j++) {
            oss << mshPtsIndxList[j] << ",";
        }
        oss.seekp(-1, oss.cur);
        oss << " ] ";
        oss << node->state;

        oss << " [ " << numOfGeoFFaces << "," << numOfGeoFEdges << "," << numOfGeoFPts;
        oss << " ] [ ";

        for (int j = 0; j < 6; j++) {
            oss << mshFaces_stateList[j] << ",";
        }
        oss.seekp(-1, oss.cur);
        oss << " ] [ ";

        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < mshFaces_ptIndxList[j].size(); k++) {
                oss << mshFaces_ptIndxList[j][k] << ",";
            }
            oss.seekp(-1, oss.cur);
            oss << " | ";
        }
        oss.seekp(-2, oss.cur);
        oss << "]";

        myFile << oss.str() << endl;
        oss.clear();
    }
    //**********output octree************
    myFile.close();

}


void cOctreeUtil::checkNodes_ErrorType1(vector<cOctNode*> nodesList)
{
    int errorNodesCount = 0;
    queue<cOctNode*> possibleErrorNodes;
    for (cOctNode* node : nodesList) {
        possibleErrorNodes.push(node);
    }

    while (possibleErrorNodes.size()!=0) {
        cOctNode* myNode = possibleErrorNodes.front();
        if (myNode->isLeafNode()) {
            if (isErrorType1(myNode)) {
                vector<cOctNode*> introducedNodes = fix_ErrorType1(myNode);
                errorNodesCount++;
                for (cOctNode* node : introducedNodes) {
                    possibleErrorNodes.push(node);
                }
            }
        }
        possibleErrorNodes.pop();
    }
   
    cout << "Check Error Type1: Number of error nodes are " << errorNodesCount <<"\n";
    octree->update_leafNodesList(&(octree->root));
}

bool cOctreeUtil::isErrorType1(cOctNode* node)
{
    //Error Type 1 : Node's nbrs exist boundary nodes and exterior nodes in the same time. 
    int faceCount = -1;
    bool flag = false;

    for (vector<cOctNode*> nbr : node->nbrsList) {
        faceCount++;
        int numOfExtNodes = 0;
        int numOfBNodes = 0;
        if (nbr.size() == 4) {
            for (cOctNode* nbrNode : nbr) {
                if (nbrNode->state == NodeState::exterior) {
                    numOfExtNodes++;
                }
                else if (nbrNode->state == NodeState::boundary) {
                    numOfBNodes++;
                }
            }
        }
        if (numOfExtNodes!=0 && numOfBNodes!=0) {
            //cout << "Error: Node-" << node->nid.c_str() << " has exterior && boundary nbrNodes at face-" << faceCount << "\n";
            flag = true;
        }
    }

    return flag;
}


vector<cOctNode*>  cOctreeUtil::fix_ErrorType1(cOctNode* node)
{
    //fix one error node would introduce at least 6 possible error nodes from 6 nbrs
    vector<cOctNode*> possibleErrorNodes;
    if (node->isLeafNode()) {
        //1. Address error node by split the node to same level of the nbr node
        octree->splitNode(node);

        //2. Setup the new split nodes
        for (cOctNode* child : node->children) {
            octree->setup_boundaryNode(child);
            octree->setup_interiorNode(child);
            octree->setLeafNodeNbr(child);
        }

        //3. Repair the nbr relationship of the split node's nbrs and introduce new possible error nodes to check
        vector<cOctNode*> nbr6Nodes = octree->getNbr6Nodes(node);
        for (cOctNode* myNode : nbr6Nodes) {
            octree->setLeafNodeNbr(myNode);
            if (myNode->isLeafNode()) {
                possibleErrorNodes.push_back(myNode);
            }
        }
    }
    return possibleErrorNodes;
}

