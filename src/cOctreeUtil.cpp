#include "cOctree.h"
#include "cOctNode.h"
#include "cFeature.h"
#include "cOctreeUtil.h"

cOctreeUtil::cOctreeUtil(){

}

cOctreeUtil::~cOctreeUtil(){

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
