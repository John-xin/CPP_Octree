/*
 * commonFunc.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */
#include "commonFunc.h"
#include "cOctree.h"



//++++++++++++++++++++++++++++++++++++++
double commonFunc::dotProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates dot product v1.v2
    double dp=0.0;
    for (unsigned int i=0; i<3; i++)
        dp += v1[i]*v2[i];
    return dp;
}

double commonFunc::distBetweenPoints(vector<double> &p1, vector<double> &p2)
{
    // Calculate the distance between points p1 and p2, |p1-p2|
    double sum=0.0;
    for (unsigned int i=0; i<3; i++)
        sum += pow((p1[i]-p2[i]),2.0);
    sum = sqrt(sum);
    return sum;
}

commonFunc::commonFunc() {
}

commonFunc::~commonFunc() {
}

commonFunc* commonFunc::getInstance() {
	static commonFunc theCommonFunc;
//	if (theCommonFunc==NULL)
//    {
//        theCommonFunc = new commonFunc();
//    }
    return &theCommonFunc;
}

vector<double> commonFunc::crossProduct(vector<double> &v1, vector<double> &v2)
{
    // Calculates cross product v1xv2
    vector<double> cp(3);
    cp[0] = v1[1]*v2[2] - v1[2]*v2[1];
    cp[1] = v1[2]*v2[0] - v1[0]*v2[2];
    cp[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return cp;
}


vector<double> commonFunc::vectAdd( vector<double> &a, vector<double> &b )
{
    // Vector addition, c=a+b
    return vectAdd(a, b, 1.0);
}

vector<double> commonFunc::vectAdd( vector<double> &a, vector<double> &b, double sf )
{
    // Vector addition and scaling, c=a+sf*b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i] + sf*b[i];
    return c;
}

vector<double> commonFunc::vectSubtract( vector<double> &a, vector<double> &b )
{
    // Vector subtraction, c=a-b
    vector<double> c(a.size());
    for (unsigned int i=0; i<a.size(); i++)
        c[i] = a[i]-b[i];
    return c;
}

string commonFunc::NumberToString( int Number )
{
    // Converts integer to string
    ostringstream ss;
    ss << Number;
    return ss.str();
}

bool commonFunc::sortNodes(const pair<cOctNode*,double>&i, const pair<cOctNode*,double>&j)
{
    // Function used to sort a vector of cOctnode,double pairs by the value of
    // the double. The double will typically represent distance from the ray
    // origin in a ray-node intersection test
    return i.second < j.second;
}

void commonFunc::output_geoFPtsList(vector<cFeaturePt*>& geoFPtsList, vector<vector<double> >& geoPts3DList)
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

void commonFunc::output_octree(const char* _fileName, cOctree* tree)
{
    tree->allNodesList.resize(0);
    tree->update_allNodesList(&(tree->root));

    ofstream myFile;
    myFile.open(_fileName, ios::ate | ios::out);
    if (!myFile.is_open()) {
        cout << "Open file failure" << endl;
    }
    //**********output octree************
    myFile << "id numOfFPts numOfFEdges numOfFFaces mshPtsIndxList mshFaces-nid mshFaces-PtsIndxList" << endl;
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
        for (int j = 0; j < 6; j++) {
            mshFaces_ptIndxList[j].resize(4);
            for (int k = 0; k < 4;k++) {
                mshFaces_ptIndxList[j][k] = -1;
            }
        }

        if (node->mshFacesList.size() != 0) {
            for (unsigned int j = 0; j < node->mshFacesList.size(); j++) {
                if (node->mshFacesList[j].ptsList.size() != 0) {
                    for (unsigned int k = 0; k < node->mshFacesList[j].ptIndxList.size(); k++) {
                        mshFaces_ptIndxList[j][k] = node->mshFacesList[j].ptIndxList[k];
                    }
                } else {
                    mshFaces_ptIndxList[j].resize(4, 0);
                }
            }
        }

        ostringstream oss;
        oss << nid << " " << numOfGeoFPts << " " << numOfGeoFEdges << " " << numOfGeoFFaces << " [";
        for (int j = 0; j < 8; j++) {
            oss << mshPtsIndxList[j] << ",";
        }
        oss << "] [";

        for (int j = 0; j < 6; j++) {
            oss << mshFaces_nidList[j] << ",";
        }
        oss << "] [";

        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 4; k++) {
                oss << mshFaces_ptIndxList[j][k] << ",";
            }
            oss << "|";
        }
        oss << "]";
        myFile << oss.str() << endl;
        oss.clear();
    }
    //**********output octree************
    myFile.close();
}




