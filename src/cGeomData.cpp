/*
 * cGeomData.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#include "cGeomData.h"
#include "cfeature.h"
#include <math.h>

//++++++++++++ utility functions
template <typename Out>
void split(const string &s, char delim, Out result) {
    istringstream iss(s);
    string item;
    while (getline(iss, item, delim)) {
        *result++ = item;
    }
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, back_inserter(elems));
    return elems;
}

double stringToDouble(const string& str)
{
istringstream iss(str);
double num;
iss >> num;
return num;
}
//++++++++++++

cGeomData::cGeomData()
{
	numOfTris=-1;
	numOfPts=-1;
	numOfPhyNames=-1;

    numOfGeoFFaces = -1;
    numOfGeoFEdges = -1;
    numOfGeoFPts = -1;

	//vertCoords3D=NULL;
	//triLabelMtx=NULL;
}

cGeomData:: ~cGeomData()
{
}

cGeomData* cGeomData::getInstance()
{
    // The only instance
    // Guaranteed to be lazy initialized
    // Guaranteed that it will be destroyed correctly
	static cGeomData theGeomData;
//    if (theGeomData==NULL)
//    {
//        theGeomData=new cGeomData();
//    }
    return &theGeomData;
}

void cGeomData::readSTLData(const char* _fileName)
{
  //prerequisite:
  //make sure that there is no empty space line in the STL file;

  string line;
  vector<string> tmp;
  int count=0;
  vector<double> v(3,0.0);

  ifstream myfile (_fileName);
  if (myfile.is_open())
  {
    while ( getline (myfile,line) )
    {
      tmp=split(line, ' ');
      if(tmp.size()<2){continue;}//the line content less than 2, skip this line

      if (tmp[0]=="solid"){
    	  phyNamesList.push_back(tmp[1]);
          count=-1;
          phyNumCount++;
      }

      if(count==2 || count==3 || count==4){
          v[0]=stringToDouble(tmp[7]);
          v[1]=stringToDouble(tmp[8]);
          v[2]=stringToDouble(tmp[9]);
          addVert(v);
      }

      if(count==6){
          count=-1;
          vLabels.push_back(phyNumCount);//put phyNumIndx at the last (pt1, pt2, pt3, pyhNumInx)
          trisList.push_back(vLabels);
          vLabels.clear();
      }

      count++;

      //for(int i=0;i<tmp.size();i++) {cout << tmp[i] << "@";}

      //cout << '\n';

    }

    /*
    cout << "num of tris " << triLabelMtx.size() << '\n';
    for(int i=0;i<triLabelMtx.size();i++)
    {
        for(int j=0;j<triLabelMtx[i].size();j++)
        {
            cout << triLabelMtx[i][j] << "@\t";
        }
        cout << '\n';
    }

    cout << "num of verts " << vertCoords3D.size() << '\n';
    for(int i=0;i<vertCoords3D.size();i++)
    {
        for(int j=0;j<vertCoords3D[i].size();j++)
        {
            cout << vertCoords3D[i][j] << "@\t";
        }
        cout << '\n';

    }
    */



    myfile.close();
  }

  else cout << "Unable to open STL file" << "\n";

  //init fields
  numOfTris=trisList.size();
  numOfPts=pts3DList.size();
  numOfPhyNames=phyNamesList.size();

  getLowUppVerts();
  output2Console();
}

void cGeomData::getLowUppVerts() {
	low.resize(3);
	upp.resize(3);
	for (unsigned int i=1; i<pts3DList.size(); i++) {
	        for (int j=0; j<3; j++) {
	            if (pts3DList[i][j] < low[j]) { low[j] = pts3DList[i][j]; }
	            if (pts3DList[i][j] > upp[j]) { upp[j] = pts3DList[i][j]; }
	        }
	    }
}

bool cGeomData::isVertSame(vector<double>& vert1, vector<double>& vert2)
{
    double tol = 1e-6;
    if ((fabs(vert1[0] - vert2[0])) > tol) return false;
    if ((fabs(vert1[1] - vert2[1])) > tol) return false;
    if ((fabs(vert1[2] - vert2[2])) > tol) return false;
    return true;

}

void cGeomData::addVert(vector<double> vert)
{
    int vertLabel=-1;
    //check if this pt3D already existed in list.
    for (unsigned int i = 0; i < pts3DList.size(); i++) {
        if (isVertSame(vert,pts3DList[i])) {
            vertLabel = i;
            break; //same pt3D found in the list, stop the loop.
        }
    }
    //this is a new pt3D -> add it into list and get the new index to vLabels
    if (vertLabel == -1) {
        pts3DList.push_back(vert);
        vLabels.push_back(countLabel);
        countLabel++;
    }else
    //not new vert -> find the same element index in list
    //put the found index into vLabels
    {
        vLabels.push_back(vertLabel);
    }

}

void cGeomData::output2Console() {
	cout << "Number of STL points is " << numOfPts  << "\n";
	cout << "Number of STL Triangles is " << numOfTris  << "\n";
	cout << "Number of STL SolidNames is " << numOfPhyNames  << "\n";
	cout << "\n";
}

// ++++++ setup feature list +++++++++++++++++++
void cGeomData::setup_geoFFacesList()
{
    vector<int> geoPtIndx_int3(3);
    int phyNameIndx;
    int geoFFaceIndx;
    geoFFacesList.reserve(numOfTris);
    for (int i=0; i< numOfTris; i++) {
        geoPtIndx_int3[0]= trisList[i][0];
        geoPtIndx_int3[1]= trisList[i][1];
        geoPtIndx_int3[2]= trisList[i][2];
        phyNameIndx= trisList[i][3];
        geoFFaceIndx=i;
        cFeatureFace* fFace=new cFeatureFace(geoFFaceIndx,phyNameIndx,geoPtIndx_int3);
        //if(fFace->indx==1060 ||fFace->indx==1061){cout<<"";} connect at two FFaces intersection line
        geoFFacesList.push_back(fFace);
    }
    std::cout<<"number of geoFFaces is " <<geoFFacesList.size()<<"\n";
}

void cGeomData::setup_geoFEdgesList()
{
    vector<int> geoPtIndx_int2;
    int geoFEdgeIndx=0;
    cFeatureFace* surface1;
	cFeatureFace* surface2;

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

void cGeomData::setup_geoFPtsList()
{

    for (unsigned int i=0; i<geoFEdgesList.size(); i++) {
    	for (int j=0; j<2; j++) {
    		cFeaturePt* fPt=new cFeaturePt();
    		fPt->indx=geoFEdgesList[i]->ptIndx_int2[j];
    		if(fPt->isNewInList(geoFPtsList)){
    			geoFPtsList.push_back(new cFeaturePt(fPt->indx));
    		}
        }
    }
    std::cout<<"number of feature points is " <<geoFPtsList.size()<<"\n\n";
}

void cGeomData::extractFeature()
{
    //extract features from geom data
    setup_geoFFacesList(); //first step
    setup_geoFEdgesList(); //find feature edge from setup_geoFFacesList
    setup_geoFPtsList();  //find feature point from setup_geoFEdgesList

    numOfGeoFFaces = geoFFacesList.size();
    numOfGeoFEdges = geoFEdgesList.size();
    numOfGeoFPts = geoFPtsList.size();
}
// ++++++ setup feature list ++++++++++++++++++++++++++
