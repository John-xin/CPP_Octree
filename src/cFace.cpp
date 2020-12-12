/*
 * cFace.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#include "cFace.h"
#include "commonFunc.h"

cFace::cFace()
{
	nid="unAssigned";
    nbr=-100;
    own=-100;
    label=-100;
    nPts=0;
    bFlag=false;
    phyName = "unAssigned";
    phyNameIndx=-100;
    mshVolIndx = -100;
    state=-100;
    node=NULL;

}

cFace::~cFace(){}

void cFace::getN()
{
    // Get cTri face normal

    vector<vector<double> > v=ptsList; //3 x 3
    vector<double> v1(3),v2(3),v3;

    // vector 1 = vert 0 - vert 1 (dir 1->0)  ---- vector 2 = vert 0 - vert 2
    for (unsigned int i=0; i<3; i++) {
        v1[i] = v[0][i] - v[1][i];
        v2[i] = v[0][i] - v[2][i]; }
    v3 = commonFunc::getInstance()->crossProduct(v1,v2);
    double v3mag = sqrt(commonFunc::getInstance()->dotProduct(v3,v3));
    N = v3;
    for (vector<double>::iterator it=N.begin(); it!=N.end(); ++it)
        *it /= v3mag;
}

double cFace::getAngle(cFace* surface)
{
	double tol=0.000001;
    getN();
    surface->getN();
	double angle;
    angle=acos(commonFunc::getInstance()->dotProduct(N,surface->N))*180/3.141592653;
    if (fabs(angle)< tol){ angle=180.0;}
    return angle;

}

vector<double> cFace::getCentroid() {
	vector<double> pt;
	pt.resize(3);
	for(unsigned int i=0; i<ptsList.size();i++){
		pt[0]+=ptsList[i][0];
		pt[1]+=ptsList[i][1];
		pt[2]+=ptsList[i][2];
	}
	pt[0]/=ptsList.size();
	pt[1]/=ptsList.size();
	pt[2]/=ptsList.size();
	return pt;
}

void cFace::findPhyName(vector<cTri*> _geoFFacesList) {
	int indx1=-1;
	int indx2=-1;
	getN();
	vector<double> centroid, pt2;
	centroid=getCentroid();
	cLine ray1(centroid,N,1);
	indx1=findGeoFaceWithMinDist(ray1,_geoFFacesList);

//	pt2.resize(3);
//	pt2[0]=centroid[0]-0.01;
//	pt2[1]=centroid[1]-0.01;
//	pt2[2]=centroid[2]-0.01;
//	cLine ray2(low,N,1);
//	indx2=findGeoFaceWithMinDist(ray2,_geoFFacesList);

	ostringstream oss;

	if( indx1==-1){//no intersections
		phyNameIndx=-100;
		phyName="unAssigned";
	}else{
		phyNameIndx=_geoFFacesList[indx1]->phyNameIndx;
		phyName=_geoFFacesList[indx1]->phyName;

//		if(_geoFFacesList[indx1].phyName !=_geoFFacesList[indx2].phyName){
//			cout<<"Warning:" << nid << "-> phyName needs checking" <<"\n";
//			oss<<"Warning:" << nid << "-> phyName maybe geoFace-";
//			oss<<oss.str() << _geoFFacesList[indx1].indx << "-" <<_geoFFacesList[indx1].phyName.c_str() <<" or geoFace-" ;
//			oss<<oss.str() << _geoFFacesList[indx2].indx << "-" <<_geoFFacesList[indx2].phyName.c_str() <<"\n\n";
//			cout << oss.str();
//		}
	}

	//	if(nid=="node-0-3-6->face-0"){
	//		cout<<"selected index - " << indx <<"\n";
	//		for(unsigned int i=0; i<_geoFFacesList.size();i++){
	//			cout<< i <<" geoFace - " <<_geoFFacesList[i].indx << "-" << _geoFFacesList[i].phyName.c_str() << "\n" ;
	//		}
	//
	//	}

}

int cFace::findGeoFaceWithMinDist(cLine _ray, vector<cTri*> _geoFFacesList) {
	cTri* geoFFace;
	vector<double> pt;
	double s=0;
	double s_min=1e10;
	int indx=-1;
	bool flag;

	for(unsigned int i=0; i<_geoFFacesList.size();i++){
		geoFFace=_geoFFacesList[i];
		flag=geoFFace->rayPlaneIntersectPoint(_ray, pt, s);

		if(flag==true){
//			if(nid=="node-0-3-6->face-0"){
//				cout<< i <<" geoFace - " <<geoFFace.indx << " - " << geoFFace.phyName.c_str() << " - " << s <<"\n" ;
//			}
			if(fabs(s)<s_min){
				s_min=fabs(s);
				indx=i;
			}
		}
	}
	//if return -1 -> no intersections


//	if(nid=="node-0-3-6->face-0"){
//	cout<<"selected index - " << indx <<"\n";}

	return indx;
}

int numOfOrderChanged=0;

void cFace::changeOrder() {
	vector<int> _ptLabelList=ptIndxList;
	int numOfPts=ptIndxList.size();
	for(int i=0; i<numOfPts;i++){
		ptIndxList[i]=_ptLabelList[numOfPts-1-i];
	}
	numOfOrderChanged++;
	//cout<<"order changed -> "<< nid <<" ("<< ptIndxList[0] <<"," << ptIndxList[1]  <<"," << ptIndxList[2]  <<"," << ptIndxList[3]<<") \n";
}

int cFace::findFaceRelationship(cFace* f1, cFace* f2) {
	//assume f1 and f2 are rectangles
	double tol=1e-6;
	if(fabs(f1->getAngle(f2)-180)<tol){
		if(fabs(f1->low[0]-f2->low[0])<tol && fabs(f1->low[1]-f2->low[1])<tol && fabs(f1->low[2]-f2->low[2])<tol &&
		   fabs(f1->upp[0]-f2->upp[0])<tol && fabs(f1->upp[1]-f2->upp[1])<tol && fabs(f1->upp[2]-f2->upp[2])<tol)
		{
			return 1;//f1 = f2

		}else if(f1->isPtInFace(f2->low) && f1->isPtInFace(f2->upp))	{
			return 3;//f1>= f2

		}else if(f2->isPtInFace(f1->low) && f2->isPtInFace(f1->upp))	{
		    return 2;//f2>= f1
		}else{
			return 0;
		}
	}
	return 0;

}


bool cFace::isPtInFace(vector<double> pt) {
	//assume face is a rectangle
	double tol=1e-6;
	if(pt.size()==3){
	    if (low[0] - pt[0]> tol) return false;
	    if (low[1] - pt[1]> tol) return false;
	    if (low[2] - pt[2]> tol) return false;
	    if (upp[0] - pt[0]< -tol) return false;
	    if (upp[1] - pt[1]< -tol) return false;
	    if (upp[2] - pt[2]< -tol) return false;
	    return true;
	}else{
		cout<<"error: pt size in <isPtInFace Func> is not 3" << "\n";
		return false;
	}
}
