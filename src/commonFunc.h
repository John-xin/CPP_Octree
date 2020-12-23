/*
 * commonFunc.h
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#ifndef COMMONFUNC_H_
#define COMMONFUNC_H_

#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm> // find, sort
#include <utility>   // pair
#include <fstream>

using namespace std;

class cOctree;
class cOctNode;
class cFeatureFace;
class cFeaturePt;

// Function prototypes

class commonFunc
{
public:
	commonFunc();
	~commonFunc();
	static commonFunc* getInstance();

	bool sortNodes(const pair<cOctNode*,double>&i, const pair<cOctNode*,double>&j);
	double distBetweenPoints(vector<double> &p1, vector<double> &p2);
	string NumberToString( int Number );
	vector<double> vectAdd( vector<double> &a, vector<double> &b);
	vector<double> vectAdd( vector<double> &a, vector<double> &b, double sf);
	vector<double> vectSubtract( vector<double> &a, vector<double> &b );
	double dotProduct( vector<double> &v1, vector<double> &v2 );
	vector<double> crossProduct( vector<double> &v1, vector<double> &v2 );

	//static commonFunc* theCommonFunc;


};



//++++++++++++++++++++++++++++++++++++++
typedef struct intersection
{
    int triLabel;
    double s;
    vector<double> p;
    intersection() { triLabel=0; p.resize(3,0.0); s=0.0; }
    intersection(int _triLabel, vector<double> _p, double _s) { triLabel=_triLabel; p=_p; s=_s;}
    // Overload operator < so we can use sort function of a std::vector
    bool operator < (const intersection& intersect) const {
        return (s < intersect.s); }
} Intersection;



#endif /* COMMONFUNC_H_ */
