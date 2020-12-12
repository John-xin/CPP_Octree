/*
 * commonFunc.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */
#include "commonFunc.h"

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



