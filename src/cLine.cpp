/*
 * cLine.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */

#include "cLine.h"
#include <math.h>

cLine::cLine()
{
    // Default constructor for cLine
    // Default line is unit vector along x-axis
    p0.resize(3,0.0); p1.resize(3,0.0); dir.resize(3,0.0);
    p1[0]=1.0; dir[0]=1.0;
}

cLine::cLine(vector<double> &_p0, vector<double> &p1_dir, int isP1orDir)
{
    // cLine constructor with p0 and p1 or dir
    // if isP1orDir==0, then p1_dir is p1
    // if isP1orDir==1, then p1_dir is dir
    p0 = _p0;
    if (isP1orDir==0) {
        p1 = p1_dir;
        getDir(); }
    else if (isP1orDir==1) {
        dir = p1_dir;
        getP1(); }
}

// cLine destructor
cLine::~cLine() {}

void cLine::getDir()
{
    // Get unit vector defining direction of cLine
    vector<double> p0p1(3); double dmag=0.0;
    for (unsigned int i=0; i<3; i++) {
        p0p1[i] = p1[i]-p0[i];
        dmag   += pow(p0p1[i],2.0);
    }
    dmag = sqrt(dmag);
    dir  = p0p1;
    for (vector<double>::iterator it=dir.begin(); it!=dir.end(); ++it)
        *it /= dmag;
}

void cLine::getP1()
{
    // Get a point on the cLine, p1, located 1.0 units away from the origin, p0
    vector<double> p1(3);
    for (unsigned int i=0; i<3; i++)
        p1[i] = p0[i]+dir[i];
}

// ------------------------------------------------------



