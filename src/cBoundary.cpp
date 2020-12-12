/*
 * cBoundary.cpp
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */
#include "cBoundary.h"

cBoundary::cBoundary()
{
    // Default cTri constructor
    phyName="unAssigned";
    phyNameIndx=-100;
    patchType="unAssigned";
    startFace=-100;
    nFaces=0;
}


// cTri destructor
cBoundary::~cBoundary() {
    //cout << "Destroying cTri" << endl;
}



