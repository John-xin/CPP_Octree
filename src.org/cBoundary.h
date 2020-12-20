#ifndef CBOUNDARY_H
#define CBOUNDARY_H

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
//#include "cGeomData.h"
#include "cOctNode.h"
#include "cFace.h"

// OpenMP headers
#ifdef _OPENMP
  #include <omp.h>
#endif

class cBoundary {
public:
    cBoundary();
    ~cBoundary();
    
    string phyName;
    int phyNameIndx;
    vector<cFace*> mshFacesList;
    string patchType;
    int startFace;
    int nFaces;
};


#endif
