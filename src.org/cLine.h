/*
 * cLine.h
 *
 *  Created on: 24 Jan, 2020
 *      Author: WangXJ
 */


#ifndef CLINE_H_
#define CLINE_H_
//#include <stdlib.h>
//#include <stdio.h>
#include <iostream>
#include <vector>
//#include <set>
//#include <math.h>
//#include <string>
//#include <sstream>
//#include <algorithm> // find, sort
//#include <utility>   // pair

using namespace std;

class cLine {
public:
    vector<double> p0,p1,dir;
    cLine();
    cLine(vector<double> &lp0, vector<double> &p1_dir, int isp1orDir);
    ~cLine();
    void getDir();
    void getP1();
};

//++++++++++++++++++++++++++++++++++++++


#endif /* CLINE_H_ */
