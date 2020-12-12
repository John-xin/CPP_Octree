#ifndef CGEOMDATA_H
#define CGEOMDATA_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>

using namespace std;

/* This Class read STL data.
 * Using pts3DList to store point coordinates <x y z>.
 * Using trisList to store STL triangle faces <ptIndx1 ptIndx2 ptIndx3 phyNameIndx>.
 * Using phyNamesList to store phyNames, i.e. solid name
 * low/upp -> two points of the bounding box of all points.
 *
 * Importance:
 * Keep one non-repeated pt3D list and mark an unique index for each pt3D.
 * Using the unique index to mark the Triangle faces.
 * Deal with empty line in the STL file.
 *
 * Solution:
 * Maintain a pts3DList. When having a pt3D, compare pt3D to each element of the list.
 * If it is existed, get same element index of the list; otherwise add it into the list and get the new index.
 *
*/
class cGeomData
{
    private:
        // Private Constructor
        cGeomData();
        // Stop the compiler generating methods of copy the object
        cGeomData(cGeomData const& copy);            // Not Implemented
        cGeomData& operator=(cGeomData const& copy); // Not Implemented
        ~cGeomData();

    public:
        static cGeomData* getInstance();
        void readSTLData(const char* _fileName);
        void addVert(vector<double> vert);
        bool isVertSame(vector<double>& vert1, vector<double>& vert2);
        void getLowUppVerts();//two points of the bounding box
        void output2Console();

        //fields
        //static cGeomData* theGeomData;
        vector<vector<double>> pts3DList; // n x 3
        vector<vector<int>> trisList; // n x 4 ->format (pt1, pt2, pt3, phyNameInx)
        vector<string> phyNamesList;
        int numOfTris;
        int numOfPts;
        int numOfPhyNames;
        vector<double> low, upp; //low & upp of the pts3DList

        //utility
        int countLabel=0;
        int phyNumCount=-1;
        vector<int> vLabels; //format (pt1, pt2, pt3, pyhNumInx)

};

//make it global var
//cGeomData* cGeomData::theGeomData;

#endif

