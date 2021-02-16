// Octree_20201212.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include "src/cOctreeApp.h"


using namespace std;

int main()
{
    const char* fileName = "./etc/APR1st_Cube_Case_combined.stl"; //or const char* fileName ="F:\\Project\\cpluplus\\Eclipse_Proj\\Octree\\Octree\\cube_1m.stl";  
 
    cOctreeApp* app = new cOctreeApp();

    //+++++++ cGeomData geom is initialized +++++++++++++
    app->geom->readSTLData(fileName);
    app->geom->extractFeature();
    //vector<double> pt3D{ 0.15,0.5,0.1 };
    vector<double> pt3D{ 7,7,7 };
    app->defineBody(pt3D);
    //+++++++ cGeomData geom is initialized +++++++++++++

    app->buildOctree();
    app->buildOFMesh();
    app->saveAsOFMesh();
    //commonFunc->output_geoFPtsList(myTree->geoFPtsList, myTree->geoPts3DList);
    ////std::cout<<"NumOfNodes is "<< myTree->countNodes(myTree->root) <<"\n";
    //extern int numOfOrderChanged;
    //cout << "number of order changed is " << numOfOrderChanged << "\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
