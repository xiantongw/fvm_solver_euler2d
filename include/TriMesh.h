#ifndef TRIMESH_H
#define TRIMESH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include "../include/utils.h"

using namespace std;

class TriMesh
{
  public:
    string gri_filename;
    double num_node;
    int num_element;
    double num_boundary;
    int curved_group;
    vector<vector<int> > E;
    vector<vector<double> > V;
    vector<vector<vector<int> > > B;
    vector<double> Area;
    vector<vector<double> > Centroid;
    vector<string> Bname;
    vector<string> name_base;
    vector<int> n_base;
    vector<vector<int> > I2E;
    vector<vector<int> > B2E;
    vector<vector<double> > In;
    vector<vector<double> > Bn;
    vector<vector<int> > NeighbourElements;
    vector<bool> isCurved;
    vector<int> CurvedElementIndex;
    vector<int> LinearElementIndex;
    vector<int> CurvedEdgeIndex;
    vector<int> LinearEdgeIndex;

    TriMesh(string &gri_filename_in);
    TriMesh(TriMesh &mesh);
    void ReadGri(string &gri_filename);
    void CalcI2E();
    void CalcB2E();
    void CalcIn();
    void CalcBn();
    void CalcArea();
    void CalcCentroid();
    void FindNeighbour();
    void FindCurvedIndex();
    void WriteGri(string &gri_filename);

};

#endif