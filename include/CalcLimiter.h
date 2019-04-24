//
// Created by Xiantong Wang on 2019-04-17.
//

#ifndef CALCLIMITER_H
#define CALCLIMITER_H

// #define strncasecmp _strnicmp
// #define strcasecmp _stricmp

#include <iostream>
#include <vector>
#include <algorithm>
#include <string.h>

#include "utils.h"
#include "Param.h"
#include "TriMesh.h"

using namespace std;

vector<vector<double> > CalcLimiter(TriMesh mesh, Param param, vector<vector<double> > state_vectors, vector<vector<vector<double> > > gradu);

vector<double> LP_Optimizer(vector<double> c, vector<vector<double> > A, vector<double> b, vector<int> W0, vector<double> x0);

#endif //CALCLIMITER_H
