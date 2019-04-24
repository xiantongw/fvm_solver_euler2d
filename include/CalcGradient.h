//
// Created by Xiantong Wang on 2019-04-10.
//

#ifndef CALCGRADIENT_H
#define CALCGRADIENT_H

#include <iostream>
#include <vector>
#include "../include/TriMesh.h"
#include "../include/utils.h"

vector<vector<vector<double> > > CalcGradient(TriMesh mesh, vector<vector<double> > state_vectors);

#endif //CALCGRADIENT_H
