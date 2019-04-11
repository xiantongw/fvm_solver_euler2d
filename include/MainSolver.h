//
// Created by Xiantong Wang on 2019-04-10.
//

#ifndef MAINSOLVER_H
#define MAINSOLVER_H

#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/CalcGradient.h"
#include "../include/CalcResidual.h"
#include "../include/TriMesh.h"
#include "../include/Param.h"

int MainSolver(TriMesh mesh, Param param, vector<vector<double> > state_vectors);

#endif //MAINSOLVER_H
