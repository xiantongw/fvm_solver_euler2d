//
// Created by Xiantong Wang on 2019-04-10.
//

#ifndef CALCRESIDUAL_H
#define CALCRESIDUAL_H

#include <iostream>
#include <vector>

#include "utils.h"
#include "ApplyBoundaryCondition.h"
#include "CalcNumericalFlux.h"
#include "Param.h"
#include "TriMesh.h"

vector<vector<double> > CalcResidual(TriMesh mesh, Param& param, vector<vector<vector<double> > > gradu,
                                     vector<vector<double> > state_vectors, vector<double>& dtA);

#endif //CALCRESIDUAL_H
