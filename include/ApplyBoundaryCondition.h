//
// Created by Xiantong Wang on 2019-04-10.
//

#ifndef APPLYBOUNDARYCONDITION_H
#define APPLYBOUNDARYCONDITION_H

#include <iostream>
#include <vector>
#include <cmath>
#include <cstring>
#include <strings.h>

#include "../include/Param.h"
#include "../include/CalcNumericalFlux.h"

vector<double> ApplyBoundaryCondition(vector<double> u, vector<double> norm,  string boundary_type, Param& cparam, double &mws);

#endif //APPLYBOUNDARYCONDITION_H
