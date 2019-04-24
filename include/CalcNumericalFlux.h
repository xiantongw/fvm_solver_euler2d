//
// Created by Xiantong Wang on 2019-04-09.
//


#ifndef CALCNUMERICALFLUX_H
#define CALCNUMERICALFLUX_H

#include <iostream>
#include <cmath>
#include <vector>
#include <string.h>

// #define strncasecmp _strnicmp
// #define strcasecmp _stricmp

using namespace std;

vector<double> projection_2d(vector<vector<double> > F, vector<double> n);
vector<vector<double> > flux_function_2d(vector<double> u, double gamma);
vector<double> CalcNumericalFlux(vector<double> uL, vector<double> uR, vector<double> n, char const *type_flux, double gamma, double& mws);

#endif //CALCNUMERICALFLUX_H
