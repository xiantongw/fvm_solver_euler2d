//
// Created by Xiantong Wang on 2019-04-10.
//

#ifndef PARAM_H
#define PARAM_H

#include <iostream>
#include <string>

// Simulation Parameter is stored in a Struct
typedef struct Param{
    double cfl;
    double mach_inf;
    double attack_angle;
    double gamma;
    double p_inf;
    double T_inf;
    double R;
    double eps;
    double h;
    int dnOutput;
    int dnSaveRestart;
    int MAXITER;
    std::string bound0;
    std::string bound1;
    std::string bound2;
    std::string bound3;
    std::string mesh_file;
    std::string name_limiter;
} Param;

#endif

