//
// Created by Xiantong Wang on 2019-04-10.
//

#include "../include/Collective.h"

using namespace std;

Param ReadParamIn(string param_filename)
{
    // variables dealing with strings
    ifstream param_file(param_filename);
    stringstream ss;
    string line, param_name, param_value;

    Param param;
    while (getline(param_file, line))
    {
        ss.clear();
        ss.str(line);
        ss >> param_name >> param_value;
        if (strcasecmp(param_name.c_str(), "gamma") == 0)
        {
            param.gamma = atof(param_value.c_str());
        } else if (strcasecmp(param_name.c_str(), "attack_angle") == 0)
        {
            param.attack_angle = atof(param_value.c_str());
        } else if (strcasecmp(param_name.c_str(), "cfl") == 0)
        {
            param.cfl = atof(param_value.c_str());
        }else if (strcasecmp(param_name.c_str(), "mach_inf") == 0)
        {
            param.mach_inf = atof(param_value.c_str());
        }else if (strcasecmp(param_name.c_str(), "p_inf") == 0)
        {
            param.p_inf = atof(param_value.c_str());
        }else if (strcasecmp(param_name.c_str(), "bound0") == 0)
        {
            param.bound0 = param_value;
        }else if (strcasecmp(param_name.c_str(), "bound1") == 0)
        {
            param.bound1 = param_value;
        }else if (strcasecmp(param_name.c_str(), "bound2") == 0)
        {
            param.bound2 = param_value;
        }else if (strcasecmp(param_name.c_str(), "bound3") == 0)
        {
            param.bound3 = param_value;
        }else if (strcasecmp(param_name.c_str(), "mesh_file") == 0)
        {
            param.mesh_file = param_value;
        }else if (strcasecmp(param_name.c_str(), "order") == 0)
        {
            param.order = int(atof(param_value.c_str()));
        }else if (strcasecmp(param_name.c_str(), "eps") == 0)
        {
            param.eps = atof(param_value.c_str());
        }else if (strcasecmp(param_name.c_str(), "MAXITER") == 0)
        {
            param.MAXITER = int(atof(param_value.c_str()));
        }else if (strcasecmp(param_name.c_str(), "T_inf") == 0)
        {
            param.T_inf = atof(param_value.c_str());
        }
        else if (strcasecmp(param_name.c_str(), "R") == 0)
        {
            param.R = atof(param_value.c_str());
        }else if (strcasecmp(param_name.c_str(), "dnOutput") == 0)
        {
            param.dnOutput = int(atof(param_value.c_str()));
        }
    }
    param_file.close();
    return param;
}
