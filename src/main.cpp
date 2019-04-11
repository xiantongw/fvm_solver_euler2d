#include <iostream>
#include "../include/TriMesh.h"
#include "../include/CalcNumericalFlux.h"
#include "../include/CalcGradient.h"
#include "../include/Collective.h"
#include "../include/MainSolver.h"

using namespace std;

int main(int argc, char* argv[])
{
    Param param;
    // Set up the param struct
    param = ReadParamIn(string(argv[1]));
    TriMesh mesh(param.mesh_file);

    /*Initialize state vectors*/
    vector<vector<double> > state_vectors(mesh.E.size(), vector<double> (4, 0.0));

    vector<double> u_free(4, 0.0);
    u_free[0] = 1.0;
    u_free[1] = param.mach_inf * cos(param.attack_angle);
    u_free[2] = param.mach_inf * sin(param.attack_angle);
    u_free[3] = 1 / (param.gamma * (param.gamma - 1)) +
                0.5 * param.mach_inf * param.mach_inf;

    for (int ielem = 0; ielem < mesh.E.size(); ielem++)
    {
        for (int istate = 0; istate < 4; istate++)
        {
            state_vectors[ielem][istate] = u_free[istate];
        }
    }

    int converged;
    converged = MainSolver(mesh, param, state_vectors);

    return 0;

}