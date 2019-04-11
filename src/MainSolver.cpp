//
// Created by Xiantong Wang on 2019-04-10.
//

#include "../include/MainSolver.h"

using namespace std;

int MainSolver(TriMesh mesh, Param param, vector<vector<double> > state_vectors)
{
    vector<vector<double> > residual(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<double> > residual_FE(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<vector<double> > > gradu(mesh.E.size(), vector<vector<double> > (4, vector<double> (2, 0.0)));
    vector<vector<double> > state_vectors_FE(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<vector<double> > > gradu_FE(mesh.E.size(), vector<vector<double> > (4, vector<double> (2, 0.0)));
    vector<double> dtA (mesh.E.size(), 0.0);
    vector<double> dtA_temp (mesh.E.size(), 0.0);
    int nstate = 4;

    /* Main Iteration*/
    if (param.order == 1)
    {
        for (int niter = 0; niter < param.MAXITER; niter++)
        {
            /* Calculate residual for the first time*/
            residual = CalcResidual(mesh, param, gradu, state_vectors, dtA);
            /*Time Marching FORWARD EULER*/
            for (int ielem = 0; ielem < mesh.E.size(); ielem++)
            {
                for (int istate = 0; istate < nstate; istate++)
                {
                    state_vectors[ielem][istate] -= dtA[ielem] * residual[ielem][istate];
                }
            }
            std::cout << "NITER: " << niter << "\t" << "Residual Norm_Inf: ";
            cout.setf(ios::scientific, ios::floatfield);
            std::cout << setprecision(10) << utils::norm_inf(residual) << std::endl;
            if (utils::norm_inf(residual) < param.eps)
            {
                return 0;
            }
        }
    } else if (param.order == 2)
    {
        for (int niter = 0; niter < param.MAXITER; niter++)
        {
            /*Calculate Gradient for the first time*/
            gradu = CalcGradient(mesh, state_vectors);
            residual = CalcResidual(mesh, param, gradu, state_vectors, dtA);
            for (int ielem = 0; ielem < mesh.E.size(); ielem++)
            {
                for (int istate = 0; istate < nstate; istate++)
                {
                    state_vectors_FE[ielem][istate] = state_vectors[ielem][istate] - dtA[ielem] * residual[ielem][istate];
                }
            }
            gradu_FE = CalcGradient(mesh, state_vectors_FE);
            residual_FE = CalcResidual(mesh, param, gradu, state_vectors_FE, dtA_temp);
            for (int ielem = 0; ielem < mesh.E.size(); ielem++)
            {
                for (int istate = 0; istate < nstate; istate++)
                {
                    state_vectors[ielem][istate] = 0.5 * (state_vectors[ielem][istate] + state_vectors_FE[ielem][istate]- dtA[ielem] * residual_FE[ielem][istate]);
                }
            }
            std::cout << "NITER: " << niter << "\t" << "Residual Norm_Inf: ";
            cout.setf(ios::scientific, ios::floatfield);
            std::cout << setprecision(10) << utils::norm_inf(residual) << std::endl;
            if (utils::norm_inf(residual) < param.eps)
            {
                return 0;
            }
        }
    }
    return -1;
}