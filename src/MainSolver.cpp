//
// Created by Xiantong Wang on 2019-04-10.
//

#include "../include/MainSolver.h"

using namespace std;

int MainSolver(TriMesh mesh, Param param, vector<vector<double> >& state_vectors)
{
    vector<vector<double> > residual(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<double> > residual_FE(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<vector<double> > > gradu(mesh.E.size(), vector<vector<double> > (4, vector<double> (2, 0.0)));
    vector<vector<double> > state_vectors_FE(mesh.E.size(), vector<double>(4, 0.0));
    vector<vector<vector<double> > > gradu_FE(mesh.E.size(), vector<vector<double> > (4, vector<double> (2, 0.0)));
    vector<double> dtA (mesh.E.size(), 0.0);
    vector<double> dtA_temp (mesh.E.size(), 0.0);
    vector<vector<double> > limiter;

    int nstate = 4; int num_elem = mesh.E.size();

    /* Main Iteration*/
    ofstream file_residual;
    file_residual.open("residual.log");
    for (int niter = 0; niter < param.MAXITER; niter++)
    {
        /*Calculate Gradient for the first time*/
        gradu = CalcGradient(mesh, state_vectors);
        limiter = CalcLimiter(mesh, param, state_vectors, gradu);
        residual = CalcResidual(mesh, param, gradu, state_vectors, limiter, dtA);
        for (int ielem = 0; ielem < num_elem; ielem++)
        {
            for (int istate = 0; istate < nstate; istate++)
            {
                state_vectors_FE[ielem][istate] =
                        state_vectors[ielem][istate] - dtA[ielem] * residual[ielem][istate];
            }
        }
        gradu_FE = CalcGradient(mesh, state_vectors_FE);
        limiter = CalcLimiter(mesh, param, state_vectors_FE, gradu_FE);
        residual_FE = CalcResidual(mesh, param, gradu, state_vectors_FE, limiter, dtA_temp);
        for (int ielem = 0; ielem < num_elem; ielem++)
        {
            for (int istate = 0; istate < nstate; istate++)
            {
                state_vectors[ielem][istate] = 0.5 *
                                               (state_vectors[ielem][istate] + state_vectors_FE[ielem][istate] -
                                                dtA[ielem] * residual_FE[ielem][istate]);
            }
        }
        std::cout << "NITER: " << niter << "\t" << "Residual Norm_Inf: ";
        cout.setf(ios::scientific, ios::floatfield);
        std::cout << setprecision(10) << utils::norm_inf(residual) << std::endl;
	    if (niter % param.dnSaveRestart == 0)
        {
            utils::SaveState(state_vectors);
        }
        file_residual << niter << "\t" << setprecision(20) << utils::norm_inf(residual) << std::endl;
        if (utils::norm_inf(residual) < param.eps)
        {
            file_residual.close();
            return 0;
        }
    }
    file_residual.close();
    return -1;
}
