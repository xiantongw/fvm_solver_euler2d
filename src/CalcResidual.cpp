//
// Created by Xiantong Wang on 2019-04-10.
//

#include "../include/CalcResidual.h"

using namespace std;

vector<vector<double> > CalcResidual(TriMesh mesh, Param& param, vector<vector<vector<double> > > gradu,
                        vector<vector<double> > state_vectors, vector<double>& dtA){

    /* Local Vars*/
    int iedge, i, ielem, index_left, index_right, index_pA, index_pB;
    vector<double> norm_vector(2, 0.0);
    vector<double> state_left(4, 0.0), state_right(4, 0.0);
    int index_edge_nodes[2];
    double dl;
    double mws;
    vector<double> num_flux(4, 0.0);
    double mid_point[2];
    vector<vector<double> > residual(mesh.E.size(), vector<double>(4, 0.0));
    vector<double> max_wave_speed_tally(mesh.E.size(), 0.0);

    /* Initialize the residual and tally to 0.0*/
    for (ielem = 0; ielem < mesh.E.size(); ielem++){
        max_wave_speed_tally[ielem] = 0.0;
        for (i = 0; i < 4; i++){
            residual[ielem][i] = 0.0;
        }
    }

    /* Loop through interior edges*/
    for (iedge = 0; iedge < mesh.I2E.size(); iedge++){

        norm_vector[0] = mesh.In[iedge][0];
        norm_vector[1] = mesh.In[iedge][1];
        index_left = mesh.I2E[iedge][0] - 1;
        index_right = mesh.I2E[iedge][2] - 1;

        dl = mesh.In[iedge][2];

        /* Get the local index of the edge points*/
        utils::mod3(mesh.I2E[iedge][1] - 1, index_edge_nodes);

        index_pA = mesh.E[index_left][index_edge_nodes[0]] - 1;
        index_pB = mesh.E[index_left][index_edge_nodes[1]] - 1;

        /* Collect numerical flux information*/

        if (param.order == 1){
            state_left = state_vectors[index_left];
            state_right = state_vectors[index_right];
        }
        else if (param.order == 2){

            // edge middle point
            mid_point[0] = 0.5 * (mesh.V[index_pA][0] + mesh.V[index_pB][0]);
            mid_point[1] = 0.5 * (mesh.V[index_pA][1] + mesh.V[index_pB][1]);

            // get left state
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i] +
                                (gradu[index_left][i][0] * (mid_point[0] - mesh.Centroid[index_left][0]) +
                                 gradu[index_left][i][1] * (mid_point[1] - mesh.Centroid[index_left][1]));
            }

            // get right state
            for (i = 0; i < 4; i++){
                state_right[i] = state_vectors[index_right][i] +
                                 (gradu[index_right][i][0] * (mid_point[0] - mesh.Centroid[index_right][0]) +
                                  gradu[index_right][i][1] * (mid_point[1] - mesh.Centroid[index_right][1]));
            }

        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }

        num_flux = CalcNumericalFlux(state_left, state_right, norm_vector, "roe", param.gamma, mws);

        for (i = 0; i < 4; i++){
            residual[index_left][i] += num_flux[i] * dl;
            residual[index_right][i] -= num_flux[i] * dl;
        }
        max_wave_speed_tally[index_left] += mws * dl;
        max_wave_speed_tally[index_right] += mws * dl;
    }

    /* Loop through boundary edges*/
    for (iedge = 0; iedge < mesh.B2E.size(); iedge++){

        norm_vector[0] = mesh.Bn[iedge][0];
        norm_vector[1] = mesh.Bn[iedge][1];
        index_left = mesh.B2E[iedge][0] - 1;

        dl = mesh.Bn[iedge][2];

        /* Get the local index of the edge points*/
        utils::mod3(mesh.B2E[iedge][1] - 1, index_edge_nodes);

        index_pA = mesh.E[index_left][index_edge_nodes[0]] - 1;
        index_pB = mesh.E[index_left][index_edge_nodes[1]] - 1;

        /* Collect numerical flux information*/
        if (param.order == 1) {
            state_left = state_vectors[index_left];
        }else if (param.order == 2){

            // edge middle point
            mid_point[0] = 0.5 * (mesh.V[index_pA][0] + mesh.V[index_pB][0]);
            mid_point[1] = 0.5 * (mesh.V[index_pA][1] + mesh.V[index_pB][1]);

            // get left state
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i] + \
                                (gradu[index_left][i][0] * (mid_point[0] - mesh.Centroid[index_left][0]) + \
                                 gradu[index_left][i][1] * (mid_point[1] - mesh.Centroid[index_left][1]));
            }
        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }

        /* Get the state vector outside*/
        switch (mesh.B2E[iedge][2] - 1)
        {
            case 0:
                num_flux = ApplyBoundaryCondition(state_left, norm_vector, param.bound0, param, mws);
                break;
            case 1:
                num_flux = ApplyBoundaryCondition(state_left, norm_vector, param.bound1, param, mws);
                break;
            case 2:
                num_flux = ApplyBoundaryCondition(state_left, norm_vector, param.bound2, param, mws);
                break;
            case 3:
                num_flux = ApplyBoundaryCondition(state_left, norm_vector, param.bound3, param, mws);
                break;
        }
        for (i = 0; i < 4; i++){
            residual[index_left][i] += num_flux[i] * dl;
        }
        max_wave_speed_tally[index_left] += mws * dl;
    }
    // Calculate dt/A
    for (ielem = 0; ielem < mesh.E.size(); ielem++)
    {
        dtA[ielem] = 2.0 * param.cfl / max_wave_speed_tally[ielem];
    }

    return residual;

}
