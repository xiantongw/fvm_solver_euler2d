//
// Created by Xiantong Wang on 2019-04-10.
//

#include "../include/CalcGradient.h"

using namespace std;

vector<vector<vector<double> > > CalcGradient(TriMesh mesh, vector<vector<double> > state_vectors){

    int i, j, iedge, ielem, index_left, index_right, index_pA, index_pB;
    int index_edge_nodes[2];
    double dl;
    double norm_vector[2], average_state[4];
    vector<vector<vector<double> > > gradu(mesh.E.size(), vector<vector<double> > (4, vector<double> (2, 0.0)));

    /* Initialize the gradu in each calculation*/
    for (ielem = 0; ielem < mesh.E.size(); ielem++){
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[ielem][i][j] = 0.0;
            }
        }
    }

    /* Calculate the gradient for all elements*/
    // Loop through interior edges
    for (iedge = 0; iedge < mesh.I2E.size(); iedge++){

        norm_vector[0] = mesh.In[iedge][0];
        norm_vector[1] = mesh.In[iedge][1];
        index_left = mesh.I2E[iedge][0] - 1;
        index_right = mesh.I2E[iedge][2] - 1;

        dl = mesh.In[iedge][2];
        /* Calculate the length of the edge*/
        utils::mod3(mesh.I2E[iedge][1] - 1, index_edge_nodes);

        index_pA = mesh.E[index_left][index_edge_nodes[0]] - 1;
        index_pB = mesh.E[index_left][index_edge_nodes[1]] - 1;

        // The average state
        for (i = 0; i < 4; i++){
            average_state[i] = 0.5 * (state_vectors[index_left][i] + state_vectors[index_right][i]);
        }


        // Apply the average state on the gradu calculation
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[index_left][i][j] +=  average_state[i] * norm_vector[j] * dl;
                gradu[index_right][i][j] -=  average_state[i] * norm_vector[j] * dl;
            }
        }
    }

    // Loop through boundary edges, using the interior state as the average state
    for (iedge = 0; iedge < mesh.B2E.size(); iedge++){

        norm_vector[0] = mesh.Bn[iedge][0];
        norm_vector[1] = mesh.Bn[iedge][1];

        dl = mesh.Bn[iedge][2];
        index_left = mesh.B2E[iedge][0] - 1;
        /* Calculate the length of the edge*/
        utils::mod3(mesh.B2E[iedge][1] - 1, index_edge_nodes);

        index_pA = mesh.E[index_left][index_edge_nodes[0]] - 1;
        index_pB = mesh.E[index_left][index_edge_nodes[1]] - 1;

        // Apply the average state on the gradu calculation
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[index_left][i][j] += state_vectors[index_left][i] * norm_vector[j] * dl;
            }
        }

    }

    for (ielem = 0; ielem < mesh.E.size(); ielem++){
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[ielem][i][j] = gradu[ielem][i][j] / mesh.Area[ielem];
            }
        }
    }

    return gradu;

}