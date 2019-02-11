void calculate_gradient(MeshParam* mparam,
                        int E[][3], double V[][2], int I2E[][4], int B2E[][3], double In[][2], double Bn[][2], double Area[],
                        double state_vectors[][4], double gradu[][4][2]){
    
    int i, j, iedge, ielem, index_left, index_right, index_pA, index_pB;
    int index_edge_nodes[2];
    double dl;
    double norm_vector[2], average_state[4];

    /* Initialize the gradu in each calculation*/
    for (ielem = 0; ielem < mparam->nelem; ielem++){
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[ielem][i][j] = 0.0;
            }
        }
    }
    
    /* Calculate the gradient for all elements*/
    // Loop through interior edeges
    for (iedge = 0; iedge < mparam->niedge; iedge++){

        norm_vector[0] = In[iedge][0]; 
        norm_vector[1] = In[iedge][1];
        index_left = I2E[iedge][0] - 1; 
        index_right = I2E[iedge][2] - 1;

        /* Calculate the length of the edge*/ 
        mod3(I2E[iedge][1] - 1, index_edge_nodes);
        
        index_pA = E[index_left][index_edge_nodes[0]] - 1;
        index_pB = E[index_left][index_edge_nodes[1]] - 1;

        dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2.0) + pow(V[index_pA][1] - V[index_pB][1], 2.0));

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
    for (iedge = 0; iedge < mparam->nbedge; iedge++){

        norm_vector[0] = Bn[iedge][0]; 
        norm_vector[1] = Bn[iedge][1];

        index_left = B2E[iedge][0] - 1;
        /* Calculate the length of the edge*/ 
        mod3(B2E[iedge][1] - 1, index_edge_nodes);

        index_pA = E[index_left][index_edge_nodes[0]] - 1;
        index_pB = E[index_left][index_edge_nodes[1]] - 1;

        dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2.0) + pow(V[index_pA][1] - V[index_pB][1], 2.0));

        // Apply the average state on the gradu calculation
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[index_left][i][j] += state_vectors[index_left][i] * norm_vector[j] * dl;
            }
        }

    }

    for (ielem = 0; ielem < mparam->nelem; ielem++){
        for (i = 0; i < 4; i++){
            for (j = 0; j < 2; j++){
                gradu[ielem][i][j] = gradu[ielem][i][j] / Area[ielem];
            }
        }
    }

}