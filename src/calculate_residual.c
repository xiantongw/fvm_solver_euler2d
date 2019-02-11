void calculate_residual(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam, 
                        int nstage, int E[][3], double V[][2], int I2E[][4], int B2E[][3],
                        double In[][2], double Bn[][2], double Area[], double Centroid[][2], double gradu[][4][2],
                        double state_vectors[][4], double residual[][4], double max_wave_speed_tally[]){

    /* Local Vars*/ 
    int niter, iedge, i, ielem, index_left, index_right, index_pA, index_pB;
    double norm_vector[2];
    double state_left[4], state_right[4];
    int index_edge_nodes[2];
    double dl;
    double mws;
    double num_flux[4];
    double mid_point[2];

    /* Initialize the residual and tally to 0.0*/
    for (ielem = 0; ielem < mparam->nelem; ielem++){
        max_wave_speed_tally[ielem] = 0.0;
        for (i = 0; i < 4; i++){
            residual[ielem][i] = 0.0;
        }
    }

    /* Loop through interior edges*/
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

        /* Collect numerical flux information*/
        
        if (nstage == 1){
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
                state_right[i] = state_vectors[index_right][i];
            }
        }
        else if (nstage == 2){

            // edge middle point
            mid_point[0] = 0.5 * (V[index_pA][0] + V[index_pB][0]);
            mid_point[1] = 0.5 * (V[index_pA][1] + V[index_pB][1]);

            // get left state
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i] + 
                                (gradu[index_left][i][0] * (mid_point[0] - Centroid[index_left][0]) + 
                                 gradu[index_left][i][1] * (mid_point[1] - Centroid[index_left][1]));
            }

            // get right state
            for (i = 0; i < 4; i++){
                state_right[i] = state_vectors[index_right][i] + 
                                (gradu[index_right][i][0] * (mid_point[0] - Centroid[index_right][0]) + 
                                 gradu[index_right][i][1] * (mid_point[1] - Centroid[index_right][1]));
            }
            
        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }
        
        roe_2d(state_left, state_right, norm_vector, cparam->gamma, num_flux, &mws);

        for (i = 0; i < 4; i++){
            residual[index_left][i] += num_flux[i] * dl;
            residual[index_right][i] -= num_flux[i] * dl;
        }
        max_wave_speed_tally[index_left] += mws * dl;
        max_wave_speed_tally[index_right] += mws * dl;
    }

    /* Loop through boundary edges*/
    for (iedge = 0; iedge < mparam->nbedge; iedge++){
        
        norm_vector[0] = Bn[iedge][0]; 
        norm_vector[1] = Bn[iedge][1];
        index_left = B2E[iedge][0] - 1;

        /* Calculate the length of the edge*/ 
        mod3(B2E[iedge][1] - 1, index_edge_nodes);

        index_pA = E[index_left][index_edge_nodes[0]] - 1;
        index_pB = E[index_left][index_edge_nodes[1]] - 1;
        dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2.0) + pow(V[index_pA][1] - V[index_pB][1], 2.0));
        
        /* Collect numerical flux information*/
        if (nstage == 1){
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
            }

        } else if (nstage == 2){

            // edge middle point
            mid_point[0] = 0.5 * (V[index_pA][0] + V[index_pB][0]);
            mid_point[1] = 0.5 * (V[index_pA][1] + V[index_pB][1]);

            // get left state
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i] + \
                                (gradu[index_left][i][0] * (mid_point[0] - Centroid[index_left][0]) + \
                                 gradu[index_left][i][1] * (mid_point[1] - Centroid[index_left][1]));   
            }
        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }

        /* Get the state vector outside*/
        switch (B2E[iedge][2] - 1)
        { 
            case 0:
                boudary_condition(state_left, norm_vector, bparam->bound0, cparam, num_flux, &mws);
                break;
            case 1:
                boudary_condition(state_left, norm_vector, bparam->bound1, cparam, num_flux, &mws);
                break;
            case 2:
                boudary_condition(state_left, norm_vector, bparam->bound2, cparam, num_flux, &mws);
                break;
            case 3:
                boudary_condition(state_left, norm_vector, bparam->bound3, cparam, num_flux, &mws);
                break;
        }
        for (i = 0; i < 4; i++){
            residual[index_left][i] += num_flux[i] * dl;
        }
        max_wave_speed_tally[index_left] += mws * dl;
    }
}
