/*
    Euler Solver: Node-Centered Finite-Volume-Method
    
    Node-centered 1st and 2nd order finite-volume method for unstructured grids
    Roe flux with an entropy fix
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "boundary_condition.c"

void mod3(int input, int res[2]);
double norm1d(double r[], int n);
double norm2d(double r[][4], int n);

void calculate_gradient(MeshParam* mparam,
                        int E[][3], double V[][2], int I2E[][4], int B2E[][3], double In[][2], double Bn[][2], double Area[],
                        double state_vectors[][4], double gradu[][4][2]);

void calculate_residual(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam, 
                        int nstage, int E[][3], double V[][2], int I2E[][4], int B2E[][3],
                        double In[][2], double Bn[][2], double Area[], double Centroid[][2], double gradu[][4][2],
                        double state_vectors[][4], double residual[][4], double max_wave_speed_tally[]);


int euler_solver_main(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam, int nstage,
                        int E[][3], double V[][2], 
                        int I2E[][4], int B2E[][3],
                        double In[][2], double Bn[][2], double Area[], double Centroid[][2],
                        double state_vectors[][4],
                        int MAXITER, double eps){
    
    /* Local Vars*/
    int niter, iedge;
    double norm_vector[2];
    int index_left, index_right;
    double state_left[4], state_right[4];
    int index_edge_nodes[2];
    double dl;
    double mws;
    double num_flux[4];
    double max_wave_speed_tally[mparam->nelem];
    double gradu[mparam->nelem][4][2]; 
    double residual[mparam->nelem][4];
    double dtA[mparam->nelem];
    int i, j, ielem, icolumn;

    /* Main Iteration*/
    for (niter=0; niter<MAXITER; niter++){
        
        /* Set temp arrays into 0 in every iteration*/
        for (ielem=0; ielem < mparam->nelem; ielem++){
            max_wave_speed_tally[ielem] = 0.0;
            for (i = 0; i < 4; i++){
                residual[ielem][i] = 0.0;
                for (j = 0; j < 4; j++){
                    gradu[ielem][i][j] = 0.0;
                }
            }
        }

        /* Calculate residual*/
        if (nstage == 1){
            ;
        }
        else if (nstage == 2){
            calculate_gradient(mparam, E, V, I2E, B2E, In, Bn, Area, state_vectors, gradu);
            for (int i=0; i<mparam->nelem; i++){
                for (int j=0;j<4;j++){
                    // printf("%10.6e %10.6e ", gradu[i][j][0], gradu[i][j][1]);
                }
                // printf("\n");
            }
        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }

        calculate_residual(cparam, mparam, bparam, nstage, E, V, I2E, B2E, In, Bn, Area, Centroid, gradu,
                           state_vectors, residual, max_wave_speed_tally);
        
        /* Time Marching*/
        int converge_flag = 1; // set convergence flag to TRUE
        if (nstage == 1){
            double row[4];
            // Time step
            for (ielem = 0; ielem < mparam->nelem; ielem++){
                    dtA[ielem] = 2.0 * cparam->cfl / max_wave_speed_tally[ielem];
            }

            for (ielem = 0; ielem < mparam->nelem; ielem++){
                for (icolumn = 0; icolumn < 4; icolumn++){
                    row[icolumn] = residual[ielem][icolumn];
                }
                if (norm1d(row, 4) > eps){
                    converge_flag = 0;
                    // Time marching on this local element
                    for (icolumn = 0; icolumn < 4; icolumn++){
                        state_vectors[ielem][icolumn] -= dtA[ielem] * residual[ielem][icolumn];
                    }
                }
            }
        }
        else if (nstage == 2){

            double row[4];
            double state_vectors_FE[mparam->nelem][4];
            double residual_FE[mparam->nelem][4];
            double gradu_FE[mparam->nelem][4][2];
            double tally_temp[mparam->nelem];
            
            // Time step
            for (ielem = 0; ielem < mparam->nelem; ielem++){
                    dtA[ielem] = 2.0 * cparam->cfl / max_wave_speed_tally[ielem];
            }

            /* Set FE state to 0.0*/
            for (ielem = 0; ielem < mparam->nelem; ielem++){
                tally_temp[ielem] = 0.0;
                for (j = 0; j < 4; j++){
                    residual_FE[ielem][j] = 0.0;
                    state_vectors_FE[ielem][j] = 0.0;
                }
            }
            
            for (ielem = 0; ielem < mparam->nelem; ielem++){
                for (i = 0; i < 4; i++){
                    state_vectors_FE[ielem][i] = state_vectors[ielem][i] - dtA[ielem] * residual[ielem][i];
                }
            }

            // Calculate gradient for FE states
            calculate_gradient(mparam, E, V, I2E, B2E, In, Bn, Area, state_vectors_FE, gradu_FE);

            // Re-calculate the residual
            calculate_residual(cparam, mparam, bparam, nstage, E, V, I2E, B2E, In, Bn, Area, Centroid, gradu_FE,
                                state_vectors_FE, residual_FE, tally_temp);

            for (ielem = 0; ielem < mparam->nelem; ielem++){
                for (j = 0; j < 4; j++){
                    row[j] = residual[ielem][j];
                }
                if (norm1d(row, 4) > eps){
                    converge_flag = 0;
                    // Time marching on this local element          
                    for (i = 0; i < 4; i++){
                        state_vectors[ielem][i] = 0.5 * (state_vectors_FE[ielem][i] + state_vectors[ielem][i] - \
                                                         dtA[ielem] * residual_FE[ielem][i]);
                    }
                }
            }
        }
        else {
            printf("Unsupported order of accuracy!\n");
            converge_flag = 0;
            abort();
        }

        // Print the simulation log
        if ((niter + 1) % 10 == 0 || niter == 0){
            printf("niter: %6d \t Res: %.18e \n", niter+1, norm2d(residual, mparam->nelem));
        }

        if (converge_flag){
            printf("Converged \n");
            return 1;
        }

    }
    return 0;
}


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

        dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + pow(V[index_pA][1] - V[index_pB][1], 2));

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

            /* Variables for second order calculation*/
            double mid_point[2];

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
        dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + \
                    pow(V[index_pA][1] - V[index_pB][1], 2));
        
        /* Collect numerical flux information*/
        if (nstage == 1){
            for (i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
            }

        } else if (nstage == 2){
            
            /* Variables for second order calculation*/
            double mid_point[2];

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

void mod3(int input, int res[2]){

    switch (input)
    {
        case 0:
            res[0] = 1;
            res[1] = 2;
            break;
        case 1:
            res[0] = 2;
            res[1] = 0;
            break;
        case 2:
            res[0] = 0;
            res[1] = 1;
            break;
        default:
            break;
    }

}

double norm1d(double r[], int n){
    double sum = 0.0;
    int i;
    for (i = 0; i < n; i++){
        sum = sum + r[i] * r[i];
    }
    return sqrt(sum);
}

double norm2d(double r[][4], int n){
    double sum = 0.0;
    int i, j;
    double row[4];
    for (i = 0; i < n; i++){
        for (j = 0; j < 4; j++){
            row[j] = r[i][j];
        } 
        sum = sum + norm1d(row, 4)*norm1d(row, 4);
    }
    return sqrt(sum);
}


