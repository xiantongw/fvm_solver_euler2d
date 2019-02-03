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

int euler_solver_first_order(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam,
                        int E[][3], double V[][2], 
                        int I2E[][4], int B2E[][3],
                        double In[][2], double Bn[][2], double Area[],
                        double state_vectors[][4],
                        int MAXITER, double eps){
    
    /* Local Vars*/
    double residual[mparam->nelem][4];
    double max_wave_speed_tally[mparam->nelem];
    int max_wave_speed_tally_counter[mparam->nelem];

    int niter, iedge;
    double norm_vector[2];
    int index_left, index_right;
    double state_left[4], state_right[4];
    int index_edge_nodes[2];
    double dl;
    double mws;
    double num_flux[4];
    
    /* Main Iteration*/
    for (niter=0; niter<MAXITER; niter++){

        /* Set temp arrays into 0 in every iteration*/
        for (int ielem=0; ielem < mparam->nelem; ielem++){
            residual[ielem][0] = 0.0; residual[ielem][1] = 0.0; residual[ielem][2] = 0.0; residual[ielem][3] = 0.0;
            max_wave_speed_tally[ielem] = 0.0;
            max_wave_speed_tally_counter[ielem] = 0;
        }

        /* Loop through interior edges*/
        for (iedge = 0; iedge < mparam->niedge; iedge++){
            norm_vector[0] = In[iedge][0]; norm_vector[1] = In[iedge][1];
            index_left = I2E[iedge][0] - 1; index_right = I2E[iedge][2] - 1;

            /* Calculate the length of the edge*/ 
            mod3(I2E[iedge][1] - 1, index_edge_nodes);
            
            int index_pA = E[index_left][index_edge_nodes[0]] - 1;
            int index_pB = E[index_left][index_edge_nodes[1]] - 1;

            dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + \
                        pow(V[index_pA][1] - V[index_pB][1], 2));

            /* Collect numerical flux information*/
            for (int i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
                state_right[i] = state_vectors[index_right][i];
            }
            
            roe_2d(state_left, state_right, norm_vector, cparam->gamma, num_flux, &mws);

            for (int i = 0; i < 4; i++){
                residual[index_left][i] += num_flux[i] * dl;
                residual[index_right][i] -= num_flux[i] * dl;
            }

            max_wave_speed_tally[index_left] += mws * dl;
            max_wave_speed_tally[index_right] += mws * dl;
            max_wave_speed_tally_counter[index_left] += 1;
            max_wave_speed_tally_counter[index_right] += 1;
        }

        /* Loop through boundary edges*/
        for (iedge = 0; iedge < mparam->nbedge; iedge++){
            
            norm_vector[0] = Bn[iedge][0]; norm_vector[1] = Bn[iedge][1];
            index_left = B2E[iedge][0] - 1;

            /* Calculate the length of the edge*/ 
            mod3(B2E[iedge][1] - 1, index_edge_nodes);

            int index_pA = E[index_left][index_edge_nodes[0]] - 1;
            int index_pB = E[index_left][index_edge_nodes[1]] - 1;
            dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + \
                        pow(V[index_pA][1] - V[index_pB][1], 2));
            
            /* Collect numerical flux information*/
            for (int i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
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
            for (int i = 0; i < 4; i++){
                residual[index_left][i] += num_flux[i] * dl;
            }
            max_wave_speed_tally[index_left] += mws * dl;
            max_wave_speed_tally_counter[index_left] += 1;
            
        }

        /* Time Marching*/
        int converge_flag = 1;
        int ielem, icolumn;
        double row[4];
        for (ielem = 0; ielem < mparam->nelem; ielem++){
            for (icolumn = 0; icolumn < 4; icolumn++){
                row[icolumn] = residual[ielem][icolumn];
            }
            if (norm1d(row, 4) > eps){
                converge_flag = 0;
                // Time marching on this local element
                for (icolumn = 0; icolumn < 4; icolumn++){
                    state_vectors[ielem][icolumn] = state_vectors[ielem][icolumn] - \
                        (2.0 * cparam->cfl * residual[ielem][icolumn]) / max_wave_speed_tally[ielem];
                }
            }
        }

        // Print the simulation log
        printf("niter: %6d \t \t residual_norm: %9.8e \n", niter+1, norm2d(residual, mparam->nelem));

        if (converge_flag){
            printf("Converged \n");
            return 1;
        }
    }
    return 0;
}


int euler_solver_second_order(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam,
                        int E[][3], double V[][2], 
                        int I2E[][4], int B2E[][3],
                        double In[][2], double Bn[][2], double Area[],
                        double state_vectors[][4],
                        int MAXITER, double eps){

    printf("%s\n", "second order");
    
    /* Local Vars*/
    double residual[mparam->nelem][4];
    double max_wave_speed_tally[mparam->nelem];
    int max_wave_speed_tally_counter[mparam->nelem];

    int niter, iedge;
    double norm_vector[2];
    int index_left, index_right;
    double state_left[4], state_right[4];
    int index_edge_nodes[2];
    double dl;
    double mws;
    double num_flux[4];
    
    /* Main Iteration*/
    for (niter=0; niter<MAXITER; niter++){

        /* Set temp arrays into 0 in every iteration*/
        for (int ielem=0; ielem < mparam->nelem; ielem++){
            residual[ielem][0] = 0.0; residual[ielem][1] = 0.0; residual[ielem][2] = 0.0; residual[ielem][3] = 0.0;
            max_wave_speed_tally[ielem] = 0.0;
            max_wave_speed_tally_counter[ielem] = 0;
        }

        /* Loop through interior edges*/
        for (iedge = 0; iedge < mparam->niedge; iedge++){
            norm_vector[0] = In[iedge][0]; norm_vector[1] = In[iedge][1];
            index_left = I2E[iedge][0] - 1; index_right = I2E[iedge][2] - 1;

            /* Calculate the length of the edge*/ 
            mod3(I2E[iedge][1] - 1, index_edge_nodes);
            
            int index_pA = E[index_left][index_edge_nodes[0]] - 1;
            int index_pB = E[index_left][index_edge_nodes[1]] - 1;

            dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + \
                        pow(V[index_pA][1] - V[index_pB][1], 2));

            /* Collect numerical flux information*/
            for (int i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
                state_right[i] = state_vectors[index_right][i];
            }
            
            roe_2d(state_left, state_right, norm_vector, cparam->gamma, num_flux, &mws);

            for (int i = 0; i < 4; i++){
                residual[index_left][i] += num_flux[i] * dl;
                residual[index_right][i] -= num_flux[i] * dl;
            }

            max_wave_speed_tally[index_left] += mws * dl;
            max_wave_speed_tally[index_right] += mws * dl;
            max_wave_speed_tally_counter[index_left] += 1;
            max_wave_speed_tally_counter[index_right] += 1;
        }

        /* Loop through boundary edges*/
        for (iedge = 0; iedge < mparam->nbedge; iedge++){
            
            norm_vector[0] = Bn[iedge][0]; norm_vector[1] = Bn[iedge][1];
            index_left = B2E[iedge][0] - 1;

            /* Calculate the length of the edge*/ 
            mod3(B2E[iedge][1] - 1, index_edge_nodes);

            int index_pA = E[index_left][index_edge_nodes[0]] - 1;
            int index_pB = E[index_left][index_edge_nodes[1]] - 1;
            dl = sqrt(pow(V[index_pA][0] - V[index_pB][0], 2) + \
                        pow(V[index_pA][1] - V[index_pB][1], 2));
            
            /* Collect numerical flux information*/
            for (int i = 0; i < 4; i++){
                state_left[i] = state_vectors[index_left][i];
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
            for (int i = 0; i < 4; i++){
                residual[index_left][i] += num_flux[i] * dl;
            }
            max_wave_speed_tally[index_left] += mws * dl;
            max_wave_speed_tally_counter[index_left] += 1;
            
        }

        /* Time Marching*/
        int converge_flag = 1;
        int ielem, icolumn;
        double row[4];
        for (ielem = 0; ielem < mparam->nelem; ielem++){
            for (icolumn = 0; icolumn < 4; icolumn++){
                row[icolumn] = residual[ielem][icolumn];
            }
            if (norm1d(row, 4) > eps){
                converge_flag = 0;
                // Time marching on this local element
                for (icolumn = 0; icolumn < 4; icolumn++){
                    state_vectors[ielem][icolumn] = state_vectors[ielem][icolumn] - \
                        (2.0 * cparam->cfl * residual[ielem][icolumn]) / max_wave_speed_tally[ielem];
                }
            }
        }

        // Print the simulation log
        printf("niter: %6d \t \t residual_norm: %9.8e \n", niter+1, norm2d(residual, mparam->nelem));

        if (converge_flag){
            printf("Converged \n");
            return 1;
        }
    }
    return 0;
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
