/*
    Euler Solver: Node-Centered Finite-Volume-Method
    
    Node-centered 1st and 2nd order finite-volume method for unstructured grids
    Roe flux with an entropy fix
*/

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "flux.c"
#include "params.h"
#include "boundary_condition.c"
#include "utils.c"
#include "calculate_gradient.c"
#include "calculate_residual.c"

int euler_solver_main(CfdParam* cparam, MeshParam* mparam, BoundaryParam* bparam, int norder,
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

    /* Simualtion Information*/

    printf("2D Euler Solver: Node-Centered Finite-Volume-Method\n");
    printf("Numerical Flux: 2D Roe flux with an entropy fix\n");
    printf("Order of Accuracy: %d\n", norder);
    printf("\n");
    printf("Mesh Parameters:\n");
    printf("Number of Elements: %d\n", mparam->nelem);
    printf("Number of Interior Edges: %d\n", mparam->niedge);
    printf("Number of Boundary Edges: %d\n", mparam->nbedge);
    printf("\n");




    printf("****************************************\n");
    printf("           Simulation Start             \n");
    printf("****************************************\n");

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
        if (norder == 1){
            ;
        }
        else if (norder == 2){
            calculate_gradient(mparam, E, V, I2E, B2E, In, Bn, Area, state_vectors, gradu);
        }
        else {
            printf("Unsupported order of accuracy!\n");
            abort();
        }

        calculate_residual(cparam, mparam, bparam, norder, E, V, I2E, B2E, In, Bn, Area, Centroid, gradu,
                           state_vectors, residual, max_wave_speed_tally);

        /* Time Marching*/
        int converge_flag = 1; // set convergence flag to TRUE
        if (norder == 1){
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
        else if (norder == 2){

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
            calculate_residual(cparam, mparam, bparam, norder, E, V, I2E, B2E, In, Bn, Area, Centroid, gradu_FE,
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
        printf("niter: %6d \t Res: %.18e \n", niter+1, norm2d(residual, mparam->nelem));

        if (converge_flag){
            printf("Converged \n");
            return 1;
        }

    }
    return 0;
}
