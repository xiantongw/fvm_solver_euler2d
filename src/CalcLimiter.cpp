//
// Created by Xiantong Wang on 2019-04-17.
//

#include "../include/CalcLimiter.h"

vector<vector<double> > CalcLimiter(TriMesh mesh, Param param, vector<vector<double> > state_vectors, vector<vector<vector<double> > > gradu)
{
    vector<vector<double> > limiter(mesh.E.size(), vector<double> (2, 1.0));

    if (strcasecmp(param.name_limiter.c_str(), "None") == 0)
    {
        for (int i = 0; i < mesh.E.size(); i++)
        {
            for (int j = 0; j < 2; j++)
            {
                limiter[i][j] = 1.0; // Do not use limiter when doing the interpolation for state
            }
        }
    }
    if (strcasecmp(param.name_limiter.c_str(), "MC") == 0)
    {
        // The LP limiter, need to solve a all-inequality optimization problem
        // minimize c * x, constraints: Ax >= b
        // first allocate the arrays for the LP
        vector<double> c(2, 0.0);
        vector<vector<double> > A(10, vector<double> (2, 0.0));
        vector<double> b(10, 0.0);
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            if (mesh.NeighbourElements[ielem][0] > 0 && mesh.NeighbourElements[ielem][1] > 0 && mesh.NeighbourElements[ielem][2] > 0)
            {
                // If this is an interior element
                // Construct the matrices for LP Optimization
                int ind1, ind2, ind3;
                ind1 = mesh.NeighbourElements[ielem][0] - 1;
                ind2 = mesh.NeighbourElements[ielem][1] - 1;
                ind3 = mesh.NeighbourElements[ielem][2] - 1;
                vector<double> u = {state_vectors[ielem][0], state_vectors[ind1][0],
                                    state_vectors[ind2][0], state_vectors[ind3][0]};
                auto result  = std::minmax_element(u.begin(), u.end());
                double u_min = *result.first;
                double u_max = *result.second;
                double u_m = u[0];
                double x1, y1, x2, y2, x3, y3, xm, ym, Dx, Dy;
                xm = mesh.Centroid[ielem][0];   ym = mesh.Centroid[ielem][1];
                x1 = mesh.Centroid[ind1][0];    y1 = mesh.Centroid[ind1][1];
                x2 = mesh.Centroid[ind2][0];    y2 = mesh.Centroid[ind2][1];
                x3 = mesh.Centroid[ind3][0];    y3 = mesh.Centroid[ind3][1];
                Dx = gradu[ielem][0][0]; Dy = gradu[ielem][0][1];
                A[0][0] = (x1 - xm) * Dx;   A[0][1] = (y1 - ym) * Dy;
                A[1][0] = (x2 - xm) * Dx;   A[1][1] = (y2 - ym) * Dy;
                A[2][0] = (x3 - xm) * Dx;   A[2][1] = (y3 - ym) * Dy;
                A[3][0] = -(x1 - xm) * Dx;  A[3][1] = -(y1 - ym) * Dy;
                A[4][0] = -(x2 - xm) * Dx;  A[4][1] = -(y2 - ym) * Dy;
                A[5][0] = -(x3 - xm) * Dx;  A[5][1] = -(y3 - ym) * Dy;
                A[6][0] = 1.0;              A[6][1] = 0.0;
                A[7][0] = 0.0;              A[7][1] = 1.0;
                A[8][0] = -1.0;             A[8][1] = 0.0;
                A[9][0] = 0.0;              A[9][1] = -1.0;

                c[0] = abs(Dx); c[1] = abs(Dy);

                b[0] = u_min - u_m; b[1] = u_min - u_m; b[2] = u_min - u_m;
                b[3] = u_m - u_max; b[4] = u_m - u_max; b[5] = u_m - u_max;
                b[6] = 0.0; b[7] = 0.0; b[8] = -1.0; b[9] = -1.0;

                vector<double> phi(2, 0.0);
                vector<int> W0 = {6, 7};
                vector<double> x0 = {0.0, 0.0};
                phi = LP_Optimizer(c, A, b, W0, x0);
                limiter[ielem][0] = phi[0];
                limiter[ielem][1] = phi[1];
            } else
            {
                // Do not apply limiter on boundary elements
                limiter[ielem][0] = 1.0; limiter[ielem][1] = 1.0;
            }
        }
    }
    if (strcasecmp(param.name_limiter.c_str(), "MINMOD") == 0)
    {
        // The LP limiter, need to solve a all-inequality optimization problem
        // minimize c * x, constraints: Ax >= b
        // first allocate the arrays for the LP
        vector<double> c(2, 0.0);
        vector<vector<double> > A(10, vector<double> (2, 0.0));
        vector<double> b(10, 0.0);
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            if (mesh.NeighbourElements[ielem][0] > 0 && mesh.NeighbourElements[ielem][1] > 0 && mesh.NeighbourElements[ielem][2] > 0)
            {
                // If this is an interior element
                // Construct the matrices for LP Optimization
                int ind1, ind2, ind3;
                ind1 = mesh.NeighbourElements[ielem][0] - 1;
                ind2 = mesh.NeighbourElements[ielem][1] - 1;
                ind3 = mesh.NeighbourElements[ielem][2] - 1;
                vector<double> u = {state_vectors[ielem][0], state_vectors[ind1][0],
                                    state_vectors[ind2][0], state_vectors[ind3][0]};
                double u1 = u[1], u2 = u[2], u3 = u[3];
                double u_m = u[0];
                double x1, y1, x2, y2, x3, y3, xm, ym, Dx, Dy;
                xm = mesh.Centroid[ielem][0];   ym = mesh.Centroid[ielem][1];
                x1 = mesh.Centroid[ind1][0];    y1 = mesh.Centroid[ind1][1];
                x2 = mesh.Centroid[ind2][0];    y2 = mesh.Centroid[ind2][1];
                x3 = mesh.Centroid[ind3][0];    y3 = mesh.Centroid[ind3][1];
                Dx = gradu[ielem][0][0]; Dy = gradu[ielem][0][1];
                A[0][0] = (x1 - xm) * Dx;   A[0][1] = (y1 - ym) * Dy;
                A[1][0] = (x2 - xm) * Dx;   A[1][1] = (y2 - ym) * Dy;
                A[2][0] = (x3 - xm) * Dx;   A[2][1] = (y3 - ym) * Dy;
                A[3][0] = -(x1 - xm) * Dx;  A[3][1] = -(y1 - ym) * Dy;
                A[4][0] = -(x2 - xm) * Dx;  A[4][1] = -(y2 - ym) * Dy;
                A[5][0] = -(x3 - xm) * Dx;  A[5][1] = -(y3 - ym) * Dy;
                A[6][0] = 1.0;              A[6][1] = 0.0;
                A[7][0] = 0.0;              A[7][1] = 1.0;
                A[8][0] = -1.0;             A[8][1] = 0.0;
                A[9][0] = 0.0;              A[9][1] = -1.0;

                c[0] = abs(Dx); c[1] = abs(Dy);

                b[0] = utils::min(u_m, u1) - u_m; b[1] = utils::min(u_m, u2) - u_m; b[2] = utils::min(u_m, u3) - u_m;
                b[3] = u_m - utils::max(u_m, u1); b[4] = u_m - utils::max(u_m, u2); b[5] = u_m - utils::max(u_m, u3);

                b[6] = 0.0; b[7] = 0.0; b[8] = -1.0; b[9] = -1.0;

                vector<double> phi(2, 0.0);
                vector<int> W0 = {6, 7};
                vector<double> x0 = {0.0, 0.0};
                phi = LP_Optimizer(c, A, b, W0, x0);
                limiter[ielem][0] = phi[0];
                limiter[ielem][1] = phi[1];
            } else
            {
                // Do not apply limiter on boundary elements
                limiter[ielem][0] = 1.0; limiter[ielem][1] = 1.0;
            }
        }
    }
    if (strcasecmp(param.name_limiter.c_str(), "BJ") == 0)
    {
        // Barth and Jespersen Limiter
        // Calculate limiter based on density
        vector<double> umax(mesh.E.size());
        vector<double> umin(mesh.E.size());
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            umax[ielem] = state_vectors[ielem][0];
            umin[ielem] = state_vectors[ielem][0];
        }
        // loop through interior elements;
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            if (mesh.NeighbourElements[ielem][0] > 0 && mesh.NeighbourElements[ielem][1] > 0 && mesh.NeighbourElements[ielem][2] > 0)
            {
                // If this is an interior element
                // Construct the matrices for LP Optimization
                int ind1, ind2, ind3;
                ind1 = mesh.NeighbourElements[ielem][0] - 1;
                ind2 = mesh.NeighbourElements[ielem][1] - 1;
                ind3 = mesh.NeighbourElements[ielem][2] - 1;
                vector<double> u = {state_vectors[ielem][0], state_vectors[ind1][0],
                                    state_vectors[ind2][0], state_vectors[ind3][0]};
                auto result  = std::minmax_element(u.begin(), u.end());
                double u_min = *result.first;
                double u_max = *result.second;
                umax[ielem] = u_max;
                umin[ielem] = u_min;
            }
        }
        int index_edge_nodes[2];
        double u_left[2], u_right[2];
        vector<vector<double> > phi_temp(mesh.E.size(), vector<double> (2, 1.0));
        for (int iedge = 0; iedge < mesh.I2E.size(); iedge++)
        {
            int index_left = mesh.I2E[iedge][0] - 1;
            int index_right = mesh.I2E[iedge][2] - 1;
            /* Get the local index of the edge points*/
            utils::mod3(mesh.I2E[iedge][1] - 1, index_edge_nodes);
            int index_pA = mesh.E[index_left][index_edge_nodes[0]] - 1;
            int index_pB = mesh.E[index_left][index_edge_nodes[1]] - 1;
            // vertex coordinates
            double coord_pA[2], coord_pB[2];
            coord_pA[0] = mesh.V[index_pA][0]; coord_pA[1] = mesh.V[index_pA][1];
            coord_pB[0] = mesh.V[index_pB][0]; coord_pB[1] = mesh.V[index_pB][1];
            // get left state
            u_left[0] = state_vectors[index_left][0] +
                                (gradu[index_left][0][0] * (coord_pA[0] - mesh.Centroid[index_left][0]) +
                                gradu[index_left][0][1] * (coord_pA[1] - mesh.Centroid[index_left][1]));
            u_left[1] = state_vectors[index_left][0] +
                                (gradu[index_left][0][0] * (coord_pB[0] - mesh.Centroid[index_left][0]) +
                                gradu[index_left][0][1] * (coord_pB[1] - mesh.Centroid[index_left][1]));
            // get right state
            u_right[0] = state_vectors[index_right][0] +
                                (gradu[index_right][0][0] * (coord_pA[0] - mesh.Centroid[index_right][0]) +
                                gradu[index_right][0][1] * (coord_pA[1] - mesh.Centroid[index_right][1]));
            u_right[1] = state_vectors[index_right][0] +
                                (gradu[index_right][0][0] * (coord_pB[0] - mesh.Centroid[index_right][0]) +
                                gradu[index_right][0][1] * (coord_pB[1] - mesh.Centroid[index_right][1]));

            for (int ivertex = 0; ivertex < 2; ++ivertex)
            {
                /* Left element */
                if (u_left[ivertex] - state_vectors[index_left][0] > 1e-11)
                {
                    phi_temp[index_left][ivertex] = min(1.0, (umax[index_left] - state_vectors[index_left][0]) / (u_left[ivertex] - state_vectors[index_left][0]));
                } else if (u_left[ivertex] - state_vectors[index_left][0] < -1e-11)
                {
                    phi_temp[index_left][ivertex] = min(1.0, (umin[index_left] - state_vectors[index_left][0]) / (u_left[ivertex] - state_vectors[index_left][0]));
                }

                /* Right element */
                if (u_right[ivertex] - state_vectors[index_right][0] > 1e-11)
                {
                    phi_temp[index_right][ivertex] = min(1.0, (umax[index_right] - state_vectors[index_right][0]) / (u_right[ivertex] - state_vectors[index_right][0]));
                } else if (u_right[ivertex] - state_vectors[index_right][0] < -1e-11)
                {
                    phi_temp[index_right][ivertex] = min(1.0, (umin[index_right] - state_vectors[index_right][0]) / (u_right[ivertex] - state_vectors[index_right][0]));
                }

            }
        }
        for (int ielem = 0; ielem < mesh.E.size(); ielem++)
        {
            limiter[ielem][0] = *min_element(phi_temp[ielem].begin(), phi_temp[ielem].end());
            limiter[ielem][1] = limiter[ielem][0];
        }
    }
    return limiter;
}

vector<double> LP_Optimizer(vector<double> c, vector<vector<double> > A, vector<double> b, vector<int> W0, vector<double> x0)
{
    vector<double> phi = x0; // the optimization result returned
    int MAXITER = 10; // The maximum iteration steps of the linear programming
    // A has the structure that the 1st, 2nd and 3rd rows are linear independent
    // 4th, 5th and 6th rows are negative of 1st, 2nd and 3rd rows
    // According to the paper, (0, 0) is a good starting point， which is in the initialize of phi
    /* Allocate the working set */
    vector<int> Wk = W0; // The initial working set: 6 and 7
    vector<vector<double> > Ak(2, vector<double> (2, 0.0));
    vector<vector<double> > inv_Ak_trans(2, vector<double> (2, 0.0));
    vector<vector<double> > inv_Ak(2, vector<double> (2, 0.0));
    vector<double> bk(2, 0.0);
    // The temp vectors in optimization process
    vector<double> lambda_k(2, 0.0);
    vector<double> p_k(2, 0.0);
    vector<int> D_k;
    int ind_q = 0;
    // Assembling the working matrices
    for (int i = 0; i < 2; i++)
    {
        Ak[i][0] = A[Wk[i]][0]; Ak[i][1] = A[Wk[i]][1];
        bk[i] = b[Wk[i]];
    }
    for (int niter = 0; niter < MAXITER; niter++)
    {
        D_k.clear();
        // Solve lagrange multipliers lambda_k
        inv_Ak_trans = utils::Invert22Matrix(utils::Trans22Matrix(Ak));
        lambda_k[0] = inv_Ak_trans[0][0] * c[0] + inv_Ak_trans[0][1] * c[1];
        lambda_k[1] = inv_Ak_trans[1][0] * c[0] + inv_Ak_trans[1][1] * c[1];
        if ((lambda_k[0] > 0 && lambda_k[1] > 0) || ((abs(lambda_k[0]) < 1e-6) && (abs(lambda_k[1]) < 1e-6)))
            break; // this time phi is optimal
        if (lambda_k[0] < lambda_k[1]) {
            ind_q = 0;
        }else {
            ind_q = 1;
        }
        // Calculate the descent direction p_k
        inv_Ak = utils::Invert22Matrix(Ak);
        p_k[0] = inv_Ak[0][0] * (1 - ind_q) + inv_Ak[0][1] * ind_q;
        p_k[1] = inv_Ak[1][0] * (1 - ind_q) + inv_Ak[1][1] * ind_q;
        // Find the set of decreasing constraints along p_k
        for (int i = 0; i < A.size(); i++)
        {
            double ai_pk;
            ai_pk = A[i][0] * p_k[0] + A[i][1] * p_k[1];
            if (ai_pk < -1e-11)
            {
                D_k.push_back(i);
            }
        }
        // Calculate the max step length gamma
        vector<double> gamma(D_k.size(), 0.0);
        double alpha_k = std::numeric_limits<double>::infinity();
        for (int i = 0; i < gamma.size(); i++)
        {
            if (D_k[i] != Wk[0] && D_k[i] != Wk[1]) {
                double ai_pk, ai_phik;
                ai_pk = A[D_k[i]][0] * p_k[0] + A[D_k[i]][1] * p_k[1];
                ai_phik = A[D_k[i]][0] * phi[0] + A[D_k[i]][1] * phi[1];
                gamma[i] = (ai_phik - b[D_k[i]]) / (-ai_pk);
                if (gamma[i] < alpha_k) {
                    alpha_k = gamma[i];
                }
            }
        }
        // Update the working set
        phi[0] = phi[0] + alpha_k * p_k[0]; phi[1] = phi[1] + alpha_k * p_k[1];
        vector<int> S_k;
        for (int i = 0; i < gamma.size(); i++)
        {
            if (abs(gamma[i] - alpha_k) < 1e-5)
            {
                S_k.push_back(D_k[i]);
            }
        }
        // choose the first element of S_k
        int t = S_k[S_k.size() - 1];
        Wk[ind_q] = t;
        // Assembling the working matrices
        for (int i = 0; i < 2; i++)
        {
            Ak[i][0] = A[Wk[i]][0]; Ak[i][1] = A[Wk[i]][1];
            bk[i] = b[Wk[i]];
        }
    }
    return phi;
}
