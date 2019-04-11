#include "../include/CalcNumericalFlux.h"

using namespace std;

vector<double> projection_2d(vector<vector<double> > F, vector<double> n){

    vector<double> F_hat(4, 0.0);

    F_hat[0] = F[0][0] * n[0] + F[0][1] * n[1];
    F_hat[1] = F[1][0] * n[0] + F[1][1] * n[1];
    F_hat[2] = F[2][0] * n[0] + F[2][1] * n[1];
    F_hat[3] = F[3][0] * n[0] + F[3][1] * n[1];

    return F_hat;

}

vector<vector<double> > flux_function_2d(vector<double> u, double gamma){

    double p, H;
    vector<vector<double> > F(4, vector<double> (2, 0.0));

    p = (gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
    H = u[3] / u[0] + p / u[0];

    F[0][0] = u[1]; F[1][0] = u[1] * u[1] / u[0] + p;
    F[2][0] = u[1] * u[2] / u[0]; F[3][0] = u[1] * H;

    F[0][1] = u[2]; F[1][1] = u[2] * u[1] / u[0];
    F[2][1] = u[2] * u[2] / u[0] + p; F[3][1] = u[2] * H;

    return F;

}

vector<double> CalcNumericalFlux(vector<double> uL, vector<double> uR, vector<double> n, char const *type_flux, double gamma, double& mws){

    if (strcasecmp(type_flux, "roe") == 0) {

        vector<vector<double> > FL(4, vector<double>(2, 0.0));
        vector<vector<double> > FR(4, vector<double>(2, 0.0));
        vector<double> FL_hat(4, 0.0), FR_hat(4, 0.0), F_hat(4, 0.0); // projected flux function on normal vector

        double v_roe_vec[2];
        double u_roe;
        double H_roe, c_roe, q_roe;
        double pL, pR, HL, HR, rhoL, rhoR;
        double lambda_1, lambda_2, lambda_3;
        double s1, s2, G1, G2, C1, C2;

        double drho, drhoE, drhov_vec[2];

        //Primitive variables
        rhoL = uL[0];
        rhoR = uR[0];
        pL = (gamma - 1) * (uL[3] - 0.5 * (uL[1] * uL[1] + uL[2] * uL[2]) / uL[0]);
        HL = uL[3] / uL[0] + pL / uL[0];
        pR = (gamma - 1) * (uR[3] - 0.5 * (uR[1] * uR[1] + uR[2] * uR[2]) / uR[0]);

        HR = uR[3] / uR[0] + pR / uR[0];
        // calculate v_roe
        v_roe_vec[0] = (sqrt(rhoL) * (uL[1] / uL[0]) + sqrt(rhoR) * (uR[1] / uR[0])) / (sqrt(rhoL) + sqrt(rhoR));
        v_roe_vec[1] = (sqrt(rhoL) * (uL[2] / uL[0]) + sqrt(rhoR) * (uR[2] / uR[0])) / (sqrt(rhoL) + sqrt(rhoR));
        u_roe = v_roe_vec[0] * n[0] + v_roe_vec[1] * n[1];
        // calculate H_roe
        H_roe = (sqrt(rhoL) * HL + sqrt(rhoR) * HR) / (sqrt(rhoL) + sqrt(rhoR));

        // velocity magnitude
        q_roe = sqrt(v_roe_vec[0] * v_roe_vec[0] + v_roe_vec[1] * v_roe_vec[1]);
        // sound speed
        c_roe = sqrt((gamma - 1.0) * (H_roe - 0.5 * (v_roe_vec[0] * v_roe_vec[0] + v_roe_vec[1] * v_roe_vec[1])));
        // eigenvalues
        double lambda[3];
        lambda[0] = u_roe + c_roe;
        lambda[1] = u_roe - c_roe;
        lambda[2] = u_roe;

        // entropy fix for lambda_1, lambda_2 and lambda_3
        double epsilon = 0.1 * c_roe;
        int i;
        for (i = 0; i < 3; i++) {
            if (lambda[i] < epsilon && lambda[i] > -epsilon) {
                lambda[i] = 0.5 * (epsilon + lambda[i] * lambda[i] / epsilon);
            }
        }
        lambda_1 = lambda[0];
        lambda_2 = lambda[1];
        lambda_3 = lambda[2];

        // calculate the delta quantities for roe flux
        drho = uR[0] - uL[0];
        drhoE = uR[3] - uL[3];
        drhov_vec[0] = uR[1] - uL[1];
        drhov_vec[1] = uR[2] - uL[2];
        // temporary vars for calculating the Roe Flux
        s1 = 0.5 * (fabs(lambda_1) + fabs(lambda_2));
        s2 = 0.5 * (fabs(lambda_1) - fabs(lambda_2));
        G1 = (gamma - 1.0) *
             (q_roe * q_roe * drho / 2 + drhoE - (v_roe_vec[0] * drhov_vec[0] + v_roe_vec[1] * drhov_vec[1]));
        G2 = -1.0 * u_roe * drho + (drhov_vec[0] * n[0] + drhov_vec[1] * n[1]);
        C1 = G1 * (s1 - fabs(lambda_3)) / pow(c_roe, 2) + G2 * s2 / c_roe;
        C2 = G1 * s2 / c_roe + (s1 - fabs(lambda_3)) * G2;

        // Calculate the Roe flux
        FL = flux_function_2d(uL, gamma);
        FR = flux_function_2d(uR, gamma);

        FL_hat = projection_2d(FL, n);
        FR_hat = projection_2d(FR, n);


        F_hat[0] = 0.5 * (FL_hat[0] + FR_hat[0]) - 0.5 * (fabs(lambda_3) * drho + C1);
        F_hat[1] =
                0.5 * (FL_hat[1] + FR_hat[1]) - 0.5 * (fabs(lambda_3) * drhov_vec[0] + C1 * v_roe_vec[0] + C2 * n[0]);
        F_hat[2] =
                0.5 * (FL_hat[2] + FR_hat[2]) - 0.5 * (fabs(lambda_3) * drhov_vec[1] + C1 * v_roe_vec[1] + C2 * n[1]);
        F_hat[3] = 0.5 * (FL_hat[3] + FR_hat[3]) - 0.5 * (fabs(lambda_3) * drhoE + C1 * H_roe + C2 * u_roe);

        mws = fabs(u_roe) + c_roe;

        return F_hat;
    }
    else{
        std::cout << "Unsupport flux name: " << type_flux << " Aborting" << std::endl;
        abort();
    }
}
