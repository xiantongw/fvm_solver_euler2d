#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h> 

void boudary_condition(double u[4], double norm[2], char* boundary_type, CfdParam* cparam, double num_flux[4], double *mws){ 

    if (strcasecmp(boundary_type, "Inflow") == 0){

        double p = (cparam->gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
        double un = (u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1];
        double c = sqrt(cparam->gamma * p / u[0]);
        double J = un + 2.0 * c / (cparam->gamma - 1); // Riemann Invariant
        double dn = cos(cparam->attack_angle)*norm[0] + sin(cparam->attack_angle)*norm[1];
        double R = 1.0;

        double Tt = 1.0 + 0.5 * (cparam->gamma - 1) * cparam->mach_inf * cparam->mach_inf;
        double pt = pow(Tt, cparam->gamma / (cparam->gamma - 1.0));

        // Solve for Mb
        // ca, cb, cc are coefficients for quadratic equation
        double tmpa = cparam->gamma * R * Tt * dn * dn - 0.5 * (cparam->gamma - 1.0) * J * J;
        double tmpb = 4.0 * cparam->gamma * R * Tt * dn / (cparam->gamma - 1.0);
        double tmpc = 4.0 * cparam->gamma * R * Tt  / pow((cparam->gamma - 1.0), 2) - J * J;
        double Mb1 = (-1.0 * tmpb - sqrt(tmpb * tmpb - 4.0 * tmpa * tmpc)) / (2 * tmpa);
        double Mb2 = (-1.0 * tmpb + sqrt(tmpb * tmpb - 4.0 * tmpa * tmpc)) / (2 * tmpa);
        double Mb;
        if (Mb1 < 0)
            Mb = Mb2;
        else
            Mb = Mb1;
        // Calculate the exterior states
        double Tb = Tt / (1.0 + 0.5*(cparam->gamma - 1.0) * Mb * Mb);
        double pb = pt * pow(Tb / Tt, cparam->gamma / (cparam->gamma - 1.0));
        double rhob = pb / (R * Tb);
        double cb = sqrt(cparam->gamma * pb / rhob);
        double vb[2] = {Mb * cb * cos(cparam->attack_angle), Mb * cb * sin(cparam->attack_angle)};
        double rhoEb = pb / (cparam->gamma - 1.0) + 0.5 * rhob * (vb[0] * vb[0] + vb[1] * vb[1]);

        double ub[4] = {rhob, rhob * vb[0], rhob * vb[1], rhoEb};
        double F_b[4][2];
        flux_function_2d(ub, cparam->gamma, F_b);
        projection_2d(F_b, norm, num_flux); // num_flux is passed out

        // Max Wave Speed
        *mws = fabs(vb[0] * norm[0] + vb[1] * norm[1]) + cb;

    } else if (strcasecmp(boundary_type, "Inviscid_Wall") == 0){

        double vb[2];
        
        vb[0] = u[1] / u[0] - ((u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1]) * norm[0];
        vb[1] = u[2] / u[0] - ((u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1]) * norm[1];

        double pb = (cparam->gamma - 1.0) * (u[3] - 0.5 * u[0] * (vb[0] * vb[0] + vb[1] * vb[1]));

        num_flux[0] = 0.0;
        num_flux[1] = norm[0] * pb;
        num_flux[2] = norm[1] * pb;
        num_flux[3] = 0.0;

        double p = (cparam->gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
        *mws = sqrt(pow(u[1] / u[0], 2) + pow(u[2] / u[0], 2)) + sqrt(cparam->gamma * p / u[0]);

    } else if (strcasecmp(boundary_type, "Subsonic_Outflow") == 0){

        /* Interior entropy*/
        double p = (cparam->gamma - 1.0) * (u[3] - 0.5 * (u[1]*u[1] + u[2]*u[2]) / u[0]);
        double S = p / pow(u[0], cparam->gamma);
        double pb = cparam->p_inf;
        double rhob = pow(pb / S,  1.0 / cparam->gamma);
        double cb = sqrt(cparam->gamma * pb / rhob);
        double un = (u[1] / u[0]) * norm[0] + (u[2] / u[0]) * norm[1];
        
        double c = sqrt(cparam->gamma * p / u[0]);
        double J = un + 2.0 * c / (cparam->gamma - 1.0); // Riemann Invariant

        double ub_n = J - 2.0 * cb / (cparam->gamma - 1.0);

        /* Solve for vb*/
        double vb[2] = {0.0, 0.0};
        vb[0] = u[1] / u[0] - norm[0] * (u[1]*norm[0] + u[2]*norm[1]) / u[0] + ub_n * norm[0];
        vb[1] = u[2] / u[0] - norm[1] * (u[1]*norm[0] + u[2]*norm[1]) / u[0] + ub_n * norm[1];

        double rhoEb = pb / (cparam->gamma - 1.0) + 0.5 * rhob * (vb[0] * vb[0] + vb[1] * vb[1]);

        double ub[4] = {rhob, rhob * vb[0], rhob * vb[1], rhoEb};
        
        double F_b[4][2];
        flux_function_2d(ub, cparam->gamma, F_b);
        projection_2d(F_b, norm, num_flux); // num_flux is passed out
        
        *mws = sqrt(pow(ub[1] / ub[0], 2) + pow(ub[2] / ub[0], 2)) + sqrt(cparam->gamma * pb / ub[0]);

    } else if (strcasecmp(boundary_type, "Free_Stream") == 0){
        
        double mws_temp;
        double u_free[4];

        u_free[0] = 1.0;
        u_free[1] = cparam->mach_inf * cos(cparam->attack_angle);
        u_free[2] = cparam->mach_inf * sin(cparam->attack_angle);
        u_free[3] = 1 / (cparam->gamma * (cparam->gamma - 1)) + 
                        0.5 * cparam->mach_inf * cparam->mach_inf;

        roe_2d(u, u_free, norm, cparam->gamma, num_flux, &mws_temp);
        *mws = mws_temp;
    
    } else{
        printf("ERROR: Unknown Boundary Condition: %s \n", boundary_type);
        abort();
    }
}
