//# include "FVM_header.h"
# include "Limiter_header.h"


vector<double> perimeter(const vector<double> Ilength, const vector<double> Blength, const vector<vector<int>> I2E, const vector<vector<int>> B2E, int nelem) {
	vector<double> Elem_perimeter(nelem);
	for (int ind = 0; ind < nelem; ind++) {
		Elem_perimeter[ind] = 0;
	}

	for (int ind = 0; ind < Ilength.size(); ind++) {
		Elem_perimeter[I2E[ind][0] - 1] = Elem_perimeter[I2E[ind][0] - 1] + Ilength[ind];
		Elem_perimeter[I2E[ind][2] - 1] = Elem_perimeter[I2E[ind][2] - 1] + Ilength[ind];
	}

	for (int ind = 0; ind < Blength.size(); ind++) {
		Elem_perimeter[B2E[ind][0] - 1] = Elem_perimeter[B2E[ind][0] - 1] + Blength[ind];
	}


	return Elem_perimeter;
}


vector<double> Roe_Flux(vector<double> u_L, vector<double> u_R, vector<double> n_elem, double Gamma, double gas_const_R) {
/*	
	inputs: the state vectors of left and right elements; 
	n_elem: the normal vector of the face from left to right;
	Gamma, gas_const_R: constants
*/
	vector<vector<double>> F_L(4,vector<double>(2));
	vector<vector<double>> F_R(4, vector<double>(2));
	
	vector <double> v_L = {u_L[1]/u_L[0], u_L[2]/u_L[0]};
	vector <double> v_R = {u_R[1] / u_R[0], u_R[2] / u_R[0] };
	
	double p_L = (Gamma - 1)*(u_L[3]-0.5*u_L[0]*(pow(v_L[0],2.0)+pow(v_L[1],2.0)));
	double p_R = (Gamma - 1)*(u_R[3] - 0.5*u_R[0] * (pow(v_R[0], 2.0) + pow(v_R[1], 2.0)));

	double H_L = (u_L[3] + p_L) / u_L[0];
	double H_R = (u_R[3] + p_R) / u_R[0];
	//mexPrintf("test: \npL, HL: %lf %lf\n:pR, HR %lf %lf\n", p_L, H_L, p_R, H_R);
	F_L = state2flux(u_L, Gamma);
	F_R = state2flux(u_R, Gamma);


	vector<double> vavg_L =  scal_multip(v_L,sqrt(u_L[0]));
	vector<double> vavg_R = scal_multip(v_R, sqrt(u_R[0]));
	vector<double> v_nom = vec_add(vavg_L, vavg_R);
	vector<double> v_avg = scal_multip( v_nom ,1/(sqrt(u_L[0]) + sqrt(u_R[0])));
	double H_avg = (sqrt(u_L[0])*H_L + sqrt(u_R[0])*H_R)/(sqrt(u_L[0])+ sqrt(u_R[0]));
	
	//mexPrintf("test: \nv_avg: %lf %lf\nH_avg: %lf\n", v_avg[0], v_avg[1], H_avg);


	double norm_vavg_square = pow(v_avg[0], 2.0) + pow(v_avg[1], 2.0);
	double c = (Gamma - 1)*(H_avg - 1.0 / 2.0 * (norm_vavg_square));
	
	if (c<=0) {
		throw std::invalid_argument("Info input leads to imaginary or 0 speed of sound");
	}

	c = abs(sqrt(c));
	double u_speed = dot_prod(v_avg, n_elem);

	//mexPrintf("normal speed: %lf,  speed of sound: %lf,\n", u_speed, c);
	double entropy_fix = 0.1*c;

	vector<double> eig_values = { u_speed + c, u_speed - c, u_speed, u_speed };
	
	// entropy fixing
	for (int ind = 0; ind < 4; ind++) {
		if (eig_values[ind]<entropy_fix && eig_values[ind]>(-1.0)*entropy_fix) {
			eig_values[ind] = (pow(entropy_fix,2.0)+pow(eig_values[ind],2)) / (2.0*entropy_fix);
		}
	}

	double s1 = 0.5*(abs(eig_values[0]) + abs(eig_values[1]));
	double s2 = 0.5*(abs(eig_values[0]) - abs(eig_values[1]));
	
	vector<double> rhou_diff = {u_R[1]-u_L[1], u_R[2] - u_L[2]};
	double rhoE_diff = u_R[3] - u_L[3];
	double rho_diff = u_R[0] - u_L[0];

	double G1 = (Gamma - 1.0) * (norm_vavg_square/2.0*rho_diff - dot_prod(v_avg,rhou_diff) + rhoE_diff);
	double G2 = (-1)*u_speed*rho_diff + dot_prod(rhou_diff, n_elem);
	
	double C1 = G1 / pow(c, 2.0)*(s1 - abs(eig_values[2])) + G2 / c*s2;
	double C2 = G1 / c*s2 + G2*(s1 - abs(eig_values[2]));

	vector<double> F_L_n = { dot_prod(F_L[0],n_elem), dot_prod(F_L[1],n_elem),dot_prod(F_L[2],n_elem),dot_prod(F_L[3],n_elem) };
	vector<double> F_R_n = { dot_prod(F_R[0],n_elem), dot_prod(F_R[1],n_elem),dot_prod(F_R[2],n_elem),dot_prod(F_R[3],n_elem) };

	vector<double> F_correct = {
		abs(eig_values[2])*rho_diff + C1,
		abs(eig_values[2])*rhou_diff[0] + C1*v_avg[0] + C2*n_elem[0],
		abs(eig_values[2])*rhou_diff[1] + C1*v_avg[1] + C2*n_elem[1],
		abs(eig_values[2])*rhoE_diff + C1*H_avg + C2*u_speed
	};
	//mexPrintf("Upwind State: %lf,%lf,%lf,%lf\n", F_L_n[0], F_L_n[1], F_L_n[2], F_L_n[3]);
	vector<double> temp1 = scal_multip(vec_add(F_L_n, F_R_n), 0.5);
	vector<double> temp2 = scal_multip(F_correct, -0.5);
	vector<double> Fhat = vec_add(temp1, temp2);
	
	vector<double> max_eigen = {abs(u_speed)+c, abs(u_speed)-c ,abs(eig_values[2]) ,abs(eig_values[3]) };
	double max_eig = find_max(max_eigen,false);
	Fhat.push_back(max_eig); // add wave speed

	return Fhat;

}


vector<double> scal_multip(vector<double> vec, double c) {
	vector<double> res_vec(vec.size());
	for (int ind = 0; ind < vec.size(); ind++) {
		res_vec[ind] = vec[ind] * c;
	}
	return res_vec;
}


vector<double> vec_add(vector<double> a, vector<double> b) {
	// check size:
	vector<double> res_vec(a.size());
	if (a.size() != b.size()) {
		throw std::invalid_argument("Vector of different size. Cannot add.");
	}
	else {
		for (int ind = 0; ind < a.size(); ind++) {
			res_vec[ind] = a[ind] + b[ind];
		}
	}

	return res_vec;
}

double dot_prod(vector<double> a, vector<double> b) {
	double res = 0;
	// check size
	if (a.size() != b.size()) {
		throw std::invalid_argument("Vector of different size. Cannot add.");
	}
	else {
		for (int ind = 0; ind < a.size(); ind++) {
			res += a[ind] * b[ind];
		}
	}

	return res;
}


vector<vector<double>> state2flux(vector<double> u, double Gamma) {
	vector<vector<double>> F(4, vector<double>(2));

	vector <double> v = { u[1] / u[0], u[2] / u[0] };

	double p = (Gamma - 1)*(u[3] - 0.5*u[0] * (pow(v[0], 2.0) + pow(v[1], 2.0)));


	double H = (u[3] + p) / u[0];

	F[0][0] = u[1]; // rho*u
	F[0][1] = u[2]; // rho*v
	F[1][0] = pow(u[1], 2.0) / u[0] + p; //rho*u**2 + p
	F[1][1] = u[1] * u[2] / u[0]; // rho*u*v
	F[2][0] = F[1][1];
	F[2][1] = pow(u[2], 2.0) / u[0] + p; //rho*v**2 + p
	F[3][0] = u[1] * H; //rho*u*H
	F[3][1] = u[2] * H; //rho*v*H
	return F;
}


void Roe_flux_tester(double Gamma) {
	mexPrintf("\nRoe Flux Tester\n\n");
	//vector<double> In = { 1 / sqrt(2.0), 1 / sqrt(2.0) };
	vector<double> In = {0,1};
	
	//vector<double> In = { 0, -1};
	double Gas_Const_R = 1;
	/* Supersonic case*/
	//vector<double> u_L = { 1,2,3,15 };
	//vector<double> u_R = { 1, 3, 2, 14 };

	//vector<double> u_L = {0.99274128, 0.60694781, -0.01611135, 2.66032805};
	//vector<double> u_R = {0.99768788, 0.59653994, -0.02375726, 2.67036239};
	
	//vector<double> u_L = { 1.000000, 0.591608, 0.000000, 2.675000 };
	//vector<double> u_R = { 1.000000 ,0.591608 ,0.000000, 2.675000 };

	vector<double> u_L = { 0.7, 0.5, 0.6, 2.0 };
	vector<double> u_R = scal_multip(u_L, 0.9);

	/*
	vector<double> test1 = Roe_Flux(u_L, u_L, In, Gamma, Gas_Const_R);
	vector<vector<double>> flux_test = state2flux(u_L, Gamma);
	vector<double> test1_rhs = { dot_prod(flux_test[0],In),dot_prod(flux_test[1],In) ,dot_prod(flux_test[2],In) ,dot_prod(flux_test[3],In) };
	mexPrintf("Test 1: Func: %lf %lf %lf %lf\n", test1[0], test1[1], test1[2], test1[3]);
	mexPrintf("Test 1: dotprod: %lf %lf %lf %lf\n", test1_rhs[0], test1_rhs[1], test1_rhs[2], test1_rhs[3]);
	*/
	//test 2:
	vector<double> test2 = Roe_Flux(u_L, u_R, In, Gamma, Gas_Const_R);
	vector<double> test2_rhs = Roe_Flux(u_R, u_L, scal_multip(In, -1), Gamma, Gas_Const_R);
	mexPrintf("Test 2: regular: %lf %lf %lf %lf\n", test2[0], test2[1], test2[2], test2[3]);
	mexPrintf("Test 2: change dir: %lf %lf %lf %lf\n", test2_rhs[0], test2_rhs[1], test2_rhs[2], test2_rhs[3]);
	mexPrintf("\n\n\n");
}


vector<double> Boundary_Flux(vector<double> u_in, vector<double> b_n, int b_type, double Gamma, double Gas_Const_R) {
	/*Note: boundary states defined here! Change the test state here*/
	// u_in: the state inside the computational domain;
	// b_n: normal vector pointing outward of the boundary face
	// byype: the type of the boundary condition. 
	double p_inf = 1.0; 
	double M_inf = 0.8; // Only used for free_stream tests, i.e., full state boundary
	double T_t = 1.0 + (Gamma - 1.0) / 2.0 * pow(M_inf,2.0);
	double p_t = pow(T_t, Gamma / (Gamma - 1));
	double alpha = 0;


	vector<double> flux_out(5);//change the memory definition

	
	if (b_type == 1) { // Full State Computation. Here only for the free_stream test
		double T_inf, rho_inf, c_inf, u_speed, rhoE;
		vector<double> u_out(4);

		T_inf = T_t / (1.0 + (Gamma - 1.0) / 2.0 * pow(M_inf, 2.0));
		rho_inf = p_inf / (Gas_Const_R*T_inf);
		c_inf = sqrt(Gamma *Gas_Const_R *T_inf);
		u_speed = c_inf * M_inf;
		rhoE = p_inf / (Gamma - 1.0) + 0.5 * rho_inf * pow(u_speed, 2.0);
		//system("PAUSE");
		u_out[0] = rho_inf; u_out[1] = rho_inf*u_speed; u_out[2] = 0; u_out[3] = rhoE; 
		//mexPrintf("Uout: %lf %lf %lf %lf\n",u_out[0], u_out[1], u_out[2], u_out[3]);
		//mexPrintf("Uin: %lf %lf %lf %lf\n", u_in[0], u_in[1], u_in[2], u_in[3]);
		//mexPrintf("bn:%lf %lf\n", b_n[0], b_n[1]);
		//system("PAUSE");

		flux_out = Roe_Flux(u_in, u_out, b_n, Gamma, Gas_Const_R);
		
	}

	else if (b_type == 2) { //invisid surface
		vector<double> v_b;
		double p_b;
		vector<double> v_in = { u_in[1] / u_in[0], u_in[2] / u_in[0] };
		v_b = vec_add(v_in, scal_multip(b_n, (-1)*dot_prod(v_in, b_n)));
		p_b = (Gamma - 1) * (u_in[3] - 0.5 * u_in[0] * (pow(v_b[0], 2.0) + pow(v_b[1], 2.0)));
		flux_out[0] = 0; flux_out[1] = p_b * b_n[0]; flux_out[2] = p_b*b_n[1]; flux_out[3] = 0; 
		flux_out[4] = 0; // no wavespeed across boundary
	}

	else if (b_type == 3) { //subsonic inflow
		vector<double> v_in = { u_in[1] / u_in[0], u_in[2] / u_in[0] };
		vector<double> n_in = { cos(alpha), sin(alpha) };
		
		double p_in = (Gamma - 1)*(u_in[3] - 0.5*u_in[0]*(pow(v_in[0],2.0) + pow(v_in[1],2.0)));
		double c_in = sqrt(Gamma*p_in / u_in[0]);

		double J_plus = dot_prod(v_in, b_n) + 2 * c_in/(Gamma - 1);

		double dn = dot_prod(n_in, b_n);

		double temp1, temp2, temp3;
		temp1 = Gamma*Gas_Const_R*T_t*pow(dn, 2.0) - (Gamma - 1) / 2.0*pow(J_plus, 2.0);
		temp2 = 4 * Gamma*Gas_Const_R*T_t*dn / (Gamma - 1);
		temp3 = 4 * Gamma*Gas_Const_R*T_t / (pow(Gamma - 1, 2.0)) - pow(J_plus,2);
		
		vector<double> M_mtx = solve_quadratic(temp1, temp2, temp3);

		double M_b;
		if (M_mtx[0] > 0 && M_mtx[0] < 1) {
			M_b = M_mtx[0];
		}
		else {
			M_b = M_mtx[1];
		}
		
		double T_b = T_t / (1.0 + (Gamma - 1.0) / 2.0 * pow(M_b, 2.0));
		double p_b = p_t*pow(T_b / T_t, Gamma / (Gamma - 1));
		double rho_b = p_b / (Gas_Const_R*T_b);

		double c_b = sqrt(Gamma*p_b / rho_b);

		//double c_b = sqrt(Gamma*Gas_Const_R*T_t/(1.0 + (Gamma - 1.0) / 2.0 * pow(M_b, 2.0)));
		
		vector<double> v_b = scal_multip(n_in, M_b*c_b);

		double rhoE = p_b / (Gamma - 1) + 0.5 * rho_b*(pow(v_b[0],2.0) + pow(v_b[1],2.0));

		vector<double> b_state = { rho_b, rho_b*v_b[0], rho_b * v_b[1], rhoE };
		vector<vector<double>> b_flux = state2flux(b_state, Gamma);
		
		flux_out[0] = dot_prod(b_flux[0], b_n);  flux_out[1] = dot_prod(b_flux[1], b_n);  flux_out[2] = dot_prod(b_flux[2], b_n);  flux_out[3] = dot_prod(b_flux[3], b_n);
		vector<double> wavespeed = {abs(dot_prod(v_in,b_n))+c_in, abs(dot_prod(v_b,b_n)) + c_b };
		flux_out[4] = find_max(wavespeed, false); // wave speed

	}
	else if (b_type == 4) { //subsonic outflow
		double p_b = p_inf;
		vector<double> v_in = { u_in[1] / u_in[0], u_in[2] / u_in[0] };
		double p_in = (Gamma - 1)*(u_in[3] - 0.5*u_in[0] * (pow(v_in[0], 2.0) + pow(v_in[1], 2.0)));
		double c_in = sqrt(Gamma*p_in / u_in[0]);

		double J_plus = dot_prod(v_in, b_n) + 2 * c_in / (Gamma - 1);
		
		double S_plus = p_in / pow(u_in[0], Gamma);
		double rho_b = pow(p_b / S_plus, 1.0 / Gamma);
		
		double c_b = sqrt(Gamma*p_b / rho_b);
		double u_b = J_plus - 2 * c_b / (Gamma - 1);

		vector<double> temp = vec_add(v_in, scal_multip(b_n, (-1) * dot_prod(v_in, b_n)));
		vector<double> v_b = vec_add(temp, scal_multip(b_n, u_b));

		double rhoE = p_b / (Gamma - 1) + 0.5 * rho_b * (pow(v_b[0],2.0)+pow(v_b[1],2.0));

		vector<double> u_bdry = { rho_b, rho_b*v_b[0], rho_b*v_b[1], rhoE};
		vector<vector<double>> b_flux = state2flux(u_bdry,Gamma);


		flux_out[0] = dot_prod(b_flux[0], b_n);  flux_out[1] = dot_prod(b_flux[1], b_n);  flux_out[2] = dot_prod(b_flux[2], b_n);  flux_out[3] = dot_prod(b_flux[3], b_n);
		vector<double> wavespeed = { abs(dot_prod(v_in,b_n)) + c_in, abs(dot_prod(v_b,b_n)) + c_b };
		flux_out[4] = find_max(wavespeed,false); // wave speed
	}
	
	return flux_out;
}


vector<double> solve_quadratic(double a, double b, double c) {
	if (pow(b, 2.0) < 4 * a*c) {
		//mexPrintf("Complex Quartic Condition Detected. Check the input.");
		throw std::invalid_argument("Complex Quartic Condition Detected. Check the input.");
	}

	double temp1 = sqrt(pow(b, 2.0) - 4 * a*c);
	vector<double> res = { ((-1.0)*b + temp1) / (2*a), ((-1.0)*b - temp1) / (2*a) };

	return res;
}


double find_max(vector<double> vec, bool absmode) {
	
	double max_val;
	if (absmode) {
		max_val = abs(vec[0]);
		for (int ind = 1; ind < vec.size(); ind++) {
			if (abs(vec[ind]) > max_val) {
				max_val = abs(vec[ind]);
			}
		}
	}
	else {
		max_val = vec[0];
		for (int ind = 1; ind < vec.size(); ind++) {
			if (vec[ind] > max_val) {
				max_val = vec[ind];
			}
		}
	}
	return max_val;
}

vector<double> time_step_comp(vector<vector<double>> wavespeed, vector<vector<int>> I2E, vector<vector<int>> B2E, vector<double> Ilength, vector<double> Blength, vector<double> perimeter, vector<double> area, double CFL) {
	vector<double> si(area.size());
	vector<double> di(area.size());
	vector<double> dt(area.size());

	//initialize
	for (int idx = 0;idx < area.size();idx++) {
		si[idx] = 0;
	}

	for (int idx = 0;idx < I2E.size();idx++) {
		double temp = wavespeed[I2E[idx][0]-1][I2E[idx][1]-1]*Ilength[idx];
		si[I2E[idx][0] - 1] = si[I2E[idx][0] - 1] + temp;
		si[I2E[idx][2] - 1] = si[I2E[idx][2] - 1] + temp;
	}
	for (int idx = 0; idx < B2E.size();idx++) {
		double temp = wavespeed[B2E[idx][0]-1][B2E[idx][1]-1] * Blength[idx];
		si[B2E[idx][0] - 1] = si[B2E[idx][0] - 1] + temp;
	}
	for (int idx = 0;idx < area.size();idx++) {
		si[idx] = si[idx] / perimeter[idx];
		di[idx] = 2.0 * area[idx] / perimeter[idx];
		dt[idx] = CFL*di[idx] / si[idx];
	}
	return dt;
}


double find_max_res(vector<vector<double>> Residual) {
	vector<double> max_res_line(Residual.size());
	for (int idx = 0;idx < Residual.size();idx++) {
		max_res_line[idx] = find_max(Residual[idx], true);
	}
	double max_res = find_max(max_res_line, false);
	return max_res;
}

Residual_builder_results residule_builder(vector<vector<double>> u_step, vector<vector<int>> I2E,
	vector<vector<int>> B2E, vector<vector<double>> In, vector<vector<double>> Bn,
	vector<double> Ilength, vector<double> Blength, int Order, double Gamma, double Gas_Const_R,
	vector<vector<double>> Centroid, vector<vector<double>> ICenter, vector<vector<double>> BCenter, 
	vector<double> area, vector<vector<int>> Assem, vector<vector<double>> Coord) {

	int nElem = (int)u_step.size();
	int b_type;
	Residual_builder_results res;

	vector<vector<double>> Residual(nElem, vector<double>(4));
	vector<vector<double>> wavespeed(nElem, vector<double>(3));
	vector<double> local_flux;
	for (int ind = 0; ind < nElem; ind++) {
		Residual[ind][0] = 0; Residual[ind][1] = 0; Residual[ind][2] = 0; Residual[ind][3] = 0;
	}

	Order = 2;
	
	// compute gradients. For boundary elements, zero gradients are implemented.
	vector<vector<vector<double>>> Element_Gradient = gradient(u_step, I2E, B2E, In, Bn, Ilength, Blength, area);
	vector<double> Phi_limit = BJ_scal_limiter(u_step, I2E, B2E, Centroid, Element_Gradient, Assem, Coord);

	// Add the limiting coefficients to the state computations

	//mexPrintf("min Phi: %g\n", *min_element(Phi_limit.begin(), Phi_limit.end()));
	// Loop over interior faces:
	vector<double> u_eL(4); // state on the left of the face
	vector<double> u_eR(4); // state on the right of the face
	vector<double> dx(2);
	vector<double> dU(4);
	for (int idx = 0; idx < I2E.size(); idx++) {
		dx = vec_add(ICenter[idx], scal_multip(Centroid[I2E[idx][0] - 1], -1));
		dU = { dot_prod(Element_Gradient[I2E[idx][0] - 1][0], dx),
			dot_prod(Element_Gradient[I2E[idx][0] - 1][1], dx),
			dot_prod(Element_Gradient[I2E[idx][0] - 1][2], dx),
			dot_prod(Element_Gradient[I2E[idx][0] - 1][3], dx) };
		u_eL = vec_add(u_step[I2E[idx][0] - 1], scal_multip(dU, Phi_limit[I2E[idx][0] - 1]));

		dx = vec_add(ICenter[idx], scal_multip(Centroid[I2E[idx][2] - 1], -1));
		dU = { dot_prod(Element_Gradient[I2E[idx][2] - 1][0], dx),
			dot_prod(Element_Gradient[I2E[idx][2] - 1][1], dx),
			dot_prod(Element_Gradient[I2E[idx][2] - 1][2], dx),
			dot_prod(Element_Gradient[I2E[idx][2] - 1][3], dx) };
		u_eR = vec_add(u_step[I2E[idx][2] - 1], scal_multip(dU, Phi_limit[I2E[idx][2] - 1]));

		//mexPrintf("Original State: \nu_L[3]: %lf u_R[3]: %lf\n", u_step[I2E[idx][0] - 1][3], u_step[I2E[idx][2] - 1][3]);
		//mexPrintf("Corrected State: \nu_L[3]: %lf u_R[3]: %lf\n\n", u_eL[3], u_eR[3]);
		local_flux = Roe_Flux(u_eL, u_eR, In[idx], Gamma, Gas_Const_R);
		wavespeed[I2E[idx][0] - 1][I2E[idx][1] - 1] = local_flux[4];
		wavespeed[I2E[idx][2] - 1][I2E[idx][3] - 1] = local_flux[4];
		local_flux.pop_back();
		Residual[I2E[idx][0] - 1] = vec_add(Residual[I2E[idx][0] - 1], scal_multip(local_flux, Ilength[idx])); //left
		Residual[I2E[idx][2] - 1] = vec_add(Residual[I2E[idx][2] - 1], scal_multip(local_flux, (-1)*Ilength[idx])); //right
	}

	// Loop over boundary faces:

	// Free_stream_test
	/*
	b_type = 1; // Full State
	for (int idx = 0; idx < B2E.size(); idx++) {
		local_flux = Boundary_Flux(u_step[B2E[idx][0] - 1], Bn[idx], b_type, Gamma, Gas_Const_R);
		wavespeed[B2E[idx][0] - 1][B2E[idx][1] - 1] = local_flux[4];
		local_flux.pop_back();
		Residual[B2E[idx][0] - 1] = vec_add(Residual[B2E[idx][0] - 1], scal_multip(local_flux, Blength[idx]));

		//mexPrintf("Boundary face info:\n B2E info: %d %d %d\n Edge Length: %lf\n normal vector: %lf %lf\n\n", B2E[idx][0], B2E[idx][1], B2E[idx][2], Blength[idx], Bn[idx][0], Bn[idx][1]);
		//mexPrintf("local Flux: %lf %lf %lf %lf\n\n", local_flux[0], local_flux[1], local_flux[2], local_flux[3]);
	}
	//mexPrintf("First Elemt Residual: %e %e %e %e\n", Residual[0][0], Residual[0][1], Residual[0][2], Residual[0][3]);
	*/
		
	// True Boundary conditions
		
	for (int idx = 0; idx < B2E.size(); idx++) {

		switch (B2E[idx][2]) {
		case 1: b_type = 2; //bottom, invisid surface
			break;
		case 2: b_type = 4; //right, subsonic outflow
			break;
		case 3: b_type = 2; //top, invisid surface
			break;
		case 4: b_type = 3; //left, subsonic inflow
			break;
		}

		local_flux = Boundary_Flux(u_step[B2E[idx][0] - 1], Bn[idx], b_type, Gamma, Gas_Const_R);
		wavespeed[B2E[idx][0] - 1][B2E[idx][1] - 1] = local_flux[4];
		local_flux.pop_back();
		Residual[B2E[idx][0] - 1] = vec_add(Residual[B2E[idx][0] - 1], scal_multip(local_flux, Blength[idx]));
	}
		
		
	

	res.Residual = Residual;
	res.wavespeed = wavespeed;

	return res;
}


vector<vector<vector<double>>> gradient(vector<vector<double>> &u_step, vector<vector<int>> &I2E, vector<vector<int>> &B2E, 
	vector<vector<double>> &In, vector<vector<double>> &Bn, vector<double> &Ilength, vector<double> &Blength, vector<double> &area) {

	vector<vector<vector<double>>> Grad_Elem(u_step.size(), vector<vector<double>>(4, vector<double>(2))); // first index: element number; second index: state component; third index: x/y components
	// Initialize:
	for (int idx = 0; idx < u_step.size(); idx++) {
		for (int idx2 = 0; idx2 < 4; idx2++) {
			Grad_Elem[idx][idx2] = {0.0, 0.0};
		}
	}
																											
	// loop over interior faces

	for (int idx = 0; idx < I2E.size(); idx++) {
		vector<double> average_u = scal_multip(vec_add(u_step[I2E[idx][0]-1], u_step[I2E[idx][2] - 1]) ,0.5);
		for (int idx2 = 0; idx2 < 4; idx2++) {
			Grad_Elem[I2E[idx][0] - 1][idx2] = vec_add(Grad_Elem[I2E[idx][0] - 1][idx2], scal_multip(In[idx], average_u[idx2]*Ilength[idx]/area[I2E[idx][0] - 1]));
			Grad_Elem[I2E[idx][2] - 1][idx2] = vec_add(Grad_Elem[I2E[idx][2] - 1][idx2], scal_multip(In[idx], -average_u[idx2] * Ilength[idx] / area[I2E[idx][2] - 1]));
		}
	}
	/*
	mexPrintf("Interior:\n");
	for (int idx = 0; idx < u_step.size(); idx++) {
		mexPrintf("Element %d Gradient:\n", idx);
		for (int idx2 = 0; idx2 < 4; idx2++) {
			mexPrintf("%lf  %lf \n", Grad_Elem[idx][idx2][0], Grad_Elem[idx][idx2][1]);
		}
	}
	*/

	// loop over boundary faces. Set the boundary value to the same as the interior cell value

	for (int idx = 0; idx < B2E.size(); idx++) {
		
		for (int idx2 = 0; idx2 < 4; idx2++) {
			Grad_Elem[B2E[idx][0] - 1][idx2] = vec_add(Grad_Elem[B2E[idx][0] - 1][idx2], scal_multip(Bn[idx], u_step[B2E[idx][0]-1][idx2] * Blength[idx] / area[B2E[idx][0] - 1]));
		}
	}

	return Grad_Elem;

}

