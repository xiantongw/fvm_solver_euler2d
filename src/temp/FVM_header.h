// FVM Header
# include <iostream>
# include <cmath>
# include <vector>
# include <stdexcept>
# include <algorithm> 
# include "mex.h"

//#define printfFnc(...) { mexPrintf(__VA_ARGS__); mexEvalString("drawnow;");}
using namespace std;

struct Residual_builder_results
{
	vector<vector<double>> Residual;
	vector<vector<double>> wavespeed;
	char *other;
};



void FVM_solver(
	vector<vector<double>> &U0,
	const vector<vector<int>> &I2E,
	const vector<vector<int>> &B2E,
	const vector<vector<double>> &In,
	const vector<vector<double>> &Bn,
	const vector<double> &area,
	const vector<double> &Ilength,
	const vector<double> &Blength,
	const double CFL,
	const double Rtol,
	const int iterlim,
	int Order,
	const double Gamma,
	const double Gas_Const_R,
	const vector<vector<double>> &Centroid,
	const vector<vector<double>> &ICenter,
	const vector<vector<double>> &BCenter,
	const vector<vector<double>> &Coord,
	const vector<vector<int>> &Assem,
	double U_e[],
	double Res_history[],
	double Res_last_step[]
	);


vector<double> perimeter(const vector<double> Ilength, const vector<double> Blength, const vector<vector<int>> I2E, const vector<vector<int>> B2E, int nelem);

vector<double> Roe_Flux(vector<double> u_L, vector<double> u_R, vector<double> n_elem, double Gamma, double gas_const_R); // the first four terms are the flux, the last term is the max wavespeed
vector<vector<double>> state2flux(vector<double> u, double Gamma);
void Roe_flux_tester(double Gamma);

vector<double> Boundary_Flux(vector<double> u_in, vector<double> b_n, int b_type, double Gamma, double Gas_Const_R);// this function is designed for this problem ONLY. Change parameters inside here for other problems.


vector<double> scal_multip(vector<double> vec, double c);
vector<double> vec_add(vector<double> a, vector<double> b);
double dot_prod(vector<double> a, vector<double> b);
vector <double> solve_quadratic(double a, double b, double c);
double find_max(vector<double> vec, bool absmode);

vector<double> time_step_comp(vector<vector<double>> wavespeed, vector<vector<int>> I2E, vector<vector<int>> B2E, 
	vector<double> Ilength, vector<double>Blength, vector<double> perimeter, vector<double> area, double CFL);

double find_max_res(vector<vector<double>> Residual);

Residual_builder_results residule_builder(vector<vector<double>> u_step, 
	vector<vector<int>> I2E, 
	vector<vector<int>> B2E, 
	vector<vector<double>> In, 
	vector<vector<double>> Bn, 
	vector<double> Ilength, 
	vector<double> Blength, 
	int Order, 
	double Gamma, 
	double Gas_Const_R,
	vector<vector<double>> Centroid,
	vector<vector<double>> ICenter,
	vector<vector<double>> BCenter,
	vector<double> area,
	vector<vector<int>> Assem,
	vector<vector<double>> Coord);

vector<vector<vector<double>>> gradient(vector<vector<double>> &u_step, vector<vector<int>> &I2E, vector<vector<int>> &B2E,
	vector<vector<double>> &In,  vector<vector<double>> &Bn,  vector<double> &Ilength, vector<double> &Blength, vector<double> &area);
