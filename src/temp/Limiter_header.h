# include "mex.h"
# include "FVM_header.h"
using namespace std;

vector<double> BJ_scal_limiter(
	vector<vector<double>> &u_step,
	vector<vector<int>> &I2E,
	vector<vector<int>> &B2E,
	vector<vector<double>> &Centroid,
	vector<vector<vector<double>>> &Gradient,
	vector<vector<int>> &Assem,
	vector<vector<double>> &Coord);

vector<double> VL_scal_limiter(vector<vector<int>> &I2E,
	vector<vector<vector<double>>> &Gradient,
	vector<vector<double>> &In);