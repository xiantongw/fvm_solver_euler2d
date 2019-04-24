# include "Limiter_header.h"


vector<double> BJ_scal_limiter(
	vector<vector<double>> &u_step,
	vector<vector<int>> &I2E,
	vector<vector<int>> &B2E,
	vector<vector<double>> &Centroid,
	vector<vector<vector<double>>> &Gradient,
	vector<vector<int>> &Assem,
	vector<vector<double>> &Coord)
{
	int nElem = (int)Assem.size();
	vector<double> Phi_lim(nElem, 1.0);
	
	// conduct limit on density
	vector<double> umax(nElem);
	vector<double> umin(nElem);
	for (int idx = 0;idx < nElem;idx++) {
		umax[idx] = u_step[idx][0];
		umin[idx] = u_step[idx][0];
	}

	// loop through interior faces;
	int nIface = (int)I2E.size();
	vector<double> states(4);
	for (int idx = 0;idx < nIface;idx++) {
		states[0] = u_step[I2E[idx][0]-1][0]; states[1] = u_step[I2E[idx][0] - 1][0];
		states[2] = umax[I2E[idx][0] - 1];   states[3] = umin[I2E[idx][0] - 1];
		umax[I2E[idx][0] - 1] = *max_element(states.begin(), states.end());
		umin[I2E[idx][0] - 1] = *min_element(states.begin(), states.end());

		states[2] = umax[I2E[idx][2] - 1]; states[3] = umin[I2E[idx][2] - 1];
		umax[I2E[idx][2] - 1] = *max_element(states.begin(), states.end());
		umin[I2E[idx][2] - 1] = *min_element(states.begin(), states.end());
	}
	
	// Loop through elements to obtain nodal values and hence determine the Phi_value
	// can be parallelized
	vector<double> node_coord(2);
	vector<double> dx(2);
	vector<double> phi_temp(3);
	double u_temp;
	
	for (int idx = 0;idx < nElem;idx++) {
		fill(phi_temp.begin(), phi_temp.end(), 1);
		for (int j = 0;j < 3;j++) {
			node_coord = Coord[Assem[idx][j] - 1];
			transform(node_coord.begin(), node_coord.end(), Centroid[idx].begin(), dx.begin(), std::minus<double>());
			u_temp = u_step[idx][0] + dot_prod(Gradient[idx][0], dx);
			
			if (u_temp > u_step[idx][0]) {
				phi_temp[j] = min(1.0, (umax[idx] - u_step[idx][0]) / (u_temp - u_step[idx][0]));
			}
			else if (u_temp < u_step[idx][0]) {
				phi_temp[j] = min(1.0, (umin[idx] - u_step[idx][0]) / (u_temp - u_step[idx][0]));
			}
			
		}
		Phi_lim[idx] = *min_element(phi_temp.begin(), phi_temp.end());
	}
	
	return Phi_lim;
}

/*
vector<double> VL_scal_limiter(vector<vector<int>> &I2E,
	vector<vector<vector<double>>> &Gradient,
	vector<vector<double>> &In)

	*/