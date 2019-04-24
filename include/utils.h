#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#include "./TriMesh.h"

namespace utils
{

    int GetFullOrderIndex(int r, int s, int order);

    std::vector<int> GetVertexIndex(std::vector<int>& element);

    bool SortByColumn0(std::vector<int> const& v1, std::vector<int> const& v2);

    bool SortByColumn2(std::vector<int> const& v1, std::vector<int> const& v2);

    template <typename T>
    std::vector<std::vector<T> > SliceByRow(std::vector<std::vector<T> > const& v, int m, int n);

    int MissingFrom012(int input1, int input2);

    void mod3(int input, int res[2]);

    template <typename T>
    T norm_inf(std::vector<std::vector<T> > vector_2d);

    std::vector<std::vector<double> > Invert22Matrix(std::vector<std::vector<double> > mat_input);

    std::vector<std::vector<double> > Trans22Matrix(std::vector<std::vector<double> > mat_input);

    double max(double x, double y);
    double min(double x, double y);

    void SaveState(std::vector<std::vector<double> > state_vectors);

}

#endif