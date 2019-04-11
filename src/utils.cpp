#include "../include/utils.h"

namespace utils
{

    using namespace std;

    int GetFullOrderIndex(int r, int s, int order)
    {
        int sum_temp = 0;
        int k;
        for (int sp = 0; sp < s; sp++)
        {
            sum_temp += order + 1 - sp;
        }
        k = sum_temp + r + 1;
        return k - 1; // return the index which is zero based
    }

    std::vector<int> GetVertexIndex(std::vector<int> &element)
    {
        int num_node_in_element = element.size();
        // calculate the order of bases based on the number of nondes in the element
        int num_order = int((sqrt(8 * num_node_in_element + 1) - 3) / 2);
        std::vector<int> vertex_index(3);
        vertex_index[0] = GetFullOrderIndex(0, 0, num_order);
        vertex_index[1] = GetFullOrderIndex(num_order, 0, num_order);
        vertex_index[2] = GetFullOrderIndex(0, num_order, num_order);
        return vertex_index;
    }

    bool SortByColumn0(vector<int> const &v1, vector<int> const &v2)
    {
        return v1[0] < v2[0];
    }

    bool SortByColumn2(vector<int> const &v1, vector<int> const &v2)
    {
        return v1[2] < v2[2];
    }

    template <typename T>
    vector<vector<T> > SliceByRow(vector<vector<T> > const &v, int m, int n)
    {
        auto first = v.cbegin() + m;
        auto last = v.cbegin() + n + 1;
        vector<vector<T> > vec(first, last);
        return vec;
    }
    template vector<vector<int> > SliceByRow<int>(vector<vector<int> > const &v, int m, int n);

    int MissingFrom012(int input1, int input2)
    {
        if ((input1 == 0 && input2 == 1) || (input1 == 1 && input2 == 0))
        {
            return 2;
        }
        else if ((input1 == 0 && input2 == 2) || (input1 == 2 && input2 == 0))
        {
            return 1;
        }
        else if ((input1 == 1 && input2 == 2) || (input1 == 2 && input2 == 1))
        {
            return 0;
        }
        else
        {
            abort();
        }
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

    template <typename T>
    T norm_inf(std::vector<std::vector<T> > vector_2d)
    {
        int size1 = vector_2d.size();
        int size2 = vector_2d[0].size();
        T norminf = -1.0 * std::numeric_limits<double>::infinity() + 1;
        for (int i = 0; i < size1; i++)
        {
            for (int j = 0; j < size2; j++)
            {
                if (norminf < vector_2d[i][j])
                {
                    norminf = vector_2d[i][j];
                }
            }
        }
        return norminf;
    }
    template double norm_inf(std::vector<std::vector<double> > vector_2d);

} // namespace utils