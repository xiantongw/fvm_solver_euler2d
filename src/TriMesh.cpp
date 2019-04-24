#include "../include/TriMesh.h"

TriMesh::TriMesh(string &gri_filename_in)
{
    gri_filename = gri_filename_in;
    ReadGri(this->gri_filename);
    this->isCurved = vector<bool> (this->E.size(), false);
    this->curved_group = -1;
    CalcArea();
    CalcCentroid();
    CalcI2E();
    CalcB2E();
    CalcIn();
    CalcBn();
    FindNeighbour();
    FindCurvedIndex();
}

TriMesh::TriMesh(TriMesh &mesh)
{
    gri_filename = mesh.gri_filename;
    num_node = mesh.num_node;
    num_element = mesh.num_element;
    num_boundary = mesh.num_boundary;
    E = mesh.E;
    V = mesh.V;
    B = mesh.B;
    Bname = mesh.Bname;
    name_base = mesh.name_base;
    n_base = mesh.n_base;
    I2E = mesh.I2E;
    B2E = mesh.B2E;
    In = mesh.In;
    Bn = mesh.Bn;
    Area = mesh.Area;
    Centroid = mesh.Centroid;
    isCurved = mesh.isCurved;
    CurvedElementIndex = mesh.CurvedElementIndex;
    CurvedEdgeIndex = mesh.CurvedEdgeIndex;
    LinearElementIndex = mesh.LinearElementIndex;
    LinearEdgeIndex = mesh.LinearEdgeIndex;
    curved_group = mesh.curved_group;
    NeighbourElements = mesh.NeighbourElements;
}

void TriMesh::ReadGri(string &gri_filename)
{

    // variables dealing with strings
    ifstream gri_file(gri_filename);
    stringstream ss;
    string line;
    // Variables
    int num_node, num_element, dim, num_boundary, n_base;
    // Read the first line, get Nn, Ne, dim
    getline(gri_file, line);
    ss.clear();
    ss.str(line);
    ss >> num_node >> num_element >> dim;

    // Read the V matrix, which are the coordinates of the nodes
    vector<vector<double> > V(num_node, vector<double>(2));
    for (int i = 0; i < num_node; i++)
    {
        getline(gri_file, line);
        ss.clear();
        ss.str(line);
        ss >> V[i][0] >> V[i][1];
    }

    // Read number of boundaries
    getline(gri_file, line);
    ss.clear();
    ss.str(line);
    ss >> num_boundary;

    // Read boundary information
    vector<string> Bname(num_boundary);
    vector<vector<vector<int> > > B;
    B.resize(num_boundary);
    for (int i_boundary = 0; i_boundary < num_boundary; i_boundary++)
    {
        // The information of a specific boundary
        int num_boundary_edge, dim;
        getline(gri_file, line);
        ss.clear();
        ss.str(line);
        ss >> num_boundary_edge >> dim >> Bname[i_boundary];
        B[i_boundary].resize(num_boundary_edge);
        for (int i = 0; i < num_boundary_edge; i++)
        {
            getline(gri_file, line);
            ss.clear();
            ss.str(line);
            B[i_boundary][i].resize(2);
            ss >> B[i_boundary][i][0] >> B[i_boundary][i][1];
        }
    }

    // Read element information
    vector<vector<int> > E(num_element);
    vector<string> name_base;
    vector<int> n_bases;

    int current_total_elements = 0, num_element_part, num_element_part_old = 0;
    int nnode, buffer_int;
    string buffer_str;

    while (current_total_elements != num_element)
    {

        getline(gri_file, line);
        ss.clear();
        ss.str(line);
        ss >> num_element_part >> n_base >> buffer_str;
        n_bases.push_back(n_base);
        name_base.push_back(buffer_str);

        for (int i_element = 0; i_element < num_element_part; i_element++)
        {
            getline(gri_file, line);
            ss.clear();
            ss.str(line);
            nnode = (n_base + 1) * (n_base + 2) / 2;

            E[i_element + num_element_part_old].resize(nnode);
            int i_local_node = 0;
            while (ss >> buffer_int)
            {
                E[i_element + num_element_part_old][i_local_node] = buffer_int;
                i_local_node++;
            }
        }
        current_total_elements += num_element_part;
        num_element_part_old = num_element_part;
    }
    gri_file.close();

    this->num_boundary = num_boundary;
    this->num_element = num_element;
    this->num_node = num_node;
    this->name_base = name_base;
    this->n_base = n_bases;
    this->B = B;
    this->Bname = Bname;
    this->E = E;
    this->V = V;
}

void TriMesh::CalcI2E()
{
    // This function is re-written from the Python version, so the names of the
    // variables are not well-organized
    vector<vector<int> > E = this->E;
    vector<vector<double> > V = this->V;
    int N = this->num_node;
    vector<vector<int> > H(N, vector<int>(N));
    vector<vector<int> > C(this->num_element * 3, vector<int>(4, 0));
    int nedge = 0;
    for (int t = 0; t < this->num_element; t++)
    {
        vector<int> vertex_index(3);
        vertex_index = utils::GetVertexIndex(E[t]);
        for (int e = 0; e < 3; e++)
        {
            int n1 = E[t][vertex_index[(e + 1) % 3]] - 1;
            int n2 = E[t][vertex_index[(e + 2) % 3]] - 1;
            int nmin = min(n1, n2);
            int nmax = max(n1, n2);
            if (H[nmin][nmax] > 0)
            {
                int tN = H[nmin][nmax];
                int eN = H[nmax][nmin];
                int t1 = t + 1;
                int t2 = tN;
                int e1 = e + 1;
                int e2 = eN;
                if (t2 < t1)
                {
                    t1 = tN;
                    t2 = t + 1;
                    e1 = eN;
                    e2 = e + 1;
                }
                C[nedge][0] = t1;
                C[nedge][1] = e1;
                C[nedge][2] = t2;
                C[nedge][3] = e2;
                nedge = nedge + 1;
            }
            else
            {
                H[nmin][nmax] = t + 1;
                H[nmax][nmin] = e + 1;
            }
        }
    }
    // Sorting
    vector<vector<int> > CC(nedge, vector<int>(4, 0));
    for (int i = 0; i < nedge; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            CC[i][j] = C[i][j];
        }
    }
    sort(CC.begin(), CC.end(), utils::SortByColumn0);
    int i0 = 0;
    for (int i = 0; i < nedge; i++)
    {
        if (CC[i][0] != CC[i0][0])
        {
            vector<vector<int> > A = utils::SliceByRow(CC, i0, i - 1);
            sort(A.begin(), A.end(), utils::SortByColumn2);
            for (int ii = i0; ii < i; ii++)
            {
                for (int jj = 0; jj < 4; jj++)
                {
                    CC[ii][jj] = A[ii - i0][jj];
                }
            }
            i0 = i;
        }
    }
    this->I2E = CC;
}

void TriMesh::CalcB2E()
{
    // This function is re-written from the Python version, so the names of the
    // variables are not well-organized
    vector<vector<int> > E = this->E;
    vector<vector<vector<int> > > B = this->B;
    vector<vector<int> > B2E;

    for (int ielem = 0; ielem < E.size(); ielem++)
    {
        vector<int> vertex_index(3), elem(3);
        vertex_index = utils::GetVertexIndex(E[ielem]);
        for (int i = 0; i < 3; i++)
        {
            elem[i] = E[ielem][i];
        }
        for (int ibg = 0; ibg < this->Bname.size(); ibg++)
        {
            vector<vector<int> > B_part = B[ibg];
            for (int ib = 0; ib < B_part.size(); ib++)
            {
                vector<int> nb = B_part[ib];
                int match_1 = distance(elem.begin(), find(elem.begin(), elem.end(), nb[0]));
                int match_2 = distance(elem.begin(), find(elem.begin(), elem.end(), nb[1]));
                if ((match_1 != elem.size()) && (match_2 != elem.size()))
                {
                    int local_ind = utils::MissingFrom012(match_1, match_2);
                    vector<int> row = {ielem + 1, local_ind + 1, ibg + 1};
                    B2E.push_back(row);
                }
            }
        }
    }
    sort(B2E.begin(), B2E.end(), utils::SortByColumn2);
    this->B2E = B2E;
}

void TriMesh::FindCurvedIndex()
{
    this->CurvedElementIndex.clear();
    this->LinearElementIndex.clear();
    this->CurvedEdgeIndex.clear();
    this->LinearEdgeIndex.clear();
    for (int i = 0; i < this->isCurved.size(); i++)
    {
        if (this->isCurved[i])
        {
            this->CurvedElementIndex.push_back(i);
        }
        else
        {
            this->LinearElementIndex.push_back(i);
        }
    }
    for (int i = 0; i < this->Bn.size(); i++)
    {
        if (this->Bn[i][3] > 0)
        {
            this->CurvedEdgeIndex.push_back(i);
        }
        else
        {
            this->LinearEdgeIndex.push_back(i);
        }
    }
}

void TriMesh::CalcIn()
{
    int num_edge = this->I2E.size();
    vector<vector<double> > In(num_edge,  vector<double> (3));
    for (int iedge = 0; iedge < num_edge; iedge++)
    {
        int ielemL = this->I2E[iedge][0] - 1;
        int ilocL = this->I2E[iedge][1] - 1;
        int ilocA = (ilocL + 1) % 3;
        int ilocB = (ilocL + 2) % 3;
        int iglobA = this->E[ielemL][ilocA] - 1;
        int iglobB = this->E[ielemL][ilocB] - 1;
        double xA = this->V[iglobA][0]; double yA = this->V[iglobA][1];
        double xB = this->V[iglobB][0]; double yB = this->V[iglobB][1];
        double dl = sqrt((xA - xB) * (xA - xB) + (yA - yB) * (yA - yB));
        In[iedge][2] = dl;
        In[iedge][0] = (yB - yA) / dl; In[iedge][1] = (xA - xB) / dl;
    }
    this->In = In;
}

void TriMesh::CalcBn()
{
    int num_bedge = this->B2E.size();
    vector<vector<double> > Bn(num_bedge, vector<double> (4));
    this->Bn.clear();
    for (int iedge = 0; iedge < num_bedge; iedge++)
    {
        int ielem = this->B2E[iedge][0] - 1;
        int iloc = this->B2E[iedge][1] - 1;
        vector<int> ind_vertex = utils::GetVertexIndex(this->E[ielem]);
        int ilocA = (iloc + 1) % 3;
        int ilocB = (iloc + 2) % 3;
        int iglobA = this->E[ielem][ind_vertex[ilocA]] - 1;
        int iglobB = this->E[ielem][ind_vertex[ilocB]] - 1;
        double xA = this->V[iglobA][0]; double yA = this->V[iglobA][1];
        double xB = this->V[iglobB][0]; double yB = this->V[iglobB][1];
        double dl = sqrt((xA - xB) * (xA - xB) + (yA - yB) * (yA - yB));
        Bn[iedge][2] = dl;
        Bn[iedge][0] = (yB - yA) / dl; Bn[iedge][1] = (xA - xB) / dl;
        if(this->B2E[iedge][2] == this->curved_group)
        {
            Bn[iedge][3] = 1;
        } else
        {
            Bn[iedge][3] = 0;
        }
    }
    this->Bn = Bn;
}

void TriMesh::CalcArea()
{
    vector<double> Area(this->num_element);
    double xA, yA, xB, yB, xC, yC, area;
    for (int i = 0; i < this->num_element; i++)
    {
        xA = this->V[this->E[i][0] - 1][0]; yA = this->V[this->E[i][0] - 1][1];
        xB = this->V[this->E[i][1] - 1][0]; yB = this->V[this->E[i][1] - 1][1];
        xC = this->V[this->E[i][2] - 1][0]; yC = this->V[this->E[i][2] - 1][1];
        area = abs(0.5 * (xA * (yB - yC) + xB * (yC - yA) + xC * (yA - yB)));
        Area[i] = area;
    }
    this->Area = Area;
}

void TriMesh::CalcCentroid()
{
    vector<vector<double> > Centroid(this->num_element, vector<double> (2, 0.0));
    double xA, yA, xB, yB, xC, yC;
    for (int i = 0; i < this->num_element; i++)
    {
        xA = this->V[this->E[i][0] - 1][0]; yA = this->V[this->E[i][0] - 1][1];
        xB = this->V[this->E[i][1] - 1][0]; yB = this->V[this->E[i][1] - 1][1];
        xC = this->V[this->E[i][2] - 1][0]; yC = this->V[this->E[i][2] - 1][1];
        Centroid[i][0] = (xA + xB + xC) / 3.0;
        Centroid[i][1] = (yA + yB + yC) / 3.0;
    }
    this->Centroid = Centroid;
}

void TriMesh::FindNeighbour()
{
    vector<vector<int> > NeighbourElements(this->num_element, vector<int> (3, -1));
    vector<int> current_ind(this->num_element, 0);
    for (int i = 0; i < this->I2E.size(); i++)
    {
        int elemL = this->I2E[i][0];
        int elemR = this->I2E[i][2];
        NeighbourElements[elemL - 1][current_ind[elemL - 1]] = elemR;
        NeighbourElements[elemR - 1][current_ind[elemR - 1]] = elemL;
        current_ind[elemL - 1]++;
        current_ind[elemR - 1]++;
    }
    for (int i = 0; i < this->B2E.size(); i++)
    {
        int elemL = this->B2E[i][0];
        NeighbourElements[elemL - 1][current_ind[elemL - 1]] = elemL;
        current_ind[elemL - 1]++;
    }
    this->NeighbourElements = NeighbourElements;
}

void TriMesh::WriteGri(string& gri_filename)
{
    ofstream outfile;
    outfile.open(gri_filename);
    outfile << this->V.size() << ' ' << this->E.size() << ' ' << this->V[0].size() << endl;
    for (int i = 0; i < this->V.size(); i++)
    {
        outfile << std::scientific << this->V[i][0] << ' ' << this->V[i][1] << endl;
    }
    outfile << this->Bname.size() << endl;
    for (int iB = 0; iB < this->Bname.size(); iB++)
    {
        outfile << this->B[iB].size() << ' ' << 2 << ' ' << this->Bname[iB] << endl;
        for (int j = 0; j < this->B[iB].size(); j++)
        {
            outfile << this->B[iB][j][0] << ' ' << this->B[iB][j][1] << endl;
        }
    }
    // Write our the curve and uncurved elements
    for (int iE = 0; iE < this->E.size(); iE++)
    {
        int nnode = this->E[iE].size();
        int p = int((sqrt(8 * nnode + 1) - 3) / 2);
        outfile << 1 << ' ' << p << ' ' << "TriLagrange" << endl;
        for (int inode = 0; inode < nnode; inode++)
        {
            outfile << this->E[iE][inode] << ' ';
        }
        outfile << endl;

    }
    outfile.close();
}