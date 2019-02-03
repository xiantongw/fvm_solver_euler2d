#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 16:47:04 2019

@author: xtwang
"""


class TriMesh:
    # variables
    E = None
    V = None
    B = None
    Bname = None
    I2E = None
    B2E = None
    In = None
    Bn = None
    Area = None
    MeshName = None
    
    # methods
    def __init__(self, filename):
        from readgri import readgri
        import numpy as np
        m = readgri(filename)
        self.E = m['E'] + 1
        self.V = m['V']
        self.B = m['B']
        self.Bname = m['Bname']
        for ibg in range(np.shape(m['B'])[0]):
            self.B[ibg] = self.B[ibg] + 1
        self.I2E = self.cal_I2E()
        self.B2E = self.cal_B2E()
        self.In = self.cal_In()
        self.Bn = self.cal_Bn()
        self.Area = self.cal_Area()
        self.MeshName = filename

        
    def info(self):
        import numpy as np
        print('================== Mesh Infomation ====================')
        print('Number of Elements: {:d}'.format(np.shape(self.E)[0]))
        print('Number of Nodes: {:d}'.format(np.shape(self.V)[0]))
        print('Number of Interior Edges: {:d}'.format(np.shape(self.I2E)[0]))
        print('Number of Boundary Edges: {:d}'.format(np.shape(self.B2E)[0]))
        print('=======================================================')
        
    def save_file(self, filename):
        import numpy as np
        E = self.E
        V = self.V
        B = self.B
        Bname = self.Bname
        f = open(filename,"w")
        f.write("{:d} {:d} {:d}\n".format(np.shape(V)[0], np.shape(E)[0], np.shape(V)[1]))
        for i in range(np.shape(V)[0]):
            f.write("{:.8E} {:.8E}\n".format(V[i, 0], V[i, 1]))
        f.write("{:d}\n".format(len(Bname)))
        for iB in range(len(Bname)):
            # write the bounday info (1 line)
            f.write("{:d} {:d} {:s}\n".format(len(B[iB]), 2, Bname[iB]))
            for j in range(len(B[iB])):
                f.write("{:d} {:d}\n".format(B[iB][j, 0], B[iB][j, 1]))
        f.write("{:d} {:d} {:s}\n".format(len(E), 1, "TriLagrange"))
        for i in range(len(E)):
            f.write("{:d} {:d} {:d}\n".format(E[i, 0], E[i, 1], E[i, 2]))
        f.close()
        print("Writing to File...Done Filename: "+filename)

    def cal_I2E(self):
        import numpy as np
        from scipy import sparse
        # make a hash table to identify  edges
        E = self.E
        V = self.V
        N = V.shape[0]
        H = sparse.lil_matrix((N,N))
        nelem = E.shape[0]
        C = np.zeros([nelem*3,4], dtype=np.int)
        nedge = 0
        for t in range(nelem):
            for e in range(3):
                n1 = E[t,(e + 1) % 3] - 1; n2 = E[t,(e + 2) % 3] - 1
                nmin = min(n1,n2); nmax = max(n1,n2)
                if(H[nmin, nmax] > 0):
                    tN = H[nmin, nmax]; eN = H[nmax, nmin]
                    t1 = t+1; t2 = tN
                    e1 = e+1; e2 = eN
                    if(t2 < t1):
                        t1 = tN; t2 = t + 1
                        e1 = eN; e2 = e + 1
                    C[nedge,:]  = [t1, e1, t2, e2]
                    nedge = nedge + 1
                else:
                    H[nmin, nmax] = t + 1
                    H[nmax, nmin] = e + 1
        # sorting
        CC = np.zeros([nedge,4], dtype=np.int)
        CC[:,:]  = C[0:nedge, :]
        CC = CC[CC[:, 0].argsort()]
        i0 = 0
        for i in range(nedge):
            if(CC[i, 0] != CC[i0, 0]):
                A = CC[i0:i, :]
                A = A[A[:, 2].argsort()]
                CC[i0:i, :] = A
                i0 = i
        return CC   
    
    def cal_B2E(self):
        import numpy as np
        E = self.E
        B = self.B
        nbgroup = np.shape(self.Bname)[0]
        B2E = np.empty((0, 3), dtype=int)
        for ibg in range(nbgroup):
            B_part = B[ibg]
            for ib in range(np.shape(B_part)[0]):
                nb = B_part[ib, :]
                for ielem in range(np.shape(E)[0]):
                    elem = E[ielem, :]
                    match_1 = np.where(nb[0] == elem)[0]
                    match_2 = np.where(nb[1] == elem)[0]
                    if match_1.size > 0 and match_2.size > 0:
                        local_ind_set = {0, 1, 2} - {match_1[0], match_2[0]}
                        for value in local_ind_set:
                            local_ind = value
                        B2E = np.vstack([B2E, [ielem + 1, local_ind + 1, ibg + 1]])
                        break
        return B2E
    
    def cal_In(self):
        import numpy as np
        E = self.E
        V = self.V
        I2E = self.I2E        
        nedge = np.shape(I2E)[0]
        In = np.zeros([nedge, 2])
        for iedge in range(nedge):
            iL = I2E[iedge, 0] - 1
            local_ind = I2E[iedge, 1] - 1
            pA = V[E[iL, (local_ind + 1) % 3] - 1]
            pB = V[E[iL, (local_ind + 2) % 3] - 1]        
            dl = np.linalg.norm(pA - pB)
            In[iedge, 0] = (pB[1] - pA[1]) / dl
            In[iedge, 1] = (pA[0] - pB[0]) / dl
        return In
    
    def cal_Bn(self):
        import numpy as np
        B2E = self.B2E
        E = self.E
        V = self.V
        Bn = np.zeros([np.shape(B2E)[0], 2], dtype=np.float)
        for ib in range(np.shape(B2E)[0]):
            iE = B2E[ib, 0] - 1
            iLocal = B2E[ib, 1] - 1
            pA = V[E[iE, (iLocal + 1) % 3] - 1]
            pB = V[E[iE, (iLocal + 2) % 3] - 1]
            dl = np.linalg.norm(pA - pB)
            Bn[ib, 0] = (pB[1] - pA[1]) / dl
            Bn[ib, 1] = (pA[0] - pB[0]) / dl
        return Bn
    
    def cal_Area(self):
        import numpy as np
        nelem = np.shape(self.E)[0]
        A = np.zeros([nelem, 1])
        for iTriangle in range(nelem):
            pA = self.V[self.E[iTriangle, 0] - 1]
            pB = self.V[self.E[iTriangle, 1] - 1]
            pC = self.V[self.E[iTriangle, 2] - 1]
            area = np.abs(0.5 * (pA[0]*(pB[1] - pC[1]) + pB[0]*(pC[1] - pA[1])
                                    + pC[0]*(pA[1] - pB[1])))
            A[iTriangle] = area
        return A  
    
    def verification(self):
        import numpy as np
        E,V,I2E,B2E,In,Bn = self.E,self.V,self.I2E,self.B2E,self.In,self.Bn
        
        vec_sum = np.zeros([np.shape(E)[0], 3])
        
        for iedge in range(np.shape(I2E)[0]):
            [elemL, faceL, elemR, faceR] = I2E[iedge, :]
            nv = In[iedge, :]
            local_node = list({0, 1, 2} - {faceL - 1})
            iA = E[elemL - 1, local_node[0]] - 1
            iB = E[elemL - 1, local_node[1]] - 1
            dl = np.linalg.norm(V[iA, :] - V[iB, :])
            vec_sum[elemL - 1, 2] += 1
            vec_sum[elemR - 1, 2] += 1
            vec_sum[elemL - 1, 0:2] += nv * dl
            vec_sum[elemR - 1, 0:2] += -1.0 * nv * dl
            
        for bedge in range(np.shape(B2E)[0]):
            [elem, face, bgroup] = B2E[bedge, :]
            nv = Bn[bedge, :]
            local_node = list({0, 1, 2} - {face - 1})
            iA = E[elem - 1, local_node[0]] - 1
            iB = E[elem - 1, local_node[1]] - 1
            dl = np.linalg.norm(V[iA, :] - V[iB, :])
            vec_sum[elem - 1, 2] += 1
            vec_sum[elem - 1, 0:2] += nv * dl    
        
        max_err = np.max(np.sqrt(vec_sum[:,0]**2 + vec_sum[:,1]**2))
        print('================== Mesh Verification ==================')
        print('Maximum magnitude of mesurement vector: {:e}'.format(max_err))
        if max_err < np.finfo(float).eps:
            print('Verification Passed')
            print('=======================================================')
        else:
            print('Verification Failed')
            print('=======================================================')
            
    def refine(self):
        import numpy as np
        from bump_func import bump_func
        import re
        
        E,V,B,I2E,B2E = self.E,self.V,self.B,self.I2E,self.B2E

        nelem = np.shape(E)[0]
        refine_cache = np.zeros([nelem, 3]) - 1
        
        for iedge in range(np.shape(I2E)[0]):
            [elemL, faceL, elemR, faceR] = I2E[iedge, :]
            local_node = list({0, 1, 2} - {faceL - 1})
            iA = E[elemL - 1, local_node[0]] - 1
            iB = E[elemL - 1, local_node[1]] - 1
            # Add mew node to matrix V and label the index of new node to refine_cache
            V = np.vstack([V, (V[iA, :] + V[iB, :])/2])
            i_new_node = np.shape(V)[0]
            refine_cache[elemL - 1, faceL - 1] = i_new_node
            refine_cache[elemR - 1, faceR - 1] = i_new_node
        
        for bedge in range(np.shape(B2E)[0]):
            [elem, face, bgroup] = B2E[bedge, :]
            local_node = list({0, 1, 2} - {face - 1})
            iA = E[elem - 1, local_node[0]] - 1
            iB = E[elem - 1, local_node[1]] - 1
            if bgroup == 1 and re.match('bump*', self.MeshName):
                new_x = (V[iA, 0] + V[iB, 0])/2
                new_y = bump_func(new_x)
                V = np.vstack([V, [new_x, new_y]])
            else:
                V = np.vstack([V, (V[iA, :] + V[iB, :])/2])
            i_new_node = np.shape(V)[0]
            refine_cache[elem - 1, face - 1] = i_new_node
            # modify B matrix
            for ib in range(np.shape(B[bgroup-1])[0]):
                if {B[bgroup-1][ib, 0], B[bgroup-1][ib, 1]} == {iA+1, iB+1}:
                    break
            B[bgroup-1] = np.delete(B[bgroup-1], ib, axis=0)
            B[bgroup-1] = np.vstack([B[bgroup-1], [i_new_node, iA+1]])
            B[bgroup-1] = np.vstack([B[bgroup-1], [i_new_node, iB+1]])    
        
        # generate new E matrix based on information in refine_cache
        new_E = np.empty((0, 3), dtype=int)
        for ielem in range(np.shape(E)[0]):
            new_E = np.vstack([new_E, [E[ielem, 0], 
                                       refine_cache[ielem, 2],refine_cache[ielem, 1]]])
            new_E = np.vstack([new_E, [E[ielem, 1], 
                                       refine_cache[ielem, 0],refine_cache[ielem, 2]]])
            new_E = np.vstack([new_E, [E[ielem, 2], 
                                       refine_cache[ielem, 1],refine_cache[ielem, 0]]])
            new_E = np.vstack([new_E, [refine_cache[ielem, 0], 
                                       refine_cache[ielem, 1],refine_cache[ielem, 2]]])
        
        self.V = V
        self.E = new_E.astype(int)
        self.B = B
        self.I2E = self.cal_I2E()
        self.B2E = self.cal_B2E()
        self.In = self.cal_In()
        self.Bn = self.cal_Bn()
        self.Area = self.cal_Area()