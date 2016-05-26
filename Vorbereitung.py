import numpy as np
Se = np.array(([1, -1], [-1, 1]))
def stiffness_matrix(e_num, N_num):
    # calculate the dimension of stiffness matrix
    S_num = e_num * N_num - e_num + 1
    S = np.zeros((S_num, S_num))
    for e in range(e_num):
        el = (N_num - 1) * e
        for i in range(N_num):
            for j in range(N_num):
                S[el + i][el + j] += Se[i][j]
    return S



# calculate nodal displacements after LR Zerlegung
from numpy import matrix
from numpy import linalg
L_mat = np.asmatrix(L)
R_mat = np.asmatrix(R)
F_mat = np.asmatrix(F)
y_mat = linalg.solve(L_mat, F_mat)
y_mat
D_mat = linalg.solve(R_mat, y_mat)
D_mat
