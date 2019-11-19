
import pytest, numpy as np
from prototypes.ep_v1 import initialize_F_matrix, initialize_U_edges


def test_initialize_F_matrix_N3():
    N = 3
    h = 1./N
    F = np.zeros((N+1,N+1)) 

    matrix_answer = np.array([[  0.  ,  51.28,  51.28,   0.  ],
                    [  0.  ,  25.64,  25.64,   0.  ],
                    [ -0.  , -25.64, -25.64,  -0.  ],
                    [ -0.  , -51.28, -51.28,  -0.  ]])

    assert np.array_equal(initialize_F_matrix(F, N), matrix_answer)
