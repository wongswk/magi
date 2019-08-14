import numpy as np


def vector(arma_vector):
    return np.ndarray(buffer=arma_vector, shape=(arma_vector.size(),))


def matrix(arma_matrix):
    return np.ndarray(buffer=arma_matrix, shape=(arma_matrix.n_rows, arma_matrix.n_cols), order="F")


def cube(arma_cube):
    return np.ndarray(buffer=arma_cube, shape=(arma_cube.n_rows, arma_cube.n_cols, arma_cube.n_slices), order="F")


def integer_vector(arma_integer_vector):
    return np.ndarray(buffer=arma_integer_vector, shape=(arma_integer_vector.size(),), dtype=int)


def integer_matrix(arma_integer_matrix):
    return np.ndarray(
        buffer=arma_integer_matrix, shape=(arma_integer_matrix.n_rows, arma_integer_matrix.n_cols), order="F", dtype=int
    )
