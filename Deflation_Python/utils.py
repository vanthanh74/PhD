# Module of general utility for matrix manipulations


from scipy import io
import os
import anytree as tree
import scipy.sparse

import numpy as np

OUTPUT_PATH = "data/"

# return matrix of form
# A1 | A2
# ---|---
# A3 | A4
def matrix_from_blocks(A1, A2, A3, A4):
    return np.concatenate((
        np.concatenate((A1, A3)),
        np.concatenate((A2, A4))
    ), axis=1)

def matrix_from_blocks_generic(blocks, start_i, start_j):
    inc = 0
    lines = []
    for bloc in blocks:
        if start_j[inc] == 0:
            lines.append([])
        lines[start_i[inc]].append(bloc)
        inc += 1
    return np.concatenate([np.concatenate(lines[i], axis=1) for i in range(len(lines))])


# Read matrix for file
# * Matrix market with extension ".mtx"
# * other extensions will be treated as CSV file
# @param name the path to the matrix file
def read_matrix(name):
    extension = name.split('.')[-1]
    #print("Extension : {}".format(extension))

    if extension == 'mtx':
        sparseA = io.mmread(name)
        return np.asarray(scipy.sparse.csc_matrix.todense(sparseA))
    else:
        return np.loadtxt(open(name, "rb"), delimiter=",")


# return (2-norm(A-B), Froebenius (norm(A-B))
def compute_error(A, B, name):
    _, s, _ = np.linalg.svd(A - B)

    #print(name + " error (2-norm) : {:.2e}".format(s[0]))
    #print(name + " error (froebenius) : {:.2e}".format(np.linalg.norm(s)))

    return s[0], np.linalg.norm(s)


# Complete A with zeroes to match given shape (null complement)
# @param right place A right to the zeroes
# @param bottom place A under the zeroes
def nc(A, shape, right=0, bottom=0):
    n = shape[0]
    p = shape[1]
    result = np.copy(A)
    if result.shape[0] < n:
        if bottom:
            result = np.concatenate((np.zeros((n - A.shape[0], result.shape[1])), result))
        else:
            result = np.concatenate((result, np.zeros((n - A.shape[0], result.shape[1]))))
    if result.shape[1] < p:
        if right:
            result = np.concatenate((np.zeros((result.shape[0], p - A.shape[1])), result), axis=1)
        else:
            result = np.concatenate((result, np.zeros((result.shape[0], p - A.shape[1]))), axis=1)
    return result


def nc_middle(A, top=0, left=0, right=0, bottom=0):
    left_upper_zeroes = nc(A, (A.shape[0] + top, A.shape[1] + left), right=1, bottom=1)
    return nc(left_upper_zeroes, (left_upper_zeroes.shape[0] + bottom, left_upper_zeroes.shape[1] + right))


def nc_indexed(A, shape, x=0, y=0):
    left_upper_zeroes = nc(A, (A.shape[0] * (x + 1), A.shape[1] * (y + 1)), right=1, bottom=1)
    a = nc(left_upper_zeroes, shape)
    return a


def block_diagonal(blocks):
    A = blocks[0]
    for block in blocks[1:]:
        width = A.shape[0] + block.shape[0]
        A = np.concatenate((nc(A, (width, 0)), nc(block, (width, 0), bottom=1)), axis=1)

    return A


# Cut a matrix of size (n,m) into k submatrices
# for x in [0,k-1], submatrix x has indexes (start_i[x], start_j[x]) in the grid
# this matrix has lines (posx[start_i[x]] to posx[start_i[x] + 1]
# on columns posy[start_j[x]] to posy[start_j[x] + 1]) of the original matrix
def split_indexes(n, m, k):
    width = int(np.sqrt(k))
    height = k // width
    return split_indexes_2D(n, m, width, height)


# Provide number of blocks in a line (M) and in a column (N)
def split_indexes_2D(n, m, height, width):
    posx = np.linspace(0, n, height + 1, dtype=int)
    posy = np.linspace(0, m, width + 1, dtype=int)
    start_i = [i for i in range(height) for j in range(width)]
    start_j = [j for i in range(height) for j in range(width)]
    return start_i, start_j, posx, posy


def get_matrix_name(file_path):

    return os.path.split(os.path.splitext(file_path)[0])[1]


def get_data_dir(file_path):
    return OUTPUT_PATH + get_matrix_name(file_path)


def make_data_dir(file_path):
    data_dir = get_data_dir(file_path)

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    return data_dir

def render_tree(root):
    for pre, _, node in tree.RenderTree(root):
        if not hasattr(node, 'increase_k'):
            node.increase_k = 0
        if not hasattr(node, 'factor'):
            node.factor = 1
        print("{}rank {} + {} / error : {:.2} / requested {:.2e} * {}".format(pre, node.rank, node.increase_k, node.error, node.requested, node.factor))

# From https://stackoverflow.com/a/4710090/6400167
def delete_numerical_line(file_name):
    with open(file_name, "r") as f:
        lines = f.readlines()
    with open(file_name, "w") as f:
        for line in lines:
            if is_number(line):
                f.write(line)

# from https://stackoverflow.com/q/354038/6400167
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
