# This file contains methods to compute the distributed QR

import numpy as np
import anytree as tree
import scipy.linalg
import matplotlib.pyplot as plt
import scipy
import anytree as tree
import matplotlib.pyplot as plt

import utils


def select_columns_lapack(A, k):
    q, r, p = scipy.linalg.qr(A, pivoting=True)
    residue = np.linalg.norm(q[:, k:].T @ A)
    selected_cols = p[-k:]
    return selected_cols, residue


def truncated_qr_lapack(A, k):
    q, r, p = scipy.linalg.qr(A, pivoting=True)
    permutation = np.eye(A.shape[0], A.shape[1])[:, p]
    r = r @ permutation.T
    #print(np.arange(A.shape[0]) @ permutation[:, :k])
    residue = np.linalg.norm(r[k:, :])

    error_node = tree.Node("a", requested=0, error=residue, factor=1, increase_k=0, rank=k)

    return q[:, -k:], r[-k:, :], error_node, p[-k:]

def split_qr_rank(A, k, nb_cores, opt, depth=0):
    # Divide matrix
    split_2d = np.fromstring(nb_cores[0], dtype=int, sep="x")

    nb_cores_cur = 0
    if len(split_2d) == 2:
        posx, posy, start_i, start_j = utils.split_indexes_2D(A.shape[0], A.shape[1], split_2d[0], split_2d[1])
        nb_cores_cur = split_2d[0] * split_2d[1]
    elif len(split_2d) == 1:
        posx, posy, start_i, start_j = utils.split_indexes(A.shape[0], A.shape[1], split_2d[0])
        nb_cores_cur = split_2d[0]
    else:
        print("split_2d", split_2d)
        raise ValueError("Specified splitting {} is incorrect. Format is 4, 4x4 or 4x4c".format(nb_cores[0]))

    a = []
    selected_columns = set()
    selected_cols_local_offset = []

    # huge error to make use of increase_k directly
    error = np.linalg.norm(A) + 1

    current_node = tree.Node("a", requested=0, increase_k=0, factor=1)

    # Split A
    for i in range(nb_cores_cur):
        a.append(A[start_i[posx[i]]:start_i[posx[i] + 1], start_j[posy[i]]:start_j[posy[i] + 1]])

    for i in range(nb_cores_cur):

        # Compute truncated pivoted QR to get submatrice columns
        if len(nb_cores) == 1:
            selected_cols_local, residue = select_columns_lapack(a[i], k)
            tree.Node("a", requested_k=k, rank=k, parent=current_node, requested=error, error=residue, posy=posy[i])
        else:
            selected_cols_local, residue, child_node = split_qr_rank(a[i], k, nb_cores[1:], opt, depth + 1)
            child_node.posy = posy[i]
            child_node.parent = current_node

        if len(selected_cols_local) == 0:
            selected_cols_local = [0]

        # Add offset to get absolute column number
        if opt.selected_cols:
            print("depth", depth, "submatrix", i, [x + start_j[posy[i]] for x in selected_cols_local])
        selected_cols_local_offset.append({x + start_j[posy[i]] for x in selected_cols_local if (x + start_j[posy[i]] not in selected_columns)})

        # Store to global selected column list
        selected_columns = selected_columns | selected_cols_local_offset[i]

    a_mini = A[:, list(selected_columns)]

    if depth > 0:
        step2_selected_cols, residue = select_columns_lapack(a_mini, k)

        current_node.rank = k
        current_node.error = residue

        return np.array(list(selected_columns))[step2_selected_cols], residue, current_node

    Q, R, residue, step2_selected_cols = truncated_qr_lapack(a_mini, k)
    current_node.error = residue.error
    current_node.rank = residue.rank

    print("final columns", np.array(list(selected_columns))[step2_selected_cols])

    return np.array(list(selected_columns))[step2_selected_cols], current_node

def randomize(A):
    k = 10
    n = A.shape[0]
    n2 = 2 @ np.floor(np.log2(n))
    delta = 0.9
    epsilon = 0.1
    s = np.ceil((10 / epsilon) @ (np.sqrt(k) + np.sqrt(8 @ np.log(n2 / delta))) @ 2 @ np.log(k / delta))
    print("computed s", s)

    P = np.eye(n)[:, :s]
    H2 = np.scipy.hadamard(n2)
    H = np.concatenate((H2, np.concatenate((np.zeroes(n - n2, n2), np.eye(n - n2)), axis=1)))
    D = np.diag((np.rands(n, 1) > 0.5) @ 2 - 1)
    return P @ H @ D @ A

