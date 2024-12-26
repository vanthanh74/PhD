#!/usr/local/bin/python
import numpy as np

def compute_error_f_K(seq_ref, seq_f_K):
    e = np.asarray(seq_ref) - np.asarray(seq_f_K)
    return [np.linalg.norm(v, ord=2) for v in e]

def compute_parareal_error(seq_ref, X, pitch):
    K = len(X)
    err = []

    for k in range(K):
        e = np.asarray(seq_ref[::pitch]) - np.asarray(X[k])
        errK = [np.linalg.norm(v, ord=2) for v in e ]
        err.append(errK)

    return err
