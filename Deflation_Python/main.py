# Generate several results of the distributed svd for a combination of parameters
# Run without argument to see the full list

from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
import os
import anytree as tree
import csv
import argparse
import scipy.linalg

import utils
import qr

def min_max(opt):

    data_dir = "data"

    fp = open(data_dir + "/min_max.out", 'w')
    fp.write("method\t"
             "inc\t"
             "nb_core\t"
             "matrix_name\t"
             "min_ratio\t"
             "max_ratio\t"
             "size\t"
             "error_fro\t"
             "requested_k\t"
             "increase_k\t"
             "increase_pct\t"
             "sub_k\t"
             "requested_error\t"
             "error_gathering\t"
             "error_to_qrcp\n")

    requested_k = opt.requested_ranks_single

    inc = 0
    for matrix_path in opt.matrix_list:
        matrix_name = os.path.split(os.path.splitext(matrix_path)[0])[1]

        A = utils.read_matrix(matrix_path)
        print("\nMatrix", matrix_path, A.shape)

        _, s, _ = np.linalg.svd(A)
        requested_error = np.linalg.norm(s[opt.requested_ranks_single + 1:])

        for nb_cores in opt.nbs_cores:
            opt.nbs_cores_single = nb_cores
            print("Run {} Parameters : requested_k {} nb_cores {}".format(inc, requested_k, nb_cores))

            selected_columns, error_tree = qr.split_qr_rank(A, requested_k, opt.nbs_cores_single, opt)
            #utils.render_tree(error_tree)
            Q1,_ = np.linalg.qr(A[:,np.setdiff1d(range(A.shape[1]), selected_columns)])
            s_approx = np.linalg.svd(A - Q1 * Q1.T *A)[1][:requested_k]

            u_approx_qr, v_approx_qr, p_qr = scipy.linalg.qr(A, pivoting=True)
            s_approx_qr = np.linalg.svd(v_approx_qr[-requested_k:,-requested_k:])[1]

            max_ratio = np.max(np.abs(s[-requested_k:] / s_approx))
            min_ratio = np.min(np.abs(s[-requested_k:] / s_approx))

            error_fro = 0

            submatrices_ranks = 0
            error_gathering = requested_error - np.sum([n.error for n in tree.PreOrderIter(error_tree) if len(n.children) == 0])

            fp.write(("{}\t" * 7 + "{}\n").format(
                inc,
                nb_cores,
                matrix_name,
                min_ratio,
                max_ratio,
                error_fro,
                requested_k + 1,
                error_gathering))

            inc += 1

    fp.close()


def single_rank(opt):
    
    data_dir = "data"

    rc('text', usetex=True)
    font = {'weight': 'bold',
            'family': 'DejaVu Sans',
            'size': 15}

    rc('font', **font)

    requested_rank = opt.requested_ranks_single
    matrix_name = utils.get_matrix_name(opt.matrix_list_single)

    A = utils.read_matrix(opt.matrix_list_single)
    print("Input matrix has size", A.shape)

    # Compute approximation
    _, s, _ = np.linalg.svd(A)
    requested_error = np.linalg.norm(s[requested_rank:])

    selected_columns, error_tree = qr.split_qr_rank(A, requested_rank, opt.nbs_cores_single, opt)
    Q1,_ = np.linalg.qr(A[:,np.setdiff1d(range(A.shape[1]), selected_columns)])
    s_approx = np.linalg.svd(A - Q1 @ Q1.T @ A)[1][:requested_rank]

    utils.render_tree(error_tree)

    sv_indices = np.array(range(len(s) - requested_rank, len(s)))

    # Compute full QRCP for comparison
    q_approx_qr, r_approx_qr, p_qr = scipy.linalg.qr(A, pivoting=True)
    s_approx_qr = np.linalg.svd(r_approx_qr[-requested_rank:,-requested_rank:])[1]
    print("qrcp columns :", p_qr[-requested_rank:])

    # Plot
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.plot(sv_indices + 1, s[sv_indices], "x-", label="SVD")
    ax.plot(sv_indices + 1, s_approx_qr, "+-", label="QRCP")
    ax.plot(sv_indices + 1, np.abs(np.diagonal(r_approx_qr)[sv_indices]), "o-", alpha=0.3, label="QRCP diag")
    ax.plot(sv_indices + 1, s_approx, "^-", label="QRTP", fillstyle="none")
    #ax.set_xlim(left=sv_indices[0] + 1, right=sv_indices[-1] + 1)
    ax.set_xlabel(r"$i$")
    ax.legend(loc="upper center")

    ax2 = ax.twinx()
    ax2.plot(sv_indices + 1, np.abs(s[sv_indices] / s_approx), "ro", label=r"$\sigma_i(A)/\sigma_i(A_k)$", alpha=0.5)
    #ax2.set_yscale("log")
    ax2.grid(True)

    ax2.set_ylabel("Error")
    ax2.legend(loc="upper right")
    ax.set_yscale("log")

    plt.savefig("pdf/" + matrix_name + "_" + str(requested_rank) + ".pdf")
    plt.show()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("matrix_list", metavar="matrix", nargs="+", help="matrix path")
    parser.add_argument("--requested-rank", "-k", nargs="+", dest="requested_ranks", default=[10], type=int)
    parser.add_argument("--cores", "-c", dest="nbs_cores", default=["4"], nargs="+", type=str, help="n for square splitting, nxm for rectangular splitting, append 'c' at then end for column then row gathering instead of nested squares gathering")
    parser.add_argument("--single-rank", dest="single_rank", action='store_true', help="single_rank")
    parser.add_argument("--disp-cols", dest="selected_cols", action='store_true', help="Display selected columns at each step")

    args = parser.parse_args()
    args.nbs_cores = [i.split(',') for i in args.nbs_cores]

    if len(args.requested_ranks) != 1:
        raise ValueError("Expected only one requested rank")
    args.requested_ranks_single = args.requested_ranks[0]

    if args.single_rank:
        if len(args.matrix_list) != 1:
            raise ValueError("Expected only one matrix for single-rank")
        args.matrix_list_single = args.matrix_list[0]

        if len(args.nbs_cores) != 1:
            raise ValueError("Expected only one number or cores (nb-cores) for single-rank")
        args.nbs_cores_single = args.nbs_cores[0]

        single_rank(args)
    else:
        if args.figure:
            ValueError("fig option not available for min_max")
        min_max(args)
