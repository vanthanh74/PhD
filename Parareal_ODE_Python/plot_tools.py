#!/usr/bin/python
# -*- coding: utf8 -*-
from __future__ import (absolute_import, print_function)
import argparse
import sys
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def plot_solution(ref
                , X
                , title = "Solution"
                , xlabel = "x"
                , ylabel = "y"
                , xlim = None
                , ylim = None
                , xscale = 'linear'
                , yscale = 'linear'
                , colorPalette=['#0054AF','#33cc33','#cc3300','#cc9900','#be41ad','#000000','#b3c6ff','#e6e600','#ff9933','#cc9900']
                , legendlocation = 'best'
                , outputfile = 'solution.pdf'):

    K = len(X)
    pos_x_ref = [v[0] for v in ref]
    pos_y_ref = [v[1] for v in ref]

    # plt.scatter(pos_x_ref, pos_y_ref)
    # plt.show()

    pos_x = []
    pos_y = []
    for k in range(K):
        pos_x.append([v[0] for v in X[k]])
        pos_y.append([v[1] for v in X[k]])
        # plt.scatter(pos_x[k], pos_y[k])

    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(0,0))
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    for k in range(K):
        ax.plot(pos_x[k],pos_y[k],label="parareal k="+str(k), linewidth=2.0,
                    marker='o', markersize=4.0,color=colorPalette[k])
    ax.plot(pos_x_ref, pos_y_ref, label="sequential", linewidth=2.0,
                    marker='o', markersize=4.0,color=colorPalette[K])

    if legendlocation != None:
        lines, labels = ax.get_legend_handles_labels()
        plt.legend(lines , labels,
               loc=legendlocation, shadow=True,fancybox=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)
    plt.clf()

def plot_conv(t
                , err_parareal
                , title = r'Error parareal $||(x,y)(T_N)-X^N_k||$'
                , xlabel = "Time t"
                , ylabel = ''
                , xlim = None
                , ylim = [1e-14,1.]
                , xscale = 'linear'
                , yscale = 'log'
                , colorPalette=['#0054AF','#33cc33','#cc3300','#cc9900','#be41ad','#000000','#b3c6ff','#e6e600','#ff9933','#cc9900']
                , legendlocation = 'best'
                , fontsize = 16
                , outputfile = 'conv.pdf'):
    K = len(err_parareal)

    fig, ax = plt.subplots()
    if title != None:
        plt.title(title, fontsize=fontsize)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.ticklabel_format(style='sci', scilimits=(0, 0))
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    plt.xticks(fontsize=fontsize)  # Font size abscissas
    plt.yticks(fontsize=fontsize)

    for k in range(K):
        ax.plot(t, err_parareal[k], label=r'$k=$' + str(k), linewidth=2.0,
                marker='o', markersize=4.0, color=colorPalette[k])

    if legendlocation != None:
        lines, labels = ax.get_legend_handles_labels()
        plt.legend(lines, labels,
                   loc=legendlocation,
                   shadow=True,
                   fancybox=True,
                   fontsize=14)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)
    plt.clf()

def plot_conv_uniform(err_parareal
                , title = "Error parareal sequence"
                , xlabel = "k"
                , ylabel = r'max err in k'
                , xlim = None
                , ylim = [1e-18,1.]
                , xscale = 'linear'
                , yscale = 'log'
                , colorPalette=['#0054AF','#33cc33','#cc3300','#cc9900','#be41ad','#000000','#b3c6ff','#e6e600','#ff9933','#cc9900']
                , legendlocation = 'best'
                , outputfile = 'conv.pdf'):
    K = len(err_parareal)

    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(0, 0))
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    maxVal = []
    for k in range(K):
        maxVal.append(np.max(err_parareal[k]))
    print(maxVal)

    ax.plot(maxVal, label="max err", linewidth=2.0,
            marker='o', markersize=4.0, color=colorPalette[k])

    if legendlocation != None:
        lines, labels = ax.get_legend_handles_labels()
        plt.legend(lines, labels,
                   loc=legendlocation, shadow=True, fancybox=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)
    plt.clf()

def plot_err_ref_vs_f_K(t
                , err_f_K
                , title = r'Error $F_K$'
                , xlabel = "t"
                , ylabel = r'err $l_2$'
                , xlim = None
                , ylim = None
                , xscale = 'linear'
                , yscale = 'log'
                , colorPalette=['#0054AF','#33cc33','#cc3300','#cc9900','#be41ad','#000000','#b3c6ff','#e6e600','#ff9933','#cc9900']
                , legendlocation = 'best'
                , outputfile = 'err-f-K.pdf'):

    fig, ax = plt.subplots()
    if title != None:
        plt.title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(style='sci', scilimits=(0, 0))
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    ax.plot(t, err_f_K, label="err", linewidth=2.0,
                marker='o', markersize=4.0, color=colorPalette[0])

    if legendlocation != None:
        lines, labels = ax.get_legend_handles_labels()
        plt.legend(lines, labels,
                   loc=legendlocation, shadow=True, fancybox=True)
    if xlim != None:
        plt.xlim(xlim)
    if ylim != None:
        plt.ylim(ylim)
    plt.savefig(outputfile)
    plt.clf()

def generate_plots(t, ref, X, err_parareal, err_f_K, refine):
    plot_solution(ref, X)
    title = ""
    outputfile = ""
    outputfile_u = ""
    if refine:
        title = "Error parareal $||(x,y)(T_N)-X^N_k||$ (adaptive $\delta t$)"
        outputfile = "conv-with-refinement.pdf"
        outputfile_u = "conv-with-refinement-uniform.pdf"
    else:
        title = "Error parareal $||(x,y)(T_N)-X^N_k||$ (fixed $\delta t$)"
        outputfile = "conv-no-refinement.pdf"
        outputfile_u = "conv-no-refinement-uniform.pdf"
    plot_err_ref_vs_f_K(t, err_f_K)
    plot_conv(t, err_parareal, title=title, outputfile=outputfile)
    plot_conv_uniform(err_parareal, title=title, outputfile=outputfile_u)
