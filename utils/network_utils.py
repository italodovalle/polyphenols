#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 1st, 2020

@author: italodovalle
"""

import sys
import os

import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import random
from random import shuffle
from scipy import stats
from scipy import sparse

from progressbar import ProgressBar
import pickle

class NetworkUtilsInputError(Exception):
    pass

def pick_random_nodes_matching_selected(network, bins, nodes_selected, n_random,
                                        degree_aware=True, connected=False,
                                        seed=None):
    """
    Get random nodes with the same size as original set
    Adapted from https://github.com/emreg00/toolbox

    Args:
        network (nx.Graph):
        bins (list): bins containing nodes with similar degree
        nodes_selected:
        n_random (int): number of random selections
        degree_aware (bool): degree-preserving random selection
        connected (bool):
        seed (int): random seed
    Returns:
        list of lists
    """
    if seed is not None:
        random.seed(seed)
    values = []
    nodes = network.nodes()
    for i in range(n_random):
        if degree_aware:
            if connected:
               raise ValueError("Not implemented!")
            nodes_random = set()
            node_to_equivalent_nodes = get_degree_equivalents(nodes_selected, bins, network)
            for node, equivalent_nodes in node_to_equivalent_nodes.items():
               #nodes_random.append(random.choice(equivalent_nodes))
               chosen = random.choice(equivalent_nodes)
               for k in range(20): # Try to find a distinct node (at most 20 times)
                   if chosen in nodes_random:
                     chosen = random.choice(equivalent_nodes)
               nodes_random.add(chosen)
            nodes_random = list(nodes_random)
        else:
            if connected:
               nodes_random = [ random.choice(nodes) ]
               k = 1
               while True:
                   if k == len(nodes_selected):
                     break
                   node_random = random.choice(nodes_random)
                   node_selected = random.choice(network.neighbors(node_random))
                   if node_selected in nodes_random:
                     continue
                   nodes_random.append(node_selected)
                   k += 1
            else:
               nodes_random = random.sample(nodes, len(nodes_selected))
        values.append(nodes_random)
    return values

def get_degree_equivalents(seeds, bins, g):

    """
    Get other nodes that match degree as seed
    Adapted from https://github.com/emreg00/toolbox
    Args:
        seeds (list): nodes
        bins (list): bins containing nodes with similar degree
        g (nx.Graph): graph object
    """

    seed_to_nodes = {}
    for seed in seeds:
        d = g.degree(seed)
        for l, h, nodes in bins:
            if l <= d and h >= d:
               mod_nodes = list(nodes)
               mod_nodes.remove(seed)
               seed_to_nodes[seed] = mod_nodes
               break
    return seed_to_nodes

def get_random_nodes(nodes, network, bins=None, n_random=1000, min_bin_size=100,
                    degree_aware=True, seed=None):
    """
    Get a random selection of nodes 
    Adapted from https://github.com/emreg00/toolbox
    Args:
        nodes (list): list of nodes
        network (nx.Graph): graph
        bins (list): bins containing nodes with similar degree
        n_random (int): number of random samples
        min_bin_size (int): minimum size of bin of nodes with similar degree 
        degree_aware (bool): whether selection is degree preserving or not
        seed (int): random seed
    """

    if bins is None:
        # Get degree bins of the network
        bins = get_degree_binning(network, min_bin_size)
    nodes_random = pick_random_nodes_matching_selected(network, bins, nodes, n_random, degree_aware, seed=seed)
    return nodes_random


def calculate_proximity(network, nodes_from, nodes_to,
                        nodes_from_random=None, nodes_to_random=None,
                        bins=None, n_random=1000, min_bin_size=100,
                        seed=452456, sp = None, node2index = None,lengths = None):

    """
    Calculate network proximity measures (d_closest and d_shortest)
    between two sets of nodes
    Args:
        network (nx.Graph)
        nodes_from (iterable): source nodes
        nodes_to (iterable): target nodes
        nodes_from_random (iterable, optional)
        nodes_to_random (iterable, optional)
        bins (iterable, optional)
        n_random (int, default = 1000): number of random iteractions
        min_bin_size(int, default = 100): minimum size for bins of nodes with similar degree
        seed (int, default = 452456): random seed
        sp (np.ndarray, optional): pre-computed pairwirse shortest path matrix
        node2id (dict, optional): node 2 index mapping in shortest path matrix
        lengths (dict, optional): 
    Returns:
        dict: with proximity measures
    """


    nodes_network = set(network.nodes())
    if len(set(nodes_from) & nodes_network) == 0 or len(set(nodes_to) & nodes_network) == 0:
        return None # At least one of the node group not in network

    d = calculate_distances(network, nodes_from, nodes_to,sp, node2index)

    if n_random:

        if bins is None and (nodes_from_random is None or nodes_to_random is None):
            bins = get_degree_binning(network, min_bin_size, lengths) # if lengths is given, it will only use those nodes
        if nodes_from_random is None:
            nodes_from_random = get_random_nodes(nodes_from, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        if nodes_to_random is None:
            nodes_to_random = get_random_nodes(nodes_to, network, bins = bins, n_random = n_random, min_bin_size = min_bin_size, seed = seed)
        random_values_list = list(zip(nodes_from_random, nodes_to_random))
        #values = np.empty(len(nodes_from_random)) #n_random
        null = []
        for i, values_random in enumerate(random_values_list):
            nodes_from, nodes_to = values_random
            res = calculate_distances(network, nodes_from, nodes_to,sp, node2index)
            null.append(res)


        null_s = []
        null_c = []
        for i in range(len(null)):
            null_s.append(null[i]['shortest'])
            null_c.append(null[i]['closest'])
        
        
        with np.errstate(divide='ignore', invalid='ignore'):

            d['avg_shortest'],d['std_shortest'] = np.mean(null_s), np.std(null_s)
            d['z_shortest'] = (d['shortest'] - d['avg_shortest'])/d['std_shortest']

            d['avg_closest'],d['std_closest'] = np.mean(null_c), np.std(null_c)
            d['z_closest'] = (d['closest'] - d['avg_closest'])/d['std_closest']

    return (d)


def get_degree_binning(g, bin_size, lengths=None):

    """
    Organize nodes in bins according to their degree
    Adapted from https://github.com/emreg00/toolbox

    Args:
        g (nx.Graph):
        bin_size (int):
        lengths (dict)
    """
    degree_to_nodes = {}
    degrees = dict(g.degree())
    for node, degree in degrees.items():
        if lengths is not None and node not in lengths:
            continue
        degree_to_nodes.setdefault(degree, []).append(node)
    values = list(degree_to_nodes.keys())
    values.sort()
    bins = []
    i = 0
    while i < len(values):
        low = values[i]
        val = degree_to_nodes[values[i]]
        while len(val) < bin_size:
            i += 1
            if i == len(values):
               break
            val.extend(degree_to_nodes[values[i]])
        if i == len(values):
            i -= 1
        high = values[i]
        i += 1
        #print low, high, len(val)
        if len(val) < bin_size:
            low_, high_, val_ = bins[-1]
            bins[-1] = (low_, high, val_ + val)
        else:
            bins.append((low, high, val))
    return bins


def parse_interactome(infile, sep='\t', header=False, columns=[], lcc = False,
                      dataframe=False):

    """
    Parse edgelist from file or pandas.dataframe into a networkx.graph object
    Args:
        infile (str or pandas.dataframe): path to table or pandas.dataframe
        sep (str): column separator if input is file
        header (bool): indicates if table has header
        columns (list): names of columns for source and target nodes
        lcc (bool): indicates whether to return only lcc or not
        dataframe (bool): if True input must be a pandas.dataframe
    Returns:
        networkx.graph object
    """

    if dataframe:
        dt = infile
    else:
        if header:
            dt = pd.read_table(infile,sep = sep)
        else:
            dt = pd.read_table(infile,sep = sep,header=None)

    if header:
        edges = zip(dt[columns[0]], dt[columns[1]])
    else:
        edges = zip(dt[0], dt[1])

    G = nx.Graph()
    G.add_edges_from(edges)

    if lcc:
        g = list(connected_component_subgraphs(G))[0]
        #print (len(g.nodes()), 'nodes')
        #print (len(g.edges()), 'edges')
        return(g)
    else:
        #print (len(G.nodes()), 'nodes')
        #print (len(G.edges()), 'edges')
        return(G)


def connected_component_subgraphs(G, copy=True):
    ## this function was removed from latest versions of networkx!!
    for c in nx.connected_components(G):
        if copy:
            yield G.subgraph(c).copy()
        else:
            yield G.subgraph(c)

def get_lcc(G,S):
    """
    Get largest connect component induced by nodes in `S` in the graph `G`
    Args:
        S (list): set of source nodes
        G (nx.Graph): interactome
    """
    if len(S) == 0:
        return (nx.Graph())
    else:
        g = nx.subgraph(G,S)
        if len(g.nodes()) > 0:
            lcc = max(connected_component_subgraphs(g), key=len)
            return (lcc)
        else:
            return(g)



def get_lcc_significance(G,seeds,n_random=1000, min_bin_size = 100,
                            lengths = None,seed=452456):

    """
    Get significance of the size of the largest connected component
    induced by nodes in `seeds` in the graph `G`. At `n_random` iterations,
    the function gets a degree preserving random selection of nodes with
    the same size as `seeds`. After all iterations, a z-score is calculated
    using the lcc sizes of the random distribution that was obtained.
    
    Args:
        G (nx.Graph): interactome
        seeds (list): nodes
        n_random (int): number of random samples
        min_bin_size (int): minimum size of bins containing nodes with similar degree
        lenghts
        seed (int): random seed
    """

    # getting all genes in the network
    all_genes = G.nodes()


    seeds = list(set(seeds) & set(G.nodes()))
    #number_of_seed_genes = len(set(seeds) & set(all_genes))

    l_list  = []

    bins = get_degree_binning(G, min_bin_size, lengths)
    nodes_random = get_random_nodes(seeds, G,
                                    bins = bins, n_random = n_random,
                                    min_bin_size = min_bin_size,
                                    seed = seed)

    # simulations with randomly distributed seed nodes
    for i, rand_seeds in enumerate(nodes_random):

        # get rand lcc
        lcc = get_lcc(G,rand_seeds)
        lcc = len(lcc.nodes())
        l_list.append(lcc)


    # get the actual value
    lcc_observed = get_lcc(G,seeds)
    lcc_observed_size = len(lcc_observed.nodes())

    # get the lcc z-score:
    l_mean = np.mean(l_list)
    l_std  = np.std(l_list)

    if l_std == 0:
        z_score = float('nan')
    else:
        z_score = (1.*lcc_observed_size - l_mean)/l_std


    return ({'lcc_size':lcc_observed_size, 'z_score':z_score,
             'avg_size':l_mean, 'std_size':l_std})



def read_edgelist (infile,sep=' ',header=False):
    if header:
        lines = open(infile, 'r').readlines()[1:]
    else:
        lines = open(infile, 'r').readlines()
    edges = []
    for line in lines:
        a, b = line.rstrip().split(sep)
        edges.append((a,b))
    g = nx.from_edgelist(edges)
    return(g)

def calculate_distances (G, nodes_from, nodes_to,
                         sp=None, node2index=None):

    """
    Calculates the closest and shortest distance from
    `nodes_from` to `nodes_to` in the graph `G`.

    Args:

        G (nx.Graph):
        nodes_from (list):
        nodes_to (list):
        sp (np.matrix, optional): matrix with shortest paths pre-calculated
        index2node (dict, optional): index to node id mapping
    """

    ds = defaultdict(dict)


    for i in nodes_from:
        for j in nodes_to:
            if i == j:
                ds[i][j] = 0
            else:
                if nx.has_path(G,i, j):
                    if sp is None:
                        ds[i][j] = nx.shortest_path_length(G,i, j)
                    else:
                        ds[i][j] = sp[node2index[i],node2index[j]]
                else:
                    ds[i][j] = float('nan')

    ds = pd.DataFrame.from_dict(ds)
    # nodes_to: rows
    # nodes_from: columns

    dic = {}

    dic['shortest'] = ds.mean().mean()
    dic['closest'] = ds.min().mean()

    return (dic)
