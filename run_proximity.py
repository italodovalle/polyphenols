#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import networkx as nx
import numpy as np
import pandas as pd
from collections import defaultdict
from multiprocessing import Pool
import random
import os
import time
from functools import partial
import re
import argparse
import sys

from toolbox.guney_code import wrappers
from toolbox.guney_code import network_utilities
import toolbox.network_utils as network_utils
import sys
import logging


def get_zscores (disease_chemical, n_random, outdir,
                     G, chemical2genes, disease2genes):

        disease, chemical = disease_chemical
        chemical = chemical.lower()
        nodes_from = set(chemical2genes[chemical]) & set(G.nodes())
        nodes_to = set(disease2genes[disease]) & set(G.nodes())
        min_bin_size = 2 * max(len(nodes_from), len(nodes_to))

        dic = network_utils.calculate_proximity_italo(G, nodes_from,nodes_to,
                                                      n_random = n_random,
                                                      min_bin_size=min_bin_size)
        table = {}
        table['disease'] = disease
        table['n_mapped_disease'] = len(set(nodes_to))
        table['n_mapped_chemical'] = len(set(nodes_from))
        table['chemical'] = chemical
        if dic:
            table['shortest'] = dic['shortest']
            table['closest'] = dic['closest']
            table['z_shortest'] = dic['z_shortest']
            table['z_closest'] = dic['z_closest']
        else:
            table['shortest'] = float('nan')
            table['closest'] = float('nan')
            table['z_shortest'] = float('nan')
            table['z_closest'] = float('nan')

        df = pd.DataFrame.from_dict(table, orient='index').T

        if outdir:
            out_disease = re.sub('[^A-Za-z0-9]+', '', disease)
            out_chemical = re.sub('[^A-Za-z0-9]+', '', chemical)
            outname = '%s_%s'%(out_disease, out_chemical)
            outdir = os.path.abspath(outdir)
            df.to_csv(outdir + '/%s.csv'%outname)

        return(df)

def run_proximity (G, disease2genes, chemical2genes, ncpus = 15,
                   n_random = 10, outdir=None, test_run = False,
                   sp = None, node2index = None):


    finished = []
    ##retrive pairs that had their calculations already done
    if outdir:
        fs = []
        for i in os.listdir(outdir):
            if i.endswith('.csv'):
                x = pd.read_csv(outdir + '/' + i, index_col = 0)
                fs.append(x)

        if len(fs) > 0:
            ds = pd.concat(fs)
            finished = [(i,str(j)) for i,j in zip(ds.disease, ds.chemical)]

    #print (len(finished))
    #print (finished[:10])
    #sys.exit()


    samples = []
    for disease in disease2genes.keys():
        for molecule in chemical2genes.keys():
            pair = (disease,molecule)
            if not pair in finished:
                samples.append((disease,molecule))
    print ('%d Chemical-disease pairs'%len(samples))


    if test_run:
        samples = samples[:ncpus]

    p = Pool(ncpus)
    res = p.map(partial(get_zscores,
                  n_random = n_random, outdir=outdir,
                  G = G, chemical2genes = chemical2genes,
                  disease2genes = disease2genes), samples)
    p.close()

    df = pd.concat(res)
    return(df)

if __name__ == '__main__':

    ### Interactome params

    description = 'Run Proximity'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', required=True, dest='interactome',
                        action='store',help='interactome file')
    parser.add_argument('-d', required=True, dest='diseasegenes',action='store',
                        help='disease genes file')
    parser.add_argument('-t', required=True, dest='targets',action='store',
                        help='target genes file')
    parser.add_argument('-c', dest='cores', action='store',
                        help='number of cores to use, default = 6',
                        default= 6)
    parser.add_argument('-r', dest='nrandom', action='store',
                        help='n random permutations, default = 1000',
                        default= 1000)
    parser.add_argument('-o', required=True, dest='outfolder', action='store',
                        help='outfolder')
    parser.add_argument('-sn', dest='source_node', action='store',
                        help='name source node, default = proteinA',
                        default= 'proteinA')
    parser.add_argument('-tn', dest='target_node', action='store',
                        help='name target node, default = proteinB',
                        default= 'proteinB')
    parser.add_argument('-lcc', dest = 'lcc', action = 'store',
                        help = 'LCC, default = True',
                        default = True)
    parser.add_argument('-sep', dest='separator', action='store',
                        help = 'separator, default = ","',
                        default = ',')
    parser.add_argument('-test', dest = 'test_run', action='store',
                        help = 'test run, default = False',
                        default = False)



    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()


    ##--------------------------------
    #Input file fastq already trimmed
    infile = args.interactome
    disease_genes_file = args.diseasegenes
    polyphenl_targets_file = args.targets
    ncpus = int(args.cores)
    n_random = int(args.nrandom)
    outdir_tmp_files = os.path.abspath(args.outfolder) + '/'
    sep = args.separator
    lcc = bool(args.lcc)
    columns = [args.source_node, args.target_node]
    test_run = bool(args.test_run)

    final_outfile = outdir_tmp_files + 'merged_zscore_proximity.csv'

    header=True




    logging.basicConfig(format='%(asctime)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',level=logging.INFO)


    G = network_utils.parse_interactome(infile, sep, header, columns, lcc=lcc)



    logging.info('loading network: %d nodes, %d edges'%(len(G.nodes()),
                                                        len(G.edges())))

    ## create dictionary
    ## disease names are keys and elements are disease genes

    disease2genes = {}
    for i in open(disease_genes_file).readlines():
        v = i.rstrip().split('\t')
        disease = v[1]
        genes = v[2:]
        if len(genes) > 19:
            disease2genes[disease] = [int(i) for i in genes]

    ## create dictionary
    ## polyphenol names are keys and elements are polyphenol targets


    polyphenol = pd.read_csv(polyphenl_targets_file,index_col = 0)
    polyphenol = polyphenol.reset_index()
    #polyphenol = polyphenol[(polyphenol.experimental > 0) | (polyphenol.database > 0)]
    polyphenol = polyphenol[(polyphenol.experimental > 0)]
    chemical2genes = defaultdict(list)
    for i in polyphenol.index:
        name = str(int(polyphenol.chemical.loc[i]))
        chemical2genes[name].append(polyphenol.entrezid.loc[i])

    for i in chemical2genes.keys():
        chemical2genes[i] = list(set(chemical2genes[i]))


    if test_run:
        logging.info('Test Run')

        s = time.time()
        df = run_proximity(G, disease2genes, chemical2genes, ncpus = ncpus,
                       n_random = n_random, outdir=outdir_tmp_files, test_run=True)
        e = time.time() - s


        estimated_time = (e * (len(disease2genes) * len(chemical2genes)))/ncpus/3600

        logging.info('estimated end time: %f hours'%estimated_time)


    logging.info('Running analysis')


    s = time.time()
    df = run_proximity(G, disease2genes, chemical2genes, ncpus = ncpus,
                       n_random = n_random, outdir=outdir_tmp_files,
                       test_run=False)
    e = time.time() - s
    e = e/s

    logging.info('Finished: %f hours'%e)

    df.to_csv(final_outfile)
