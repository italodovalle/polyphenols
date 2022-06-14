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

import utils.network_utils as network_utils
import sys
import logging


def get_zscores (disease_chemical, n_random, outdir,
                     G, chemical2genes, disease2genes):

        disease, chemical = disease_chemical
        if type(chemical) == str:
            chemical = chemical.lower()
        else:
            chemical = str(int(chemical))
        nodes_from = set(chemical2genes[chemical]) & set(G.nodes())
        nodes_to = set(disease2genes[disease]) & set(G.nodes())
        min_bin_size = 2 * max(len(nodes_from), len(nodes_to))

        dic = network_utils.calculate_proximity(G, nodes_from,nodes_to,
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
            table['avg_shortest'] = dic['avg_shortest']
            table['std_shortest'] = dic['std_shortest']
            table['avg_closest'] = dic['avg_closest']
            table['std_closest'] = dic['std_closest']
        else:
            table['shortest'] = float('nan')
            table['closest'] = float('nan')
            table['z_shortest'] = float('nan')
            table['z_closest'] = float('nan')
            table['avg_shortest'] = float('nan')
            table['std_shortest'] = float('nan')
            table['avg_closest'] = float('nan')
            table['std_closest'] = float('nan')

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
            print ('%d pre-calculated chemical disease pairs loaded'%len(finished))
            


    samples = []
    for disease in disease2genes.keys():
        for molecule in chemical2genes.keys():
            pair = (disease,molecule)
            if not pair in finished:
                samples.append((disease,molecule))
    print ('%d reamining chemical-disease pairs'%len(samples))
    
    if len(samples) == 0:
        print ('No chemical-disease pairs for calculation!')
        sys.exit()


    if test_run:
        if len(samples) > ncpus:
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

    description = 'Calculate Network Proximity'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i', required=True, dest='interactome',
                        action='store',help='interactome file')
    parser.add_argument('-d', required=True, dest='diseasegenes',action='store',
                        help='disease genes file')
    parser.add_argument('-t', required=True, dest='targets',action='store',
                        help='chemical protein interactions')
    parser.add_argument('-c', dest='cores', action='store',
                        help='number of cores to use, default = 6',
                        default= 6)
    parser.add_argument('-r', dest='nrandom', action='store',
                        help='n random permutations, default = 1000',
                        default= 1000)
    parser.add_argument('-o', dest='outfolder', action='store',
                        help='outfolder', default='out')
    parser.add_argument('-sn', dest='source_node', action='store',
                        help='name source node Interactome, default = proteinA',
                        default= 'proteinA')
    parser.add_argument('-tn', dest='target_node', action='store',
                        help='name target node Interactome, default = proteinB',
                        default= 'proteinB')
    parser.add_argument('-lcc', dest = 'lcc', action = 'store',
                        help = 'Consider only Interactome LCC, default = True',
                        default = True)
    parser.add_argument('-sep', dest='separator', action='store',
                        help = 'separator Interactome, default = ","',
                        default = ',')
    parser.add_argument('-test', dest = 'test_run', action='store',
                        help = 'test run')
    parser.add_argument('-chem', dest = 'chemical_col', action = 'store',
                        default = 'chemical', help = 'chemical column in chemical-protein interaction file')
    parser.add_argument('-prot', dest = 'protein_col', action = 'store',
                        default = 'entrezid', help = 'chemical column in chemical-protein interaction file')
    parser.add_argument('-tmp', dest = 'tmp_dir', action = 'store', 
                        help = 'temporary folder', default='tmp')




    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()


    infile = args.interactome
    disease_genes_file = args.diseasegenes
    polyphenl_targets_file = args.targets
    ncpus = int(args.cores)
    n_random = int(args.nrandom)
    outdir_tmp_files = os.path.abspath(args.tmp_dir) + '/'
    sep = args.separator
    lcc = bool(args.lcc)
    columns = [args.source_node, args.target_node]
    test_run = bool(args.test_run)
    chemical_col = args.chemical_col
    protein_col = args.protein_col
    outdir = os.path.abspath(args.outfolder) + '/'
    header=True

    final_outfile = outdir + 'merged_zscore_proximity.csv'


    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(outdir_tmp_files):
        os.mkdir(outdir_tmp_files)
    
    

    logging.basicConfig(format='%(asctime)s - %(message)s',
                    datefmt='%d-%b-%y %H:%M:%S',level=logging.INFO)


    G = network_utils.parse_interactome(infile, sep, header, columns, lcc=lcc)



    logging.info('loading network: %d nodes, %d edges'%(len(G.nodes()),
                                                        len(G.edges())))

    ## create dictionary
    ## disease names are keys and elements are disease genes

    disease2genes = defaultdict(list)
    
    dg = pd.read_csv(disease_genes_file)
    
    for i in dg.index:
        disease2genes[dg.disease.loc[i]].append(dg.entrez_id.loc[i])


    ## create dictionary
    ## polyphenol names are keys and elements are polyphenol targets


    polyphenol = pd.read_csv(polyphenl_targets_file)
    polyphenol[chemical_col] = polyphenol[chemical_col].astype(str)

    chemical2genes = defaultdict(list)
    for i in polyphenol.index:
        name = polyphenol[chemical_col].loc[i]
        chemical2genes[name].append(polyphenol[protein_col].loc[i])

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
