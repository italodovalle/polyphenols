#!/usr/bin/env python3

import pandas as pd
import numpy as np

from collections import defaultdict
import tqdm
import re
from multiprocessing import Pool
import time

import scipy
import networkx as nx
import re

import os
import random
import statsmodels.api as sm

from statsmodels.stats.multitest import multipletests
import utils.network_utils as network_utils
import logging
import argparse
import sys


def parse_gct(infile):
    
    """
    Parse .gct file obtained from https://clue.io/
    Args:
        infile (str): file path
    """
    
    dt = pd.read_table(infile, sep = '\t', skiprows=2, low_memory=False)
    cols = list(dt.columns[9:])
    desc = dt.loc[:7, cols]
    cols.append('id')
    exp_matrix = dt[cols]
    exp_matrix = dt.iloc[8:]
    attributes = defaultdict(dict)
    desc = dt.loc[:7, cols]
    for i in desc.columns[:-1]:
        for j in desc.index:
            attributes[i][desc['id'].loc[j]] = desc[i].loc[j]
    attributes = pd.DataFrame.from_dict(attributes, orient='index').reset_index()
    attributes['distil_cc_q75'] = attributes['distil_cc_q75'].astype(float)
    return(exp_matrix, attributes)

def get_gsea(geneset, nongeneset, df,dataframe=False):
    dic = {}
    N_H = len(geneset)
    N = df.shape[0]
    n_r = np.sum(df['value'][df.gene.isin(geneset)].abs())
    for rank in range(N):
        dx = df.iloc[:rank]
        s = list(set(geneset) & (set(dx['gene'])))
        p_hit = np.sum(dx[dx.gene.isin(s)]['value'].abs()/n_r)
        p_miss = 1.0/(N - N_H) * dx[~dx.gene.isin(geneset)].shape[0]
        dic[rank] = p_hit - p_miss
    dic = pd.DataFrame.from_dict(dic, orient='index')
    if dataframe:
        return(dic)
    else:
        return(dic[0].abs().max())
    
    
def get_gsea_performance(geneset, nongeneset, df):
    df = df.reset_index()
    ranks_mapped = list(df[df.gene.isin(geneset)].index)
    ranks = []
    for i in ranks_mapped:
        ranks.append(i - 1)
        ranks.append(i)
        ranks.append(i + 1)
    dic = {}
    N_H = len(geneset)
    N = df.shape[0]
    n_r = np.sum(df['value'][df.gene.isin(geneset)].abs())
    for rank in ranks:
        dx = df.iloc[:rank]
        s = list(set(geneset) & (set(dx['gene'])))
        p_hit = np.sum(dx[dx.gene.isin(s)]['value'].abs()/n_r)
        p_miss = 1.0/(N - N_H) * dx[~dx.gene.isin(geneset)].shape[0]
        dic[rank] = p_hit - p_miss
    dic = pd.DataFrame.from_dict(dic, orient='index')
    return(dic[0].abs().max())


def get_null(input):
    np.random.seed()
    n, df = input
    geneset = df['gene'].sample(n)
    nongeneset = list(set(df['gene']) - set(geneset))
    #es = get_gsea(geneset, nongeneset, df)
    es = get_gsea_performance(geneset, nongeneset, df)
    return(es)



if __name__ == '__main__':

    ### Interactome params

    description = 'Calculate Enrichment Score'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-indir', required=True, dest = 'indir',help='folder with .gct files')
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
                        help='outfolder', default='out_gsea')
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
    parser.add_argument('-chem', dest = 'chemical_col', action = 'store',
                        default = 'chemical', help = 'chemical column in chemical-protein interaction file')
    parser.add_argument('-prot', dest = 'protein_col', action = 'store',
                        default = 'entrezid', help = 'chemical column in chemical-protein interaction file')
    parser.add_argument('-tmp', dest = 'tmp_dir', action = 'store', 
                        help = 'temporary folder', default='tmp_gsea')
    parser.add_argument('-cell', dest='cell_line', action='store', help='restrict to cell line, default = MCF7', 
                        default='MCF7')




    if len(sys.argv) <= 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    indir = os.path.abspath(args.indir) + '/'
    infile = args.interactome
    disease_genes_file = args.diseasegenes
    polyphenl_targets_file = args.targets
    ncpus = int(args.cores)
    n_random = int(args.nrandom)
    outdir_tmp_files = os.path.abspath(args.tmp_dir) + '/'
    sep = args.separator
    lcc = bool(args.lcc)
    columns = [args.source_node, args.target_node]
    chemical_col = args.chemical_col
    protein_col = args.protein_col
    outdir = os.path.abspath(args.outfolder) + '/'
    cell_line = args.cell_line
    header=True
    
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

        
        
    infiles = os.listdir(indir)
    infiles = [i for i in infiles if i.endswith('.gct')]
    
    for fi in infiles:
        
    
        ## parse gct file
        dt, attributes = parse_gct(indir + fi)

        ## select experimental instances to run gsea
        selected = attributes[(attributes.cell_id == cell_line)]
        timepoints = ['24 h', '6 h']
        treats = list(set(attributes['pert_idose']))
        samples = []
        for tp in timepoints:
            for treat in treats:
                buf = selected[(selected['pert_itime'] == tp) & (selected['pert_idose'] == treat)]
                ### if multiple replicates for the same treatment
                if buf.shape[0] == 1:
                    samples.append(buf['index'].iloc[0])
                if buf.shape[0] > 1:
                    sample = buf[buf['distil_cc_q75'] == buf['distil_cc_q75'].max()]['index'].iloc[0]
                    samples.append(sample)
                    
                    
        sample2attributes = {}
        for sample in samples:
            dic = {}
            buf = attributes[attributes['index'] == sample][['cell_id', 'name', 'pert_idose', 'pert_itime', 'index']]
            sample2attributes[sample] = buf.to_dict('records')[0]

        ## prepare tables
        cols = ['id'] + samples
        data = dt[cols].copy()
        data.loc[:,'id'] = data.loc[:,'id'].astype(int)
        data = data[data['id'].isin(list(G.nodes()))]
        for col in samples:
            data[col] = data[col].astype(float)
            
        logging.warning('Processing %d treatment conditions in: %s'%(len(samples),fi))
            
        ## run gsea for each experimental instance
        
        
        to_merge_diseases = []
        
        
        for e, test in enumerate(samples):

            dquery = data[[test, 'id']].copy()
            dquery = dquery.sort_values(by = test)
            dquery.columns = ['value', 'gene']

            dfs = []
            dic_null = {}
            for disease in tqdm.tqdm(disease2genes.keys()):
                res = defaultdict(dict)
                c = 0
                geneset = list(set(disease2genes[disease]) & set(G.nodes()))
                ns = list(set(dquery['gene']) - set(geneset))
                es = get_gsea_performance(geneset, ns, dquery)
                if not len(geneset) in dic_null.keys():
                    p = Pool(ncpus)
                    randomsets = np.repeat(len(geneset),n_random)
                    randomsets = [[i,dquery] for i in randomsets]
                    null = p.map(get_null, randomsets)
                    p.close()
                    dic_null[len(geneset)] = null
                res[c]['disease'] = disease
                res[c]['ES'] = es
                res[c]['Z_ES'] = (es - np.mean(dic_null[len(geneset)]))/np.std(dic_null[len(geneset)])
                res[c]['p_empirical'] = len([i for i in dic_null[len(geneset)] if i > es])/len(dic_null[len(geneset)])
                for att in sample2attributes[sample]:
                    res[c][att] = sample2attributes[sample][att]
                res = pd.DataFrame.from_dict(res, orient='index')
                ### save tmp
                dfs.append(res)
            
            res = pd.concat(dfs)
            to_merge_diseases.append(res)
        
        if len(to_merge_diseases) > 0:
            final = pd.concat(to_merge_diseases)
        else:
            final = pd.DataFrame()
        outname = fi.split('.')[0]
        outname = outname.replace(' ', '_')
        final.to_csv(outdir + '%s.csv'%(outname))