#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from progressbar import ProgressBar
import networkx as nx
import os
import sys
import toolbox.databases_utils as databases_utils


if __name__ == '__main__':
    description = 'Expand associations in MeSH tree'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-m', required=True, dest='mesh',
                        action='store',help='MeSH File')
    parser.add_argument('-i', required=True, dest='infile',action='store',
                        help='association table')
    parser.add_argument('-o', dest='outfile', action='store',
                        help='outfile',
                        default= 'expanded.csv')
    
    if len(sys.argv) <= 1: 
        parser.print_help() 
        sys.exit(1) 
    else: 
        args = parser.parse_args()
    
    meshfile = args.mesh
    infile = args.infile
    outfile = args.outfile

    ## open mesh tree
    #meshfile = '/home/italodovalle/flavonoids/data/databases/mesh/mtrees2018.bin'
    g, disease2code, code2name = databases_utils.load_mesh_graph(meshfile)

    curation = {'MESH:D056486':['C06.552.195','C25.100.562','C25.723.260'],
                'MESH:C535575':['C04.557.470.200.400','C04.557.470.700.400'],
                'MESH:C538231':['C04.557.470.200.025'],
                'MESH:C537171':['C08.381.742'],
                'MESH:D024821':['C18.452.394.968.500.570','C18.452.625'],
                'MESH:D000860':['C23.888.852.079'],
                'MESH:C536572':['C15.378.071.085','C15.378.190.196','C15.378.190',
                               'C15.378.071.141.560','C15.378.190.625.460'],
                'MESH:C537948':['C10.574.500.550','C16.320.400.600','C16.320.565.398.641.509',
                 'C18.452.584.687.509','C18.452.648.398.641.509'],
                'MESH:C563161':['C12.777.419.331','C13.351.968.419.331','C14.907.489.631','C12.777.419.570',
                 'C13.351.968.419.570'],
                'MESH:C562839':['C04.588.894.797.520','C08.381.540','C08.785.520',
                               'C04.557.470.035.510','C04.557.470.660.510'],
                'MESH:C538339':['C04.557.470.200','C04.588.443.665.710.650',
                                'C07.550.350.650','C07.550.745.650',
                                'C09.647.710.650','C09.775.350.650','C09.775.549.650'],
                'MESH:C538339':['C04.588.443.665.710.650','C07.550.350.650',
                                'C07.550.745.650','C09.647.710.650',
                                'C09.775.350.650','C09.775.549.650',
                                'C04.557.470.200'],
                'MESH:C564616':['C10.228.140.252.190.530','C10.228.140.252.700.700',
                                'C10.228.854.787.875','C10.574.500.825.700',
                                'C10.597.350.090.500.530','C16.320.400.780.875'],
                'MESH:C536914':['C04.588.322.894','C04.588.443.915',
                                'C19.344.894','C19.874.788',
                               'C04.557.465.625.650.240','C04.557.470.200.025.370',
                                'C04.557.580.625.650.240']
                }
    
    ctd = pd.read_csv(infile,index_col = 0)
    print ('%d Chemical-Disease Associations'%ctd.shape[0])
    ctd['disease'] = [i.lower() for i in ctd['disease']]
    
    ## format cols in association file
    ### ['chemical', 'disease', 'DirectEvidence'] 
    ctd = ctd[~ctd.DirectEvidence.isnull()]
    ctd = ctd[ctd.DirectEvidence == 'therapeutic']
    
    ## all possible chemical disease associations
    table = defaultdict(dict)
    c = 0
    for i in set(ctd.chemical):
        for j in set(disease2code.keys()):
            table[c]['chemical'] = i
            table[c]['disease'] = j
            table[c]['DirectEvidence'] = float('nan')
            c = c + 1
    table = pd.DataFrame.from_dict(table, orient='index')
    
    merged = pd.merge(ctd, table,on = ['chemical', 'disease'],how='outer')
    merged['DirectEvidence'] = float('nan')
    merged.loc[(~merged['DirectEvidence_x'].isnull()), 'DirectEvidence'] = 'therapeutic'
    merged = merged[['chemical', 'disease', 'DirectEvidence']]
    
    print (merged[~merged.DirectEvidence.isnull()].shape[0], 'interactions DirectEvidence')
    print (len(set(merged[~merged.DirectEvidence.isnull()].chemical)), 'chemicals DirectEvidence')
    print (len(set(merged[~merged.DirectEvidence.isnull()].disease)), 'diseases DirectEvidence')
    
    print ('\nCreating chemical2code matrix\n')
    chemical2code = defaultdict(dict)
    pbar = ProgressBar()
    for i in pbar(merged.index):
        disease = merged.disease.loc[i]
        if disease in disease2code.keys():
            chemical = merged.chemical.loc[i]
            ## all possible codes for the disease
            codes = disease2code[disease]
            for code in codes:
                if merged.DirectEvidence.loc[i] == 'therapeutic':
                    chemical2code[chemical][code] = 1
                else:
                    chemical2code[chemical][code] = 0
    chemical2code = pd.DataFrame.from_dict(chemical2code)
    chemical2code = chemical2code.fillna(0)
    ## chemical2code is a matrix with all codes and all chemicals
    ## 1 if the chemical is associated with a disease code, zero otherwise
    
    print ('\nRetrieving Implicit from Explicit associations')
    merged['therapeutic'] = float('nan')
    diseases = set(merged.disease)
    pbar = ProgressBar()
    for disease in pbar(diseases):
        if disease in disease2code.keys():
            codes = disease2code[disease]
            ### of all codes associated with the disease, get those that have association with a chemical
            codes = list(set(codes) & set(chemical2code.index))
            ## chemicals explicitly associated with the disease
            ### x is matrix codes x chemicals
            x = chemical2code.loc[codes]
            x = x.T
            explicit = list(x[x>0].dropna(how='all').index)
            ## chemicals implicitly associated with the disease
            descendants = []
            for code in codes:
                descendants.extend(list(nx.descendants(g, code)))
            descendants = list(set(descendants))
            descendants = list(set(descendants) & set(chemical2code.index))
            x = chemical2code.loc[descendants]
            x = x.T
            implicit = list(x[x>0].dropna(how='all').index)
            evidence = list(set(explicit) | set(implicit))
            merged.loc[(merged.disease == disease) & (merged.chemical.isin(evidence)), 'therapeutic'] = 1
           
    print ('Writing output')
    merged.to_csv(outfile)
