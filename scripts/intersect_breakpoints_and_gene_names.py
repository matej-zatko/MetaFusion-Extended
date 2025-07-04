#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import pygeneann_MetaFusion as pygeneann
import pybedtools.bedtool as bedtools
import itertools
import sequtils
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('cff', action='store', help='CFF file')

args = parser.parse_args()
cff=args.cff

#INTERSECT FUSIONS BY BREAKPOINTS
def intersect_fusions_by_breakpoints():
    lines=[line for line in open(cff, "r")]
    fusion=pygeneann.CffFusion(lines[0])
    header=fusion.zone1_attrs + fusion.zone2_attrs + fusion.zone3_attrs + fusion.zone4_attrs
    df_cff=pd.read_csv(cff, sep='\t', keep_default_na=False, index_col=False, names=header)
    
    # Combine sample name and chr to allow for only same-sample intersections
    df_cff['chr1'] = df_cff['chr1'] + "_" + df_cff['sample_name']
    df_cff['chr2'] = df_cff['chr2'] + "_" + df_cff['sample_name']
    
    #create BedTools object with appropriate column names
    print >> sys.stderr, "create BedTools object with appropriate column names"
    df_bed=df_cff[['chr1','pos1','pos1','chr2','pos2','pos2', 'fusion_id']]
    df_bed.columns=['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id']
    df_bed.loc[:,['pos1_2','pos2_2']] +=1
    df_bed=bedtools.BedTool.from_dataframe(df_bed)
    
    #Intersect fusions: NOTE: only keeps fusions that intersect
    print >> sys.stderr, "Intersect fusions: NOTE: rdn=False, keeps self-intersections"
    df_intersect=df_bed.pair_to_pair(df_bed, slop=100, rdn=False)
    df=df_intersect.to_dataframe(header=None).iloc[:,0:14]
    df.columns = ['chr1','pos1','pos1_2','chr2','pos2','pos2_2', 'fusion_id', 'chr1_1','pos1_1','pos1_2_1','chr2_1','pos2_1','pos2_2_1', 'fusion_id_lst'] 
    df=df[['fusion_id','fusion_id_lst']]
    #write paired F_IDs to tsv
    return df

df = intersect_fusions_by_breakpoints()
df.to_csv(sys.stdout,header=True,index=True, sep="\t")

# CLUSTER GENES

def intersect_fusions_by_genes(cff_file):
    fusion_dict = {}
    # cluster fusions by gene pairs AND sample name
    for line in open(cff_file, "r"):
        if line.startswith("#"):
            continue
        fusion = pygeneann.CffFusion(line)
        if fusion.t_gene1 == "NA" or fusion.t_gene2 == "NA":
            continue
        else:
            # Include sample_name in the key to ensure same-sample intersections only
            key = fusion.sample_name + "::" + ",".join(sorted([fusion.t_gene1 + "|" + fusion.chr1, fusion.t_gene2+ "|" + fusion.chr2])) 
            fusion_dict.setdefault(key, []).append(fusion.fusion_id)
    return fusion_dict

fusion_dict = intersect_fusions_by_genes(cff)

count = df.shape[0] + 1 
for key in fusion_dict.keys():
    lst=fusion_dict[key]
    edges=list(itertools.permutations(lst, 2))
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1

# CLUSTER by ENSEMBL IDs

def intersect_fusions_by_ids(cff_file):
    fusion_dict = {}
    # cluster fusions by Ensembl gene ID pairs AND sample name
    for line in open(cff_file, "r"):
        if line.startswith("#"):
            continue
        fusion = pygeneann.CffFusion(line)
        if fusion.t_gene_id1 == "NA" or fusion.t_gene_id2 == "NA":
            continue
        else:
            # Include sample_name in the clustering key
            key = fusion.sample_name + "::" + ",".join(sorted([fusion.t_gene_id1, fusion.t_gene_id2])) 
            fusion_dict.setdefault(key, []).append(fusion.fusion_id)
    return fusion_dict

fusion_dict = intersect_fusions_by_ids(cff)

for key in fusion_dict.keys():
    lst=fusion_dict[key]
    edges=list(itertools.permutations(lst, 2))
    for edge in edges:
        print("\t".join([str(count)] + list(edge)))
        count += 1