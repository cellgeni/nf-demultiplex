#!/usr/bin/env python

import pandas as pd
import argparse

my_parser = argparse.ArgumentParser()
my_parser.add_argument("-i", "--input", default=None, help="shared sample cluster tsv")
my_parser.add_argument("-s", "--sample", default=None, help="sample id")
args = my_parser.parse_args()

#Read in shared samples tsv of cluster comparions
adata = pd.read_csv(args.input, sep='\t')
#Add sample id to each cluster comparison
adata['shared_samples'] = args.sample
#Get list of cluster ids
cluster_ids = list(adata['experiment1_cluster'].unique())
#Create empty dataframe
adata_top_clusters = pd.DataFrame()
#Find the cluster comparison with the largest loss value for each cluster
for cluster in cluster_ids:
    top_cluster = adata.loc[adata["experiment1_cluster"] == cluster].sort_values(by=['loss'], ascending=False).iloc[[0]]
    #Generate dataframe with cluster comparions with largest loss function for each cluster id
    adata_top_clusters = pd.concat([adata_top_clusters, top_cluster])
#Save values to table
adata_top_clusters.to_csv(f"loss-tables/{args.sample}.tsv", sep='\t', index=False, header=False)
