import os
import pandas as pd
import numpy as np
import argparse
from itertools import combinations

def load_exon_list(i):

    processdir = "."

    upstream_alus = os.path.join(processdir, f"temp_{i}_upstream.tsv")
    downstream_alus = os.path.join(processdir, f"temp_{i}_downstream.tsv")

    bed6_cols = ['chr', 'start', 'end']

    alu_bed_cols = ['chr', 'start', 'end', 'subfamily', 'strand', 'id']
    closest_circRNA_header = ["alu_{}".format(i) for i in alu_bed_cols] + \
                            ["circRNA_{}".format(i) for i in bed6_cols] + ["dist"]
    cols = [0,1,2,3,5,14,15,16,17,18]


    up_df = pd.read_csv(upstream_alus, sep='[\t|]', engine='python', header=None, names=closest_circRNA_header, usecols=cols)
    down_df = pd.read_csv(downstream_alus, sep='[\t|]', engine='python', header=None, names=closest_circRNA_header, usecols=cols)

    up_df['alu_type'] = "upstream"
    down_df['alu_type'] = "downstream"

    up_df = up_df.drop(up_df[up_df['circRNA_chr'] == '.'].index) #delete unassigned alus
    down_df = down_df.drop(down_df[down_df['circRNA_chr'] == '.'].index) #delete unassigned alus

    merged_df = pd.concat([up_df, down_df])
    merged_df = merged_df.sort_values(["alu_chr","alu_start"]).reset_index(drop=True)
    merged_df = merged_df.drop_duplicates()
    
    hominoid_specific_Alu_df = pd.read_csv("liftover_total.tsv", header = None, names = ["alu_id"], usecols = [4], sep = "\t")
    merged_df["hominoid"] = "no"
    merged_df.loc[merged_df["alu_id"].isin(hominoid_specific_Alu_df["alu_id"]), "hominoid"] = "hominoid"
    return merged_df

def create_ir_df(window, exon_df):
    window_df = exon_df[exon_df['dist'].between(-window,window)].dropna(how='any')  # dropna doesn't change anything once window is applied

    keep_circRNA_cols = ['circRNA_chr','circRNA_start','circRNA_end']
    grouped_df = window_df.groupby(keep_circRNA_cols, sort=False)[['alu_strand', 'alu_id', 'alu_start', 'alu_end']].apply(lambda x: list(combinations(x.values,2))).apply(pd.Series).stack().reset_index(name='strand_pairs')
    grouped_df[['pair_1', 'pair_2']] = pd.DataFrame(grouped_df['strand_pairs'].tolist(), columns=['pair_1', 'pair_2'])

    grouped_df[['pair_1_strand', 'pair_1_id', 'pair_1_start', 'pair_1_end']] = pd.DataFrame(grouped_df['pair_1'].tolist(), columns=['pair_1_strand', 'pair_1_id', 'pair_1_start', 'pair_1_end'])
    grouped_df[['pair_2_strand', 'pair_2_id', 'pair_2_start', 'pair_2_end']] = pd.DataFrame(grouped_df['pair_2'].tolist(), columns=['pair_2_strand', 'pair_2_id', 'pair_2_start', 'pair_2_end'])
    grouped_df = grouped_df.astype({'circRNA_start': 'int', 'circRNA_end': 'int'})

    flank_df = grouped_df.copy()

    cols = ['hominoid', 'dist']

    flank_df = flank_df.merge(window_df[['alu_id'] + cols]
                        .rename(columns={'alu_id': 'pair_1_id', 'hominoid': 'upstream_hominoid', 'dist': 'upstream_dist'}),
                        on='pair_1_id', how='left')
    flank_df = flank_df.merge(window_df[['alu_id'] + cols]
                        .rename(columns={'alu_id': 'pair_2_id', 'hominoid': 'downstream_hominoid', 'dist': 'downstream_dist'}),
                        on='pair_2_id', how='left')

    flank_df = flank_df[(flank_df['upstream_dist']>0) & (flank_df['downstream_dist']<0)]
    ir_df = flank_df[(flank_df.pair_1_strand != flank_df.pair_2_strand)]
    ir_df = ir_df.drop_duplicates(['circRNA_chr','circRNA_start','circRNA_end','pair_1_id','pair_2_id'])
    
    no_no_ir_df = ir_df[(ir_df["upstream_hominoid"] == "no") & (ir_df["downstream_hominoid"] == "no")]
    no_hominoid_ir_df = ir_df[(ir_df["upstream_hominoid"] != ir_df["downstream_hominoid"])]
    hominoid_hominoid_ir_df = ir_df[(ir_df["upstream_hominoid"] == "hominoid") & (ir_df["downstream_hominoid"] == "hominoid")]

    return len(ir_df.groupby(['circRNA_chr','circRNA_start','circRNA_end'])), len(no_no_ir_df.groupby(['circRNA_chr','circRNA_start','circRNA_end'])), len(no_hominoid_ir_df.groupby(['circRNA_chr','circRNA_start','circRNA_end'])), len(hominoid_hominoid_ir_df.groupby(['circRNA_chr','circRNA_start','circRNA_end']))

parser = argparse.ArgumentParser()
parser.add_argument("-i", type = str)
parser.add_argument("-w", type = int)
args = parser.parse_args()

a,b,c,d = create_ir_df(args.w, load_exon_list(args.i))
print(f"{a}\t{b}\t{c}\t{d}")