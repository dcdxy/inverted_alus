#!/usr/bin/env python
import pandas as pd
from itertools import combinations

def generate_dfs(skipped_exon_df, control_df, ir=True, header="bedops", exon_db=False, verbose=False):
    """
    Window size is fixed to +/- 5000 bp.
    skipped_exon_df: previously closest_exon_df
    control_df: generally closest_const_exon_df which are the constitutive exons
    """
    #This section below is not included in the output
    w = 5000
    window_df = skipped_exon_df[skipped_exon_df['dist'].between(-w, w)].dropna(how='any')
    window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})
    control_window_df = control_df[control_df['dist'].between(-w, w)].dropna(how='any')  # added 230428
    control_window_df = control_window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})  # added 230428

    if verbose:
        window_df_counts = window_df['alu_subfamily'].value_counts()
        control_window_df_counts = control_window_df['alu_subfamily'].value_counts()
        print(window_df_counts.head())
        print(control_window_df_counts.head())

    # Create inputs for main function (these do not depend on the window_dfs generated above)
    flank_df = create_flanking_df(w, skipped_exon_df, header=header, exon_db=exon_db)
    control_flank_df = create_flanking_df(w, control_df, header=header, exon_db=exon_db)
    
    if ir:
        ir_flank_df = filter_inverted_pair(flank_df.copy())
        control_ir_flank_df = filter_inverted_pair(control_flank_df.copy())
        return flank_df, control_flank_df, ir_flank_df, control_ir_flank_df
    else:
        return flank_df, control_flank_df

def create_flanking_df(window, exon_df, header="bedops", exon_db=False, discrete=False): # 24-04-23 added discrete
    """
    From BEDOPS closest-feature outputs (with alu file first and exon file second),
    find Alus that are both within a certain distance to exons and that form flanking pairs.
    
    -window: length upstream and downstream of exon to search within (usually 5000 bp)
    -exon_df: closest_exon_df or closest_const_exon_df
    """
    # Filter on +/- window around exons
    if not discrete:
        window_df = exon_df[exon_df['dist'].between(-window, window)].dropna(how='any')  # dropna doesn't change anything once window is applied
    elif discrete:
        range_condition = (exon_df['dist'].between(-window, -window + 500, inclusive='neither') |
                           exon_df['dist'].between(window - 500, window, inclusive='neither'))
        window_df = exon_df[range_condition].dropna(how='any')
    window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})
    display(window_df.head())

    # All possible combinations of alu pairings
    if header == "bedops":
        keep_exon_cols = ["exon_chr", "exon_start", ]
    if header=="bedops":
        keep_exon_cols = ['exon_chr', 'exon_start','exon_end','exon_gene', 'exon_score', 'exon_strand']
    elif header=="exonskipdb":
        keep_exon_cols = ['exon_chr','exon_start','exon_end', 'exon_gene', 'exon_strand', 'exon_type']
    if exon_db: # Added 23-04-28
        keep_exon_cols = keep_exon_cols + ["exon_db"]
    
    grouped_df = (
        window_df
        .groupby(keep_exon_cols, sort=False)[['alu_strand', 'idx']]
        .apply(lambda x: list(combinations(x.values,2)))
        .apply(pd.Series, dtype="object")
        .stack()
        .reset_index(name='strand_pairs'))  # dtype='float'
    grouped_df[['upstream', 'downstream']] = pd.DataFrame(grouped_df['strand_pairs'].tolist(), columns=['upstream', 'downstream'])
    grouped_df[['upstream_strand', 'upstream_id']] = pd.DataFrame(grouped_df['upstream'].tolist(), columns=['upstream_strand', 'upstream_idx'])
    grouped_df[['downstream_strand', 'downstream_id']] = pd.DataFrame(grouped_df['downstream'].tolist(), columns=['downstream_strand', 'downstream_idx'])
    grouped_df = grouped_df.astype({'exon_start': 'int', 'exon_end': 'int'})
    
    # filter on pairs that also flank exons
    cols = ['alu_chr', 'alu_start', 'alu_end', 'alu_subfamily', 'alu_sw_score', 'alu_strand', 
            'alu_percent_substitution', 'alu_percent_deleted', 'alu_percent_inserted', 'alu_num_bases_past_end',
            'alu_id', 'dist']
    rename_upstream_cols = dict(zip(cols,["upstream_{}".format(i) for i in cols]))
    rename_downstream_cols = dict(zip(cols,["downstream_{}".format(i) for i in cols]))
    flank_df = grouped_df.copy()
    flank_df = flank_df.merge(window_df[['idx'] + cols].rename(columns={'idx': 'upstream_id'} | rename_upstream_cols),
                              on='upstream_id', how='left')
    flank_df = flank_df.merge(window_df[['idx'] + cols].rename(columns={'idx': 'downstream_id'} | rename_downstream_cols),
                              on='downstream_id', how='left')
    flank_df = flank_df[(flank_df['upstream_dist']>0) & (flank_df['downstream_dist']<0)] # but only correct orientation if same BEDOPS setup
    return flank_df


def filter_inverted_pair(flanking_df):
    """Filter rows where the upstream strand does not match the downstream strand (ie. +/- or -/+).
    """
    return flanking_df[(flanking_df.upstream_strand != flanking_df.downstream_strand)]


def filter_non_inverted_pair(flanking_df):
    """Filter rows where the upstream strand matches the downstream strand (ie. +/+ or -/-).
    """
    return flanking_df[(flanking_df.upstream_strand == flanking_df.downstream_strand)]


def groupby_exons(exons_df: pd.DataFrame):
    exon_cols = ["exon_chr", "exon_start", "exon_end", "exon_gene"]
    exons_df = exons_df.copy()
    exons_df = (exons_df
                .groupby(by=exon_cols, group_keys=True)
                .count()
                .reset_index())
    exons_df[["exon_start", "exon_end"]] = exons_df[["exon_start", "exon_end"]].astype(int)
    return exons_df


def create_df_subsets(full_df: pd.DataFrame):
    """Split a dataframe according to its 'type' and 'inversion' columns.
    """
    skipped_inverted = full_df[(full_df['type'] == 'skipped') & (full_df['inversion'] == 'inverted')]
    skipped_non_inverted = full_df[(full_df['type'] == 'skipped') & (full_df['inversion'] == 'non-inverted')]
    constitutive_inverted = full_df[(full_df['type'] == 'constitutive') & (full_df['inversion'] == 'inverted')]
    constitutive_non_inverted = full_df[(full_df['type'] == 'constitutive') & (full_df['inversion'] == 'non-inverted')]
    skipped = full_df[full_df['type'] == 'skipped']
    constitutive = full_df[full_df['type'] == 'constitutive']
    inverted = full_df[full_df['inversion'] == 'inverted']
    non_inverted = full_df[full_df['inversion'] == 'non-inverted']
    return {"skipped_inverted": skipped_inverted,
            "skipped_non_inverted": skipped_non_inverted,
            "constitutive_inverted": constitutive_inverted,
            "constitutive_non_inverted": constitutive_non_inverted,
            "skipped": skipped,
            "constitutive": constitutive,
            "inverted": inverted,
            "non_inverted": non_inverted}
    #return (skipped_inverted, skipped_non_inverted,
    #        constitutive_inverted, constitutive_non_inverted,
    #        skipped, constitutive,
    #        inverted, non_inverted)