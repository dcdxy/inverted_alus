#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def add_seq_length_col(df, category):
    new_df = df.copy()
    new_df = new_df[pd.to_numeric(new_df['mfe'], errors='coerce').notnull()]
    new_df['type'] = category
    new_df['mfe'] = new_df.mfe.astype(float)
    new_df['seq_length'] = new_df['stop'] - new_df['start']
    new_df['seq_corrected_mfe'] = (new_df['mfe'] / new_df['seq_length'])*100
    return new_df


def add_seq_length_col_updated(df, category):
    new_df = df.copy()
    new_df = new_df[pd.to_numeric(new_df['mfe'], errors='coerce').notnull()]
    new_df['type'] = category
    new_df['mfe'] = new_df.mfe.astype(float)
    new_df['alu1_start'] = pd.to_numeric(new_df['alu1_start'], errors='coerce')
    new_df['alu2_end'] = pd.to_numeric(new_df['alu2_end'], errors='coerce')
    new_df['seq_length'] = new_df['alu2_end'] - new_df['alu1_start']
    new_df['seq_corrected_mfe'] = (new_df['mfe'] / new_df['seq_length'])*100
    return new_df


def add_alu_cols(df):
    new_df = df.copy()
    #new_df = new_df[pd.to_numeric(new_df['mfe'], errors='coerce').notnull()]
    #new_df['mfe'] = new_df.mfe.astype(float)
    new_df['alu1_strand'] = new_df['alu1'].str[-1]
    new_df['alu1_id'] = new_df['alu1'].str[:-1]
    new_df['alu2_strand'] = new_df['alu2'].str[-1]
    new_df['alu2_id'] = new_df['alu2'].str[:-1]
    new_df["inversion"] = np.where(new_df['alu1_strand'] != new_df['alu2_strand'],
                                   "inverted", "non-inverted")
    return new_df


def merge_and_add_indicator_skip_const_dfs(df_skip, df_const):
    """
    df_skip: esdb_df_dct["skipped"] or hexevent_df_dct["skipped"]
    df_const: esdb_df_dct["constitutive"] or hexevent_df_dct["constitutive"]
    """
    keep_cols = ["chr", "start", "stop", "alu1", "alu2"]
    df_intersect_indicator = df_skip.merge(df_const, on=keep_cols, how="outer", indicator=True)
    return df_intersect_indicator #.query('_merge != "both"')


def create_mfe_dfs(w, cohens_d_dct, w_open_dct, w_open_dct_random, plot=False, print_stats=True, ax=None):
    from general_plots import generate_kdeplot_4sets, generate_kdeplot_5sets, generate_kdeplot_6sets
    from statistics_test import run_stats_tests_2groups
    cohens_d_dct[w] = {}
    df_const_noninv = add_seq_length_col_updated(w_open_dct[w]["constitutive_noninverted"],
                                                 "constitutive_noninverted")
    df_const_noninv = add_alu_cols(df_const_noninv)
    
    df_const_inv = add_seq_length_col_updated(w_open_dct[w]["constitutive_inverted"],
                                              "constitutive_inverted")
    df_const_inv = add_alu_cols(df_const_inv)
    
    df_skip_noninv = add_seq_length_col_updated(w_open_dct[w]["skippable_noninverted"],
                                                "skippable_noninverted")
    df_skip_noninv = add_alu_cols(df_skip_noninv)
    
    df_skip_inv = add_seq_length_col_updated(w_open_dct[w]["skippable_inverted"],
                                             "skippable_inverted")
    df_skip_inv = add_alu_cols(df_skip_inv)

    ##### new 24-04-11: Random
    rand_df_const_noninv = add_seq_length_col_updated(w_open_dct_random[w]["constitutive_noninverted"],
                                                      "constitutive_noninverted")
    rand_df_const_noninv = add_alu_cols(rand_df_const_noninv)
    rand_df_const_inv = add_seq_length_col_updated(w_open_dct_random[w]["constitutive_inverted"],
                                                   "constitutive_inverted")
    rand_df_const_inv = add_alu_cols(rand_df_const_inv)
    rand_df_skip_noninv = add_seq_length_col_updated(w_open_dct_random[w]["skippable_noninverted"],
                                                     "skippable_noninverted")
    rand_df_skip_noninv = add_alu_cols(rand_df_skip_noninv)
    rand_df_skip_inv = add_seq_length_col_updated(w_open_dct_random[w]["skippable_inverted"],
                                                  "skippable_inverted")
    rand_df_skip_inv = add_alu_cols(rand_df_skip_inv)

    rand_df_skip = pd.concat([rand_df_skip_inv, rand_df_skip_noninv], ignore_index=True)
    rand_df_const = pd.concat([rand_df_const_inv, rand_df_const_noninv], ignore_index=True)
    rand_df_inv = pd.concat([rand_df_skip_inv, rand_df_const_inv], ignore_index=True)
    rand_df_noninv = pd.concat([rand_df_skip_noninv, rand_df_const_noninv], ignore_index=True)
    #####
    
    # merge
    df_skip = pd.concat([df_skip_inv, df_skip_noninv], ignore_index=True)
    df_const = pd.concat([df_const_inv, df_const_noninv], ignore_index=True)

    df_inv = pd.concat([df_skip_inv, df_const_inv], ignore_index=True)
    df_noninv = pd.concat([df_skip_noninv, df_const_noninv], ignore_index=True)

    # Updated July 26 for manuscript figure
    if plot:
        #print(f"{w}:")
        #generate_kdeplot_4sets(df_skip=df_skip, df_const=df_const)
        #generate_kdeplot_4sets(df_skip=rand_df_skip, df_const=rand_df_const)
        #generate_kdeplot_6sets(df_skip=df_skip, df_const=df_const, 
        #                       df_randskip=rand_df_skip, df_randconst=rand_df_const, title="")
        generate_kdeplot_5sets(df_inv=df_inv, df_noninv=df_noninv, 
                               df_randinv=rand_df_inv, df_randnoninv=rand_df_noninv, title="", ax=ax) 
        # no easy way to change colours in the current setup above for 5sets, since hue is skip/const

    # stats
    if print_stats:
        #print(f"Skippable vs Constitutive")
        cohens_d_dct[w]["skippable_vs_constitutive"] = run_stats_tests_2groups(df_skip, df_const)
        #print(f"Inverted vs Non-inverted")
        cohens_d_dct[w]["inverted_vs_noninverted"] = run_stats_tests_2groups(df_inv, df_noninv)
        #print(f"Skippable_inverted vs Skippable_non-inverted")
        cohens_d_dct[w]["skippable_inverted_vs_skippable_noninverted"] = (
            run_stats_tests_2groups(df_skip_inv, df_skip_noninv))
        #print(f"Constitutive_inverted vs Constitutive_non-inverted")
        cohens_d_dct[w]["constitutive_inverted_vs_constitutive_noninverted"] = (
            run_stats_tests_2groups(df_const_inv, df_const_noninv))
        #print(f"Skippable_inverted vs Constitutive_inverted")
        cohens_d_dct[w]["skippable_inverted_vs_constitutive_inverted"] = (
            run_stats_tests_2groups(df_skip_inv, df_const_inv))
        #print(f"Skippable_non-inverted vs Constitutive_non-inverted")
        cohens_d_dct[w]["skippable_noninverted_vs_constitutive_noninverted"] = (
            run_stats_tests_2groups(df_skip_noninv, df_const_noninv))
        #print(f"Skippable_inverted vs Random")
        cohens_d_dct[w]["skippable_inverted_vs_random"] = (
            run_stats_tests_2groups(df_skip_inv, rand_df_const_noninv))
        #print(f"Skippable_noninverted vs Random")
        cohens_d_dct[w]["skippable_noninverted_vs_random"] = (
            run_stats_tests_2groups(df_skip_noninv, rand_df_const_noninv))
        #print(f"Constitutive_inverted vs Random")
        cohens_d_dct[w]["constitutive_inverted_vs_random"] = (
            run_stats_tests_2groups(df_const_inv, rand_df_const_noninv))
        #print(f"Constitutive_noninverted vs Random")
        cohens_d_dct[w]["constitutive_noninverted_vs_random"] = (
            run_stats_tests_2groups(df_const_noninv, rand_df_const_noninv))
    
    mfe_dct = {"w": w,
               "df_skip": df_skip, 
               "df_const": df_const, 
               "df_inv": df_inv, 
               "df_noninv": df_noninv}
    return cohens_d_dct, mfe_dct


def create_mfe_dfs_suppl(w, w_open_dct, ax, plot=False, print_stats=True, legend=True):
    from general_plots import generate_kdeplot_4sets
    from statistics_test import run_stats_tests_2groups
    cohens_d_dct = {}
    cohens_d_dct[w] = {}
    df_const_noninv = (
        add_seq_length_col_updated(w_open_dct[w]["constitutive_noninverted"],
                                   "constitutive_noninverted"))
    df_const_noninv = add_alu_cols(df_const_noninv)
    df_const_inv = (
        add_seq_length_col_updated(w_open_dct[w]["constitutive_inverted"],
                                   "constitutive_inverted"))
    df_const_inv = add_alu_cols(df_const_inv)
    df_skip_noninv = (
        add_seq_length_col_updated(w_open_dct[w]["skippable_noninverted"],
                                   "skippable_noninverted"))
    df_skip_noninv = add_alu_cols(df_skip_noninv)
    df_skip_inv = (
        add_seq_length_col_updated(w_open_dct[w]["skippable_inverted"],
                                   "skippable_inverted"))
    df_skip_inv = add_alu_cols(df_skip_inv)
    #####
    
    # merge
    df_skip = pd.concat([df_skip_inv, df_skip_noninv], ignore_index=True)
    df_const = pd.concat([df_const_inv, df_const_noninv], ignore_index=True)

    df_inv = pd.concat([df_skip_inv, df_const_inv], ignore_index=True)
    df_noninv = pd.concat([df_skip_noninv, df_const_noninv], ignore_index=True)

    # plot
    # Updated July 29 for manuscript figure
    if plot:
        print(f"{w}:")
        generate_kdeplot_4sets(df_inv=df_inv, df_noninv=df_noninv, 
                               ax=ax, show_legend=legend)

    # stats
    if print_stats:
        print(f"Skippable vs Constitutive")
        cohens_d_dct[w]["skippable_vs_constitutive"] = run_stats_tests_2groups(df_skip, df_const)
        print(f"Inverted vs Non-inverted")
        cohens_d_dct[w]["inverted_vs_noninverted"] = run_stats_tests_2groups(df_inv, df_noninv)
        print(f"Skippable_inverted vs Skippable_non-inverted")
        cohens_d_dct[w]["skippable_inverted_vs_skippable_noninverted"] = (
            run_stats_tests_2groups(df_skip_inv, df_skip_noninv))
        print(f"Constitutive_inverted vs Constitutive_non-inverted")
        cohens_d_dct[w]["constitutive_inverted_vs_constitutive_noninverted"] = (
            run_stats_tests_2groups(df_const_inv, df_const_noninv))
        print(f"Skippable_inverted vs Constitutive_inverted")
        cohens_d_dct[w]["skippable_inverted_vs_constitutive_inverted"] = (
            run_stats_tests_2groups(df_skip_inv, df_const_inv))
        print(f"Skippable_non-inverted vs Constitutive_non-inverted")
        cohens_d_dct[w]["skippable_noninverted_vs_constitutive_noninverted"] = (
            run_stats_tests_2groups(df_skip_noninv, df_const_noninv))
    
    mfe_dct = {"w": w,
               "df_skip": df_skip, 
               "df_const": df_const, 
               "df_inv": df_inv, 
               "df_noninv": df_noninv}
    return cohens_d_dct, mfe_dct
