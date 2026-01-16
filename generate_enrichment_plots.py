#!/usr/bin/env python
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from itertools import combinations
from prettytable import PrettyTable
from pprint import pprint
from scipy.stats import fisher_exact, chi2_contingency, power_divergence
from statsmodels.stats.multitest import multipletests

from generate_exon_alu_df import (create_flanking_df, filter_inverted_pair, 
                                  filter_non_inverted_pair, groupby_exons)
from general_plots import remove_shared_exons


def generate_dfs_windows(window, skipped_exon_df, control_df, ir=True, header="bedops"):
    """
    Modified from generate_exon_alu_df.generate_dfs(). Here, the window size is NOT fixed.
    
    If "ir==True", returns a dictionary which includes window boundary dataframe.
    Otherwise, returns the same as generate_dfs() ie. flank_df and control_flank_df,
    which are not based on the window parameter.

    skipped_exon_df: previously closest_exon_df
    control_df: generally closest_const_exon_df which are the constitutive exons
    """
    w = window
    window_df = skipped_exon_df[skipped_exon_df['dist'].between(-w, w)].dropna(how='any')
    window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})
    
    flank_df = create_flanking_df(w, skipped_exon_df, header=header)
    control_flank_df = create_flanking_df(w, control_df, header=header)
    if ir:
        ir_flank_df = filter_inverted_pair(flank_df.copy())
        control_ir_flank_df = filter_inverted_pair(control_flank_df.copy())
        return {"skippable": flank_df, 
                "constitutive": control_flank_df,
                "skippable_inverted": ir_flank_df,
                "constitutive_inverted": control_ir_flank_df,
                "window": window_df}
    else:
        return flank_df, control_flank_df


def generate_dfs_windows_subset(window, 
                                hexevent_skipped_exon_df, 
                                exonskipdb_skipped_exon_df,
                                control_df,
                                skippable_set, constitutive_set,
                                ir=True, header="bedops",
                                report_counts_only=False,
                                discrete=False):
    """
    Added 24-01-09 to filter out skipped/constitutive overlapping exons.
    Skippable_set and constitutive_set are the exon sets that we have set.
    Exon_db should be either hexevent or exonskipdb for the skipped_exon_df.
    
    Modified from generate_exon_alu_df.generate_dfs(). Here, the window size is NOT fixed.
    
    If "ir==True", returns a dictionary which includes window boundary dataframe.
    Otherwise, returns the same as generate_dfs() ie. flank_df and control_flank_df,
    which are not based on the window parameter.

    skipped_exon_df: previously closest_exon_df
    control_df: generally closest_const_exon_df which are the constitutive exons
    """
    w = window
    #window_df = skipped_exon_df[skipped_exon_df['dist'].between(-w, w)].dropna(how='any')
    #window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})
    
    # Skippable set can be from "hexevent" or "exonskipdb"
    hexevent_skipped_exon_df_subset = (
        remove_shared_exons(
            hexevent_skipped_exon_df, # bedops_skippable
            skippable_set,
            exon_type="skippable",
            exon_db="hexevent",
            verbose=True
        )
    )
    exonskipdb_skipped_exon_df_subset = (
        remove_shared_exons(
            exonskipdb_skipped_exon_df, # bedops_exonskipdb_skippable
            skippable_set,
            exon_type="skippable",
            exon_db="exonskipdb",
            verbose=True
        )
    )
    # Control set will always be from "hexevent"
    control_df_subset = (
        remove_shared_exons(
            control_df, # bedops_constitutive
            constitutive_set,
            exon_type="constitutive", 
            exon_db="hexevent",
            verbose=True
        )
    )
    if report_counts_only:
        return
    
    def create_window_df(skipped_exon_df_subset):
        window_df_subset = (
            skipped_exon_df_subset[skipped_exon_df_subset['dist'].between(-w, w)]
            .dropna(how='any')
        )
        window_df_subset = (
            window_df_subset.rename_axis('idx').reset_index().astype({'dist': 'int'})
        )
        return window_df_subset
    
    
    def merge_skippable_dfs(hexevent_df, exonskipdb_df, window=False):
        key_cols = ["exon_chr", "exon_start", "exon_end"]
        if not window:
            key_cols += ["upstream_id", "downstream_id"]
        merged_df = pd.concat([hexevent_df.drop("exon_score", axis=1), exonskipdb_df])
        merged_df["key"] = merged_df[key_cols].apply(lambda row: '-'.join(row.values.astype(str)), axis=1)
        merged_df.drop_duplicates(subset=["key"], keep="first", inplace = True)
        merged_df["exon_score"] = 1
        return merged_df
    
    
    hexevent_window_df_subset = create_window_df(hexevent_skipped_exon_df_subset)
    exonskipdb_window_df_subset = create_window_df(exonskipdb_skipped_exon_df_subset)
    window_df_subset = merge_skippable_dfs(hexevent_window_df_subset, exonskipdb_window_df_subset, window=True)
    skippable_subset = merge_skippable_dfs(hexevent_skipped_exon_df_subset, exonskipdb_skipped_exon_df_subset, window=True)
    print(len(skippable_subset))
    
    hexevent_flank_df = create_flanking_df(w, hexevent_skipped_exon_df_subset, header="bedops", exon_db=True, discrete=discrete) # skipped_exon_df
    exonskipdb_flank_df = create_flanking_df(w, exonskipdb_skipped_exon_df_subset, header="exonskipdb", exon_db=True, discrete=discrete) # skipped_exon_df
    control_flank_df = create_flanking_df(w, control_df_subset, header="bedops", exon_db=True, discrete=discrete) # control_df
    flank_df = merge_skippable_dfs(hexevent_flank_df, exonskipdb_flank_df)
    
    if ir:
        ir_flank_df = filter_inverted_pair(flank_df.copy())
        nonir_flank_df = filter_non_inverted_pair(flank_df.copy())
        control_ir_flank_df = filter_inverted_pair(control_flank_df.copy())
        control_nonir_flank_df = filter_non_inverted_pair(control_flank_df.copy())
        
        return {"skippable": flank_df, 
                "constitutive": control_flank_df,
                "skippable_inverted": ir_flank_df,
                "constitutive_inverted": control_ir_flank_df,
                "skippable_noninverted": nonir_flank_df,
                "constitutive_noninverted": control_nonir_flank_df,
                "skippable_hexevent": hexevent_flank_df,
                "skippable_exonskipdb": exonskipdb_flank_df,
                "skippable_hexevent_single_preflank": hexevent_skipped_exon_df_subset,
                "skippable_exonskipdb_single_preflank": exonskipdb_skipped_exon_df_subset,
                "skippable_merged_single_preflank": skippable_subset,
                "constitutive_single_preflank": control_df_subset,
                "window": window_df_subset}
    else:
        return flank_df, control_flank_df


def generate_contingency_table(r1c1, r1c2, r2c1, r2c2, r1_label, r2_label, 
                               verbose=True, alu=None, chi=None):
    # r: row, c: col
    table = np.array([[r1c1, r1c2],
                      [r2c1, r2c2]])
    t = PrettyTable(['', 'Flanking skipped exon', 'Flanking constitutive exon'])
    t.add_row([r1_label, "{:,}".format(table[0,0]), "{:,}".format(table[0,1])])
    t.add_row([r2_label, "{:,}".format(table[1,0]), "{:,}".format(table[1,1])])
    if verbose:
        print(t)
    stats_dct = test_table_significance(table, alu, chi)
    return stats_dct


def test_table_significance(table, alu=None, chi=True):
    #oddsr, p = fisher_exact(table, alternative='two-sided') # two-sided, greater, less
    oddsr_less, p_less = fisher_exact(table, alternative='less')
    oddsr_greater, p_greater = fisher_exact(table, alternative='greater')
    #chi2, chi_p, chi_dof, chi_expctd = chi2_contingency(table)
    stats_dct = {"fisher_greater_p": p_greater,
                 "fisher_less_p": p_less,
                 "fisher_odds_ratio_greater": oddsr_greater,
                 "fisher_odds_ratio_less": oddsr_less,
                 #"fisher_p_value": p,
                 #"fisher_odds_ratio": oddsr,
                 #"chi2": chi2,
                 #"chi_p_value": chi_p,
                 #"dof": chi_dof,
                 #"expected": np.around(chi_expctd),
                 "table": table,
                 "r1c1": table[0,0],
                 "r1c2": table[0,1],
                 "r2c1": table[1,0],
                 "r2c2": table[1,1]}  # added r2c1, r2c2
    if chi:
        chi2, chi_p, chi_dof, chi_expctd = chi2_contingency(table)
        stats_dct["chi_p_value"] = chi_p
        stats_dct["dof"] = chi_dof
        stats_dct["expected"] = np.around(chi_expctd)
        stats_dct["chi2"] = chi2
    if alu:
        stats_dct['alu_subfamily'] = alu
    return stats_dct


def test_table_across_windows(flank_df, control_flank_df, 
                              ir_flank_df, control_ir_flank_df):
    """
    Modified from test_table() in 02_ir_alus_plots_windows.ipynb
    
    # Specific case: alu=False, pair=False, inverted=False, comb=False (IR vs non-IR; all subfamilies together)
    """
    # Also generate dictionaries?
    r1c1 = len(ir_flank_df)
    r1c2 = len(control_ir_flank_df)
    
    #r2c1 = len(flank_df) - len(ir_flank_df)
    #r2c2 = len(control_flank_df) - len(control_ir_flank_df)
    
    non_ir_flank_df = filter_non_inverted_pair(flank_df.copy())
    control_non_ir_flank_df = filter_non_inverted_pair(control_flank_df.copy())
    r2c1 = len(non_ir_flank_df) # number of Alu pairs with at least 1 flanking non-IR Alu pair
    r2c2 = len(control_non_ir_flank_df) # number of Alu pairs with at least 1 flanking non-IR Alu pair
    
    result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                            r1_label="Inverted pair",
                                            r2_label="Non-inverted pair",
                                            verbose=True)
    return result_dct


def test_table_across_windows_exon(flank_df, control_flank_df, 
                                   ir_flank_df, control_ir_flank_df):
    """
    Modified from test_table() in 02_ir_alus_plots_windows.ipynb
    
    This function differs from test_table_across_windows() in that it counts the number of 
    exons rather than the number of Alu-Alu pairs.
    """
    r1c1 = len(groupby_exons(ir_flank_df)) # number of exons with at least 1 flanking IR Alu pair
    r1c2 = len(groupby_exons(control_ir_flank_df)) # number of exons with at least 1 flanking IR Alu pair
    
    # change this -- don't subtract from total, instead calculate separately (so there will be double-counting
    # between skipped and constitutive sets
    #r2c1 = len(groupby_exons(flank_df)) - len(groupby_exons(ir_flank_df)) # previously; biased towards IR Alus
    #r2c2 = len(groupby_exons(control_flank_df)) - len(groupby_exons(control_ir_flank_df))
    
    # Added on 23-01-19
    non_ir_flank_df = filter_non_inverted_pair(flank_df.copy())
    control_non_ir_flank_df = filter_non_inverted_pair(control_flank_df.copy())
    r2c1 = len(groupby_exons(non_ir_flank_df)) # number of exons with at least 1 flanking non-IR Alu pair
    r2c2 = len(groupby_exons(control_non_ir_flank_df)) # number of exons with at least 1 flanking non-IR Alu pair
    
    result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                            r1_label="Inverted pair",
                                            r2_label="Non-inverted pair",
                                            verbose=True)
    return result_dct


def bar_plot_subfamilies(stats_df, window, r1c1_label, r1c2_label, title):
    fig, ax = plt.subplots()
    ax.set_xlabel("Alu subfamily")
    ax.set_ylabel("Frequency")
    
    stats_df = stats_df.sort_values(by=['r1c2'], ascending=False)
    stats_idx = list(stats_df.set_index('alu_subfamily').index)
    
    #bonferroni_reject_list = stats_df[stats_df['bonferroni_reject']]['alu_subfamily'].tolist()
    bonferroni_reject_greater = stats_df[stats_df['bonferroni_reject_greater']]['alu_subfamily'].tolist()
    bonferroni_reject_less = stats_df[stats_df['bonferroni_reject_less']]['alu_subfamily'].tolist()
    bonferroni_reject = stats_df[stats_df['bonferroni_reject_less'] | stats_df['bonferroni_reject_greater']]['alu_subfamily'].tolist()
    normal_color = 'steelblue'
    greater_color = 'crimson'
    less_color = 'green'
    color_dict = {x: greater_color for x in bonferroni_reject_greater} | {x: less_color for x in bonferroni_reject_less} ##FF0000,#0000FF

    
    # (1) constitutive exons (background)
    ax1 = stats_df.set_index('alu_subfamily')['r1c2'].plot(kind='bar',
                                                           log=False, 
                                                           figsize=(22,6), #14,5
                                                           color='grey',
                                                           label=r1c2_label)
    # (2) skipped exons (foreground)
    ax2 = stats_df.set_index('alu_subfamily')['r1c1'].plot(kind='bar',
                                                           log=False,
                                                           figsize=(22,6),
                                                           color=[color_dict.get(x, normal_color) for x in stats_idx],
                                                           label=r1c1_label)
    for idx,p in enumerate(ax1.patches):
        #if idx in [window_idx.index(x) for x in bonferroni_reject_list]:
        if idx in [stats_idx.index(x) for x in bonferroni_reject_greater]:
            ax.annotate('{:0.2e}'.format(stats_df[stats_df['alu_subfamily'] == stats_idx[idx]]["bonferroni_p_greater"].tolist()[0]),
                        (p.get_x() * 1.005, (p.get_height()+1) * 1.02), rotation=45)  #str(p.get_height()) #+1, *1.005
        if idx in [stats_idx.index(x) for x in bonferroni_reject_less]:
            ax.annotate('{:0.2e}'.format(stats_df[stats_df['alu_subfamily'] == stats_idx[idx]]["bonferroni_p_less"].tolist()[0]),
                        (p.get_x() * 1.005, (p.get_height()+1) * 1.02), rotation=45)  #str(p.get_height()) #+1, *1.005
    
    #reference: https://stackoverflow.com/questions/39500265/
    leg = ax.legend()
    leg.legendHandles[1].set_color(normal_color)
    handles, labels = ax.get_legend_handles_labels()
    patch = mpatches.Patch(color=greater_color, label='Reject null hypothesis for FE test alternative "greater"')
    patch2 = mpatches.Patch(color=less_color, label='Reject null hypothesis for FE test alternative "less"')
    handles.append(patch)
    handles.append(patch2)
    plt.legend(handles=handles, loc='upper right')
    
    plt.title(u"Window=\u00B1{:,} bp; {}".format(window,title))
    plt.show()


def print_and_plot(stats_df, window, r1c1_label, r1c2_label, title, plot=True):  # added plot=True
    """
    (1) Apply Bonferroni correction
    (2) Plot p-values (histogram)
    (3) Plot barplot across subfamilies and indicating enrichment
    
    Originally:
    stats_df['bonferroni_reject'], stats_df['bonferroni_p'],_,_ = multipletests(stats_df['fisher_p_value'],
                                                                                method='bonferroni')
    """
    stats_df['bonferroni_reject_greater'], stats_df['bonferroni_p_greater'],_,_ = multipletests(stats_df['fisher_greater_p'], method='bonferroni')
    stats_df['bonferroni_reject_less'], stats_df['bonferroni_p_less'],_,_ = multipletests(stats_df['fisher_less_p'], method='bonferroni')
    #hist_p_vals(stats_df)
    if plot:
        bar_plot_subfamilies(stats_df, window, r1c1_label, r1c2_label, title)
    return stats_df


def test_table_across_windows_subfamily(window, flank_df, control_flank_df, ir_flank_df,
                                        control_ir_flank_df, window_df, alu=True):
    window_df_counts = window_df['alu_subfamily'].value_counts()
    inverted=True
    pair=False
    comb=False
    
    non_ir_flank_df = filter_non_inverted_pair(flank_df.copy())
    control_non_ir_flank_df = filter_non_inverted_pair(control_flank_df.copy())
    r2c1 = len(non_ir_flank_df)  # number of Alu pairs with at least 1 flanking non-IR Alu pair
    r2c2 = len(control_non_ir_flank_df)
    
    if alu:
        if (inverted == True) and (pair == False) and (comb == False):  # single member
            stats_dct_list1 = []
            stats_dct_list2 = []
            for ind_alu in list(set(window_df_counts.index)):
                ir_flank_df_single_alu = ir_flank_df[(ir_flank_df['upstream_alu_subfamily'] == ind_alu) | (ir_flank_df['downstream_alu_subfamily'] == ind_alu)]
                control_ir_flank_df_single_alu = control_ir_flank_df[(control_ir_flank_df['upstream_alu_subfamily'] == ind_alu) | (control_ir_flank_df['downstream_alu_subfamily'] == ind_alu)]
                
                non_ir_flank_df_single_alu = non_ir_flank_df[(non_ir_flank_df['upstream_alu_subfamily'] == ind_alu) | (non_ir_flank_df['downstream_alu_subfamily'] == ind_alu)]
                control_non_ir_flank_df_single_alu = control_non_ir_flank_df[(control_non_ir_flank_df['upstream_alu_subfamily'] == ind_alu) | (control_non_ir_flank_df['downstream_alu_subfamily'] == ind_alu)]
                            
                r1c1 = len(ir_flank_df_single_alu)
                r1c2 = len(control_ir_flank_df_single_alu)
                
                #r2c1 = len(ir_flank_df) - len(ir_flank_df_single_alu)  # previously 
                r2c1 = len(non_ir_flank_df_single_alu)
                #r2c2 = len(control_ir_flank_df) - len(control_ir_flank_df_single_alu)  # previously
                r2c2 = len(control_non_ir_flank_df_single_alu)
                
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted pair & 1+ {}".format(ind_alu),
                                                        #r2_label="Inverted pair & {} not in pair".format(ind_alu),
                                                        r2_label="Non-inverted pair & 1+ {}".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list1.append(result_dct)
            
                # r1c1 as above
                # r1c2 as above
                r2c1 = len(flank_df) - len(ir_flank_df_single_alu)
                r2c2 = len(control_flank_df) - len(control_ir_flank_df_single_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted pair & 1+ {}".format(ind_alu),
                                                        r2_label="(non-inverted) or ({} not in pair)".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list2.append(result_dct)
            
            stats_df1 = pd.DataFrame(stats_dct_list1)
            plot1_r1c1_label = "At least one member of inverted pair flanking skipped exons"
            plot1_r1c2_label = "At least one member of inverted pair flanking constitutive exons"
            plot1_title = "Among inverted pairs, at least one Alu subfamily member in pair vs not" 
            stats_df1 = print_and_plot(stats_df1, window, plot1_r1c1_label, plot1_r1c2_label, plot1_title)
            
            stats_df2 = pd.DataFrame(stats_dct_list2)
            plot2_r1c1_label = plot1_r1c1_label
            plot2_r1c2_label = plot1_r1c2_label
            plot2_title = "At least one Alu subfamily in inverted pair vs (non-inverted or not in pair)"
            stats_df2 = print_and_plot(stats_df2, window, plot2_r1c1_label, plot2_r1c2_label, plot2_title)
            
            stats_df1['case'] = 1
            stats_df2['case'] = 2
            return pd.concat([stats_df1, stats_df2], ignore_index=True)


def run_across_windows(hexevent_skipped_exon_df, 
                       exonskipdb_skipped_exon_df,
                       control_df,
                       exon_sets,
                       header="bedops",
                       exon_centric=False,
                       discrete=False):
    """
    Use filtered exon sets by running generate_dfs_windows_subset()
    Modified from: 02_ir_alus_plots_windows.ipynb
    Wrapper to run the main enrichment function across range of window sizes.
    """
    results = []
    plot_dct = {}
    for window in np.arange(500,5001,500):
        print(window)
        window_dct = (
            generate_dfs_windows_subset(
                window, 
                hexevent_skipped_exon_df,
                exonskipdb_skipped_exon_df, 
                control_df,
                exon_sets["skippable"], 
                exon_sets["constitutive"],
                ir=True,
                header=header, 
                discrete=discrete)
            )

        if not exon_centric:
            results1 = test_table_across_windows(flank_df=window_dct["skippable"],
                                                 control_flank_df=window_dct["constitutive"],
                                                 ir_flank_df=window_dct["skippable_inverted"],
                                                 control_ir_flank_df=window_dct["constitutive_inverted"])
        if exon_centric:
            results1 = test_table_across_windows_exon(flank_df=window_dct["skippable"],
                                                      control_flank_df=window_dct["constitutive"],
                                                      ir_flank_df=window_dct["skippable_inverted"],
                                                      control_ir_flank_df=window_dct["constitutive_inverted"])
        results1['window'] = window
        results.append(results1)
        plot_dct[window] = {"r1c1": results1["r1c1"],
                            "r1c2": results1["r1c2"],
                            "r2c1": results1["r2c1"],
                            "r2c2": results1["r2c2"],
                            "log10p": -np.log10(results1["fisher_greater_p"])}
        print(plot_dct[window])
        results1 = None
    
    #return {k: v for d in results for k, v in d.items()}
    return results, plot_dct


def run_subfamily_enrichment(df_skip_hex, df_skip_esdb, df_const, exon_sets, pair=False, comb=False):
    results = []
    plot_dct = {}
    for window in np.arange(500,5001,500):
        print(window)
        window_dct = generate_dfs_windows_subset(
            window, 
            hexevent_skipped_exon_df   = df_skip_hex, 
            exonskipdb_skipped_exon_df = df_skip_esdb,
            control_df                 = df_const, 
            skippable_set              = exon_sets["skippable"], 
            constitutive_set           = exon_sets["constitutive"],
            header                     = "bedops", 
            ir                         = True)
        # compare to generate_dfs_windows
        if not pair:
            params_dct = {"alu":True, "pair":False, "inverted":True, "comb": False}
        elif pair:
            if not comb:
                params_dct = {"alu":True, "pair":True, "inverted":True, "comb": False}
            elif comb:
                params_dct = {"alu":True, "pair":True, "inverted":False, "comb":True}
        print(params_dct)
        
        results_window = test_table(window,
                                    window_dct["skippable"],
                                    window_dct["constitutive"],
                                    window_dct["skippable_inverted"],
                                    window_dct["constitutive_inverted"],
                                    window_dct["window"],
                                    alu=params_dct["alu"],
                                    pair=params_dct["pair"],
                                    inverted=params_dct["inverted"],
                                    comb=params_dct["comb"])
        results_window['window'] = window
        results.append(results_window)
        
        plot_dct[window] = {"r1c1": results_window["r1c1"], 
                            "r1c2": results_window["r1c2"],
                            "r2c1": results_window["r2c1"],
                            "r2c2": results_window["r2c2"],
                            "log10p": -np.log10(results_window["fisher_greater_p"])}
        
        print(plot_dct[window])
    
    results_df_singlealuirpair_df = pd.concat(results, ignore_index=True)
    return results, plot_dct, results_df_singlealuirpair_df
    #swarmplot_p_values(results_df_singlealuirpair_df, title=None, log=False)


def run_across_windows_subfamilies(skipped_exon_df: pd.DataFrame,
                                   control_df: pd.DataFrame,
                                   scenario: str,
                                   header="bedops"):
    """
    Modified from: 02_ir_alus_plots_windows.ipynb
    Wrapper to run the main enrichment function across range of window sizes.
    
    Scenarios:
    - "single_alu"
    - "single_alu_IR_pair"
    - "two_alus_IR_pair"
    - "two_alus"
    """
    results = []
    plot_dct = {}
    for window in np.arange(500,5001,500):
        print(window)
        window_dct = generate_dfs_windows(window, skipped_exon_df, control_df, ir=True, header=header)
        
        if scenario == "single_alu":
            params_dct = {"alu":True, "pair":False, "inverted":False, "comb": False}
        elif scenario == "single_alu_IR_pair":
            params_dct = {"alu":True, "pair":False, "inverted":True, "comb": False}
        elif scenario == "two_alus_IR_pair":
            params_dct = {"alu":True, "pair":True, "inverted":True, "comb": False}
        elif scenario == "two_alus":
            params_dct = {"alu":True, "pair":True, "inverted":False, "comb":True}
        
        results_window = test_table(window,
                                    window_dct["skippable"],
                                    window_dct["constitutive"],
                                    window_dct["skippable_inverted"],
                                    window_dct["constitutive_inverted"],
                                    window_dct["window"],
                                    alu=params_dct["alu"],
                                    pair=params_dct["pair"],
                                    inverted=params_dct["inverted"],
                                    comb=params_dct["comb"])
        results_window['window'] = window
        results.append(results_window)
        
        plot_dct[window] = {"r1c1": results_window["r1c1"], 
                            "r1c2": results_window["r1c2"],
                            "r2c1": results_window["r2c1"],
                            "r2c2": results_window["r2c2"],
                            "log10p": -np.log10(results_window["fisher_greater_p"])}
        print(plot_dct[window])
    #return {k: v for d in results for k, v in d.items()}
    return results, plot_dct


def filter_inclusion_alu_pair(df: pd.DataFrame, alu_pair: list):
    """Filter based on Alu pair (downstream/upstream order agnostic)

    Args:
        df (pd.DataFrame): flank_df, ir_flank_df, etc 
        alu_pair (list): len=2 of Alu subfamily names

    Returns:
        _type_: pd.DataFrame
    """
    df_alu_pair = (df[((df["upstream_alu_subfamily"] == alu_pair[0]) & (df["downstream_alu_subfamily"] == alu_pair[1])) | 
                      ((df["upstream_alu_subfamily"] == alu_pair[1]) & (df["downstream_alu_subfamily"] == alu_pair[0]))])
    return df_alu_pair


def filter_inclusion_alu_single(df: pd.DataFrame, alu: str):
    """Filter based on single Alu (downstream/upstream agnostic)

    Args:
        df (pd.DataFrame): flank_df, ir_flank_df, etc 
        alu (str): Alu subfamily name

    Returns:
        _type_: pd.DataFrame
    
    """
    df_alu_single = df[(df["upstream_alu_subfamily"] == alu) | (df["downstream_alu_subfamily"] == alu)]
    return df_alu_single


def filter_exact_alu_duplicate_pair(df: pd.DataFrame, alu: str):
    """Filter based on single Alu repeated twice (duplicate pairings)

    Args:
        df (pd.DataFrame): flank_df, ir_flank_df, etc 
        alu (str): Alu subfamily name

    Returns:
        _type_: pd.DataFrame
    
    """
    df_alu_dup = df[(df["upstream_alu_subfamily"] == alu) & (df["downstream_alu_subfamily"] == alu)]
    return df_alu_dup


# This function was copied directly from a previous notebook. I've extracted the "if not alu" block
# to generate dfs/dcts used in enrichment calculations and plots.
def test_table(window, flank_df, control_flank_df,
               ir_flank_df, control_ir_flank_df, window_df,
               alu=False, pair=False, inverted=False, comb=False):
    """
    Test 2x2 contigency table
    Columns are always: flanking skipped exons vs flanking constitutive exons
    
    window: length upstream and downstream of exon to search within
    
    # regenerating the four df params below takes a while
    flank_df = from create_flanking_df(window, closest_exon_df)
    control_flank_df = create_flanking_df(window, closest_const_exon_df)
    ir_flank_df = filter_inverted_pair(flank_df)
    control_ir_flank_df = filter_inverted_pair(control_flank_df)
    
    -alu: name of Alu subfamily (default: None)
    -pair: set True if want to compare pairs of Alu subfamilies, otherwise (single member of pair) leave default: False
    -inverted: set True if Alu member/pair should be inverted
    """
    #window_df = closest_exon_df[closest_exon_df['dist'].between(-window, window)].dropna(how='any')  # dropna doesn't change anything once window is applied
    #window_df = window_df.rename_axis('idx').reset_index().astype({'dist': 'int'})
    window_df_counts = window_df['alu_subfamily'].value_counts()
    
    # if no flags are set (alu, pair, inverted): IR vs non-IR (all subfamilies together)
    if not alu:
        r1c1 = len(ir_flank_df) # number of IR Alus flanking skipped exons
        r1c2 = len(control_ir_flank_df) # number of IR Alus flanking constitutive exons
        r2c1 = len(flank_df) - len(ir_flank_df) # number of non-IR Alus flanking skipped exons
        r2c2 = len(control_flank_df) - len(control_ir_flank_df)
        result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                r1_label="Inverted pair", r2_label="Non-inverted pair",
                                                verbose=True)
        #pprint(result_dct)
        return result_dct
    
    # next create all entries in 2x2 table according to flags/params
    elif alu:
        non_ir_flank_df = filter_non_inverted_pair(flank_df.copy())
        control_non_ir_flank_df = filter_non_inverted_pair(control_flank_df.copy())
        
        if (inverted == False) and (pair == False) and (comb == False): #single member
            stats_dct_list = []
            for ind_alu in list(set(window_df_counts.index)):
                # Replace with function on 23-12-11
                flank_df_single_alu = filter_inclusion_alu_single(flank_df, ind_alu)
                control_flank_df_single_alu = filter_inclusion_alu_single(control_flank_df, ind_alu)
                #flank_df_single_alu = flank_df[(flank_df['upstream_alu_subfamily'] == ind_alu) | (flank_df['downstream_alu_subfamily'] == ind_alu)]
                #control_flank_df_single_alu = control_flank_df[(control_flank_df['upstream_alu_subfamily'] == ind_alu) | (control_flank_df['downstream_alu_subfamily'] == ind_alu)]

                r1c1 = len(flank_df_single_alu)
                r1c2 = len(control_flank_df_single_alu)
                r2c1 = len(flank_df) - len(flank_df_single_alu)
                r2c2 = len(control_flank_df) - len(control_flank_df_single_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="At least 1 {} in pair".format(ind_alu),
                                                        r2_label="{} not in pair".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list.append(result_dct)
            stats_df = pd.DataFrame(stats_dct_list)
            plot_r1c1_label = "At least one member of pair flanking skipped exons"
            plot_r1c2_label = "At least one member of pair flanking constitutive exons"
            plot_title = "Member of pair vs not"
            stats_df = print_and_plot(stats_df, window, plot_r1c1_label, plot_r1c2_label, plot_title)
            return stats_df
        
        elif (inverted == True) and (pair == False) and (comb == False): #single member
            stats_dct_list1 = []
            stats_dct_list2 = []
            stats_dct_list3 = []
            for ind_alu in list(set(window_df_counts.index)):
                ir_flank_df_single_alu = filter_inclusion_alu_single(ir_flank_df, ind_alu)
                control_ir_flank_df_single_alu = filter_inclusion_alu_single(control_ir_flank_df, ind_alu)
                
                non_ir_flank_df_single_alu = filter_inclusion_alu_single(non_ir_flank_df, ind_alu)
                control_non_ir_flank_df_single_alu = filter_inclusion_alu_single(control_non_ir_flank_df, ind_alu)
                
                r1c1 = len(ir_flank_df_single_alu)
                r1c2 = len(control_ir_flank_df_single_alu)
                r2c1 = len(non_ir_flank_df_single_alu)
                r2c2 = len(control_non_ir_flank_df_single_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted pair & 1+ {}".format(ind_alu),
                                                        #r2_label="Inverted pair & {} not in pair".format(ind_alu),  # previous
                                                        r2_label="Non-inverted pair & 1+ {}".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list1.append(result_dct)
            
                # r1c1 as above
                # r1c2 as above
                r2c1 = len(flank_df) - len(ir_flank_df_single_alu)
                r2c2 = len(control_flank_df) - len(control_ir_flank_df_single_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted pair & 1+ {}".format(ind_alu),
                                                        r2_label="(Non-inverted) or ({} not in pair)".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list2.append(result_dct)
                
                # Original case 1
                r2c1 = len(ir_flank_df) - len(ir_flank_df_single_alu)
                r2c2 = len(control_ir_flank_df) - len(control_ir_flank_df_single_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted pair & 1+ {}".format(ind_alu),
                                                        r2_label="Inverted pair & {} not in pair".format(ind_alu),
                                                        alu=ind_alu, chi=False)
                stats_dct_list3.append(result_dct)
            
            stats_df1 = pd.DataFrame(stats_dct_list1)
            plot1_r1c1_label = "Flanking skipped exons"
            plot1_r1c2_label = "Flanking constitutive exons"
            plot1_title = ("At least one copy of Alu subfamily in inverted pair vs at least " +
                           "one copy of Alu subfamily in non-inverted pair (Alu counts)")
            stats_df1 = print_and_plot(stats_df1, window, plot1_r1c1_label, plot1_r1c2_label, plot1_title)
            
            stats_df2 = pd.DataFrame(stats_dct_list2)
            plot2_r1c1_label = plot1_r1c1_label
            plot2_r1c2_label = plot1_r1c2_label
            plot2_title = "At least one Alu subfamily in inverted pair vs (non-inverted or not in pair)"
            stats_df2 = print_and_plot(stats_df2, window, plot2_r1c1_label, plot2_r1c2_label, plot2_title)
            
            stats_df3 = pd.DataFrame(stats_dct_list3)
            plot3_r1c1_label = plot1_r1c1_label
            plot3_r1c2_label = plot1_r1c2_label
            plot3_title = "Among inverted pairs, at least one Alu subfamily member in pair vs not" 
            stats_df3 = print_and_plot(stats_df3, window, plot3_r1c1_label, plot3_r1c2_label, plot3_title)
            
            stats_df1['case'] = 1
            stats_df2['case'] = 2
            stats_df3['case'] = 3
            return pd.concat([stats_df1, stats_df2, stats_df3], ignore_index=True)

        elif (inverted == True) and (pair == True) and (comb == False):  # both alus in pair
            stats_dct_list1 = []
            stats_dct_list2 = []
            stats_dct_list3 = []
            for ind_alu in list(set(window_df_counts.index)):
                flank_df_dup_alu = filter_exact_alu_duplicate_pair(flank_df, ind_alu)
                control_flank_df_dup_alu = filter_exact_alu_duplicate_pair(control_flank_df, ind_alu)

                r1c1 = len(flank_df_dup_alu)
                r1c2 = len(control_flank_df_dup_alu)
                r2c1 = len(flank_df) - len(flank_df_dup_alu)
                r2c2 = len(control_flank_df) - len(control_flank_df_dup_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="{} and {}".format(ind_alu, ind_alu),
                                                        r2_label="All other pairs",
                                                        alu=ind_alu, chi=False)
                stats_dct_list1.append(result_dct)

                ir_flank_df_dup_alu = filter_exact_alu_duplicate_pair(ir_flank_df, ind_alu)
                control_ir_flank_df_dup_alu = filter_exact_alu_duplicate_pair(control_ir_flank_df, ind_alu)
            
                r1c1 = len(ir_flank_df_dup_alu)
                r1c2 = len(control_ir_flank_df_dup_alu)
                r2c1 = len(ir_flank_df) - len(ir_flank_df_dup_alu)
                r2c2 = len(control_ir_flank_df) - len(control_ir_flank_df_dup_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted {} and {}".format(ind_alu, ind_alu),
                                                        r2_label="All other inverted pairs",
                                                        alu=ind_alu, chi=False)
                stats_dct_list2.append(result_dct)
            
                # r1c1 as above
                # r1c2 as above
                r2c1 = len(flank_df) - len(ir_flank_df_dup_alu)
                r2c2 = len(control_flank_df) - len(control_ir_flank_df_dup_alu)
                result_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted {} and {}".format(ind_alu, ind_alu),
                                                        r2_label="All other inverted and non-inverted pairs",
                                                        alu=ind_alu, chi=False)
                stats_dct_list3.append(result_dct)
            
            stats_df1 = pd.DataFrame(stats_dct_list1)
            plot1_r1c1_label = "Identical Alus form pair flanking skipped exons"
            plot1_r1c2_label = "Identical Alus form pair flanking constitutive exons"
            plot1_title = "Two Alus form pair vs all other pairs (inverted and non-inverted)"
            stats_df1 = print_and_plot(stats_df1, window, plot1_r1c1_label, plot1_r1c2_label, plot1_title)
            
            stats_df2 = pd.DataFrame(stats_dct_list2)
            plot2_r1c1_label = "Identical Alus form inverted pair flanking skipped exons"
            plot2_r1c2_label = "Identical Alus form inverted pair flanking constitutive exons"
            plot2_title = "Inverted identical Alus vs all other inverted pairs"
            stats_df2 = print_and_plot(stats_df2, window, plot2_r1c1_label, plot2_r1c2_label, plot2_title)

            stats_df3 = pd.DataFrame(stats_dct_list3)
            plot3_r1c1_label = plot2_r1c1_label
            plot3_r1c2_label = plot2_r1c2_label
            plot3_title = "Inverted identical Alus vs all other pairs (inverted and non-inverted)"
            stats_df3 = print_and_plot(stats_df3, window, plot3_r1c1_label, plot3_r1c2_label, plot3_title)
        
            stats_df1['case'] = 1
            stats_df2['case'] = 2
            stats_df3['case'] = 3
            return pd.concat([stats_df1, stats_df2, stats_df3], ignore_index=True)
        
        elif (pair == True) and (comb == True):
            result_dct_list = []
            result_dct_list2 = []
            result_dct_list3 = []
            for alu_pair in list(combinations(list(set(window_df_counts.index)),2)):
                
                # Replaced long chained statements below with function on 23-12-08
                flank_df_alu_pair = filter_inclusion_alu_pair(flank_df, alu_pair)
                control_flank_df_alu_pair = filter_inclusion_alu_pair(control_flank_df, alu_pair)
                ir_flank_df_alu_pair = filter_inclusion_alu_pair(ir_flank_df, alu_pair)
                control_ir_flank_df_alu_pair = filter_inclusion_alu_pair(control_ir_flank_df, alu_pair)
                
                non_ir_flank_df_alu_pair = filter_inclusion_alu_pair(non_ir_flank_df, alu_pair)
                control_non_ir_flank_df_alu_pair = filter_inclusion_alu_pair(control_non_ir_flank_df, alu_pair)
                                
                # Find total num of pairs (upstream/downstream agnostic)
                # Find num of inverted pairs (upstream/downstream agnostic)
                result_dct_list.append({"pair": alu_pair, "flank": len(flank_df_alu_pair), "control_flank": len(control_flank_df_alu_pair), "ir_flank": len(ir_flank_df_alu_pair), "ir_control_flank": len(control_ir_flank_df_alu_pair)})       
                
                r1c1 = len(ir_flank_df_alu_pair)
                r1c2 = len(control_ir_flank_df_alu_pair)
                r2c1 = len(ir_flank_df) - len(ir_flank_df_alu_pair)
                r2c2 = len(control_ir_flank_df) - len(control_ir_flank_df_alu_pair)
                results_dct = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                        r1_label="Inverted {}".format(alu_pair),
                                                        r2_label="All other inverted pairs",
                                                        alu=alu_pair, chi=False, verbose=False)
                result_dct_list2.append(results_dct)
                
                r2c1 = len(non_ir_flank_df_alu_pair)
                r2c2 = len(control_non_ir_flank_df_alu_pair)
                results_dct3 = generate_contingency_table(r1c1, r1c2, r2c1, r2c2,
                                                          r1_label="Inverted {}".format(alu_pair),
                                                          r2_label="Non-inverted {}".format(alu_pair),
                                                          alu=alu_pair, chi=False, verbose=False)
                result_dct_list3.append(results_dct3)
                
            results_df = pd.DataFrame(result_dct_list)
            results_df2 = pd.DataFrame(result_dct_list2)
            results_df3 = pd.DataFrame(result_dct_list3)
            
            #stats_df1 = pd.DataFrame(stats_dct_list1)
            plot1_r1c1_label = "Flanking skipped exons"
            plot1_r1c2_label = "Flanking constitutive exons"
            plot1_title = ("Alu pair in inverted pair vs " +
                           "Alu pair in non-inverted pair (Alu counts)")
            results_df3 = print_and_plot(results_df3, window, plot1_r1c1_label, plot1_r1c2_label, plot1_title, plot=False)
            
            plot_alu_pairs_combo(results_df, window, title="Alu pairs")
            return results_df3  # previously results_df2
    
    bonferroni_significant_p = 0.05/51  #usually alpha/n_tests == 0.05/51 but repeated twice (less and greater)
    print("Bonferroni significant p-value for one-tailed tests across Alu subfamilies: {:0.2e}".format(bonferroni_significant_p))
    bonferroni_significant_p_two = 0.05/(51*2)  #usually alpha/n_tests == 0.05/51 but repeated twice (less and greater)
    print("Bonferroni significant p-value for 2*one-tailed tests across Alu subfamilies: {:0.2e}".format(bonferroni_significant_p_two))


def plot_alu_pairs_combo(results_df: pd.DataFrame, window: str, title: str):
    sns.set_theme()
    sns.set(font_scale=0.9)
    fig, ax = plt.subplots()
    ax.set_xlabel("Alu subfamily pair combinations")
    ax.set_ylabel("Frequency")
    
    results_df = results_df.copy().sort_values(by=['flank'], ascending=False)  # added copy
    results_idx = list(results_df.set_index('pair').index)
            
    ax1 = (results_df.head(100)
           .set_index('pair')['control_flank']
           .plot(kind='bar', log=False, figsize=(16,5), color='grey', label="constitutive flank"))
    
    ax2 = (results_df.head(100)
           .set_index('pair')['flank']
           .plot(kind='bar', log=False, figsize=(16,5), color='steelblue', label="skipped flank"))
    
    ax3 = (results_df.head(100)
           .set_index('pair')['ir_control_flank']
           .plot(kind='bar', log=False, figsize=(16,5), color='purple', label="IR constitutive flank"))
    
    ax4 = (results_df.head(100)
           .set_index('pair')['ir_flank']
           .plot(kind='bar', log=False, figsize=(16,5), color='pink', label="IR skipped flank"))
    
    #ax2.get_xaxis().set_visible(False)
    plt.legend(loc='upper right')
    title="alu pairs"
    plt.title(u"Window=\u00B1{:,} bp; {}".format(window,title))
    plt.show()


def generate_melt_df_for_enrichment_plot(enrichment_dct):
    #r1c1 = len(ir_flank_df)
    #r1c2 = len(control_ir_flank_df)
    #r2c1 = len(flank_df) - len(ir_flank_df)
    #r2c2 = len(control_flank_df) - len(control_ir_flank_df)
    
    # df_melt = pd.melt(pd.DataFrame.from_dict(enrichment_dct, orient="index")
    #                   .reset_index()
    #                   .rename(columns={"index": "window",
    #                                    "r1c1": "Inverted flanking skippable",
    #                                    "r1c2": "Inverted flanking constitutive",
    #                                    "r2c1": "Non-inverted flanking skippable",
    #                                    "r2c2": "Non-inverted flanking constitutive"}),
    #                   id_vars=["window", "log10p"], var_name="group")
    df_premelt = (pd.DataFrame
                  .from_dict(enrichment_dct, orient="index")
                  .reset_index()
                  .rename(columns={"index": "window",
                                   "r1c1": "Inverted flanking skippable",
                                   "r1c2": "Inverted flanking constitutive",
                                   "r2c1": "Non-inverted flanking skippable",
                                   "r2c2": "Non-inverted flanking constitutive"}))
    df_premelt["p-value"] = 10**(-df_premelt["log10p"])
    df_premelt["bonferroni_reject"], df_premelt["bonferroni_p"],_,_ = multipletests(df_premelt["p-value"], method='bonferroni')
    df_premelt["bonferroni_neglog10_p"] = -np.log10(df_premelt["bonferroni_p"])
    df_melt = pd.melt(df_premelt, 
                      id_vars=["window", "log10p", "p-value", "bonferroni_reject", "bonferroni_p", "bonferroni_neglog10_p"], 
                      var_name="group")
    return df_melt


def generate_enrichment_plot(df_melt, ylabel, xlabel=False, legend=True, ax_legend=None, legend_pos=None, ax1=None):
    """
    Plot a line plot (p-value) on top of a staggered bar plot (number of Alu pairs).
    
    Visualization of FE p-values based off of 2x2 contingency tables.
    
    Re-uses functions and results from last summer's FE analyses.
    """

    hue_order=["Inverted flanking skippable", "Non-inverted flanking skippable",
               "Inverted flanking constitutive", "Non-inverted flanking constitutive"]

    blue="#0077BB"
    light_blue="#86B4CF"
    orange="#EE7733"
    light_orange="#EEC5AE"

    customPalette = [blue, light_blue, orange, light_orange]
    # previously: reversed(sns.color_palette("Paired")[6:10])

    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2
    
    fontsize_min = 5

    sns.set_theme() #font_scale=1.4
    sns.set_style('ticks')
    plt.rcParams['font.family'] = 'DejaVu Sans'

    if not ax1:
        fig, ax1 = plt.subplots(figsize=(14,8))

    sns.barplot(data=df_melt, x="window", y="value", hue="group", hue_order=hue_order,
                alpha=1, palette=customPalette, ax=ax1, width=0.9, edgecolor='none') #edgecolor="black"
    
    ax1.set_xlim([-0.5, len(df_melt['window'].unique()) - 0.5])
    
    ax2 = ax1.twinx()
    ax2.plot(ax1.get_xticks(), df_melt["bonferroni_neglog10_p"].iloc[0:10], color="g", linewidth=0.25)
    ax1.xaxis.grid(True, linewidth=tick_params_width) 

    if not xlabel:
        ax1.set_xlabel("") # Window
    elif xlabel:
        ax1.set_xlabel("Flanking window (bp)", fontsize=fontsize_min, labelpad=axis_label_pad)
    
    ax1.set_ylabel(f"{ylabel}", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax2.set_ylabel('-log10(p_adj)', color='g', fontsize=fontsize_min, labelpad=axis_label_pad)
    handles, labels = ax1.get_legend_handles_labels()

    ##ax1.legend(handles=handles[0:], labels=labels[0:], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # ref: https://stackoverflow.com/questions/4700614
    #ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.30), #-0.25
    #      ncol=2, fancybox=True, shadow=False, borderaxespad=0.)
    if legend:
        if legend_pos == "right":
            ax1.legend(loc='lower center', bbox_to_anchor=(0.8, 0.77), #-0.25
                       ncol=1, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
        else:
            if ax_legend:
                #ax_legend.legend(loc='lower center', bbox_to_anchor=(0.5, -0.1), #bbox_to_anchor=(0.27, 0.55)
                #       ncol=2, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
                ax_legend.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.1),
                                 ncol=2, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
                ax_legend.axis('off')
                ax1.get_legend().remove()
            else:
                ax1.legend(loc='lower center', bbox_to_anchor=(1, 0.55), #bbox_to_anchor=(0.27, 0.55)
                       ncol=2, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
            #ax1.legend(loc='lower center', bbox_to_anchor=(0.27, 0.55), #(0.23, 0.62) for font_scale=1.2 #-0.25
            #           ncol=1, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
    else:
        ax1.get_legend().remove()
    
    # ref: https://stackoverflow.com/questions/52392855
    for tick_label in ax2.axes.get_yticklabels():
        tick_label.set_color("g")
    
    # ref: https://stackoverflow.com/questions/39376888
    ax1.get_yaxis().set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))
    
    plt.ylim(0, None)

    ax1.tick_params(axis='x', labelsize=fontsize_min, 
                    width=tick_params_width, length=tick_params_length, pad=tick_params_pad) #pad=-2
    ax1.tick_params(axis='y', labelsize=fontsize_min, 
                    width=tick_params_width, length=tick_params_length, pad=tick_params_pad) #pad=0
    
    for label in ax1.get_xticklabels():
        label.set_rotation(45)
    
    ax2.axhline(y=-np.log10(0.05), color='g', linestyle='--', linewidth=0.5)
    ax2.tick_params(axis='y', labelsize=fontsize_min, 
                    width=tick_params_width, length=tick_params_length, pad=tick_params_pad) #pad=0.1
    
    for spine in ax1.spines.values():
        spine.set_linewidth(0.25)
    for spine in ax2.spines.values():
        spine.set_linewidth(0.25) 

    #plt.title(f"{title}")
    #plt.show()


def swarmplot_p_values(results_df, title=None, log=False, use_string=False, 
                       legend=True, legend_pos="bottom", ylabel=True, ax=None, ax_legend=None):
    """ Plot (log scale) p-values where bonferroni_reject_greater is True 
    (reached statistical significance of one-tailed FE test after Bonferonni correction).
    """
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2
    
    fontsize_min = 5

    sns.set_theme() #font_scale=1.2
    sns.set_style('ticks')
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    colorblind_palette = sns.color_palette('colorblind')
    #bidirectional_palette = sns.color_palette('Spectral') #Spectral,coolwarm
    
    cmap = plt.get_cmap('Spectral')
    n_colors = 10  # Number of discrete colors
    colors = [cmap(i) for i in np.linspace(0, 1, n_colors)]
    bidirectional_palette = sns.color_palette(colors)
    # selected_colors = colorblind_palette[2:8]
    # selected_colors = selected_colors[1:] + [selected_colors[0]]
    
    #bidirectional_palette = sns.color_palette('coolwarm', as_cmap=True)
    #n_colors = 6
    #selected_colors = sns.color_palette(bidirectional_palette, n_colors=n_colors)
    
    if not ax:
        fig, ax = plt.subplots(figsize=(8,4))
    
    results_df = results_df.copy()
    results_df['window'] = results_df['window'].astype('category')
    results_df["bonferroni_p_greater_neglog10"] = -np.log10(results_df["bonferroni_p_greater"])
    
    if not use_string:
        x_col = "alu_subfamily"
    elif use_string:
        x_col = "alu_subfamily_str"
    sns.swarmplot(x=x_col, y="bonferroni_p_greater_neglog10", hue="window", s=2, dodge=False,
                  data=results_df[results_df['bonferroni_reject_greater']], 
                  palette=bidirectional_palette, ax=ax) #palette=colorblind_palette; y="bonferroni_p_greater" #s=7
    #ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    
    ax.xaxis.grid(True)
    ax.set_ylim(0, 12)
    ax.axhline(y=-np.log10(0.05), color='steelblue', linestyle='--', linewidth=0.5)
    
    if ylabel:
        if log:
            #ax.set_yscale('log')
            #ax.set_ylabel("Bonferroni-corrected p-value \n(alternative: greater)") # log(bonferroni_p_greater)
            ax.set_ylabel("-log10(p_adj)", fontsize=fontsize_min) # log(bonferroni_p_greater)
        else:
            ax.set_ylabel("Bonferroni-corrected p-value \n(alternative: greater)", fontsize=fontsize_min)
    else:
        ax.set_ylabel(None)
    
    #ax.set_xlabel("Alu subfamily")
    ax.set_xlabel(None)
    if title:
        plt.title("{}".format(title))
    
    handles, labels = ax.get_legend_handles_labels()
    if legend:
        if legend_pos == "bottom":
            # ref: https://stackoverflow.com/questions/4700614
            ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.52), #(0.5, -0.4)
                      ncol=5, fancybox=True, shadow=False, borderaxespad=0., 
                      fontsize=fontsize_min, title="Window (bp)")
        elif legend_pos == "right":
            # ax.legend(loc='lower center', bbox_to_anchor=(1.24, -0.1), # (1.21, 0.1),(1.13, 0.1)
            #           ncol=1, fancybox=True, shadow=False, borderaxespad=0., 
            #           fontsize=fontsize_min, title="Window (bp)", title_fontsize=fontsize_min)
            ax_legend.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.5),
                             ncol=5, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min, #ncol=10
                             title="Window (bp)", title_fontsize=fontsize_min, 
                             columnspacing=0.5, borderpad=0.5, handletextpad=0.1, handlelength=0.6)
            ax_legend.axis('off')
            ax.get_legend().remove()
    else:
        ax.get_legend().remove()
    
    ax.tick_params(axis='x', labelsize=fontsize_min, 
                   width=tick_params_width, length=tick_params_length, pad=tick_params_pad)
    for label in ax.get_xticklabels():
        label.set_rotation(45)  #25
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   width=tick_params_width, length=tick_params_length, pad=tick_params_pad)

    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    #plt.show()


def plot_or(df, ax=None):
    #sns.set_style("whitegrid")
    fontsize_min = 5

    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    sns.set_theme() #font_scale=1.4
    sns.set_style('ticks')
    plt.rcParams['font.family'] = 'DejaVu Sans'

    if not ax:
        fig, ax = plt.subplots(figsize=(7, 3))

    ax.set_xticks(range(0, int(max(df["window"])) + 1, 500))
    ax.set_xlim(50, max(df["window"]) + 450)
    ax.set_ylim(1, 1.5)
    
    ax.xaxis.grid(True, linewidth=0.25) 
    ax.errorbar(df["window"], df["odds_ratio"], 
                 yerr=(df["odds_ratio"] - df["ci_low"],
                       df["ci_high"] - df["odds_ratio"]), 
                 fmt='o', markersize=0.5, color="black", capsize=0.5, elinewidth=0.5) #markersize=5
    
    ax.tick_params(axis='x', labelsize=fontsize_min, 
                   width=tick_params_width, length=tick_params_length, pad=tick_params_pad)
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   width=tick_params_width, length=tick_params_length, pad=tick_params_pad)

    for label in ax.get_xticklabels():
        label.set_rotation(45)
    
    ax.set_xlabel("Flanking window (bp)", fontsize=fontsize_min, labelpad=axis_label_pad) #labelpad=20
    ax.set_ylabel("Odds ratio", fontsize=fontsize_min, labelpad=axis_label_pad)
    
    plt.title("")
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
