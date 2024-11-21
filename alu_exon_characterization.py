#!/usr/bin/env python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests
import numpy as np
import zepid
import matplotlib as mpl


def add_label_col(df, label_col):
    df = df.copy()
    df["label"] = label_col
    return df


def add_alu_len_exon_len_cols(df, pairs):
    df["exon_len"] = df["exon_end"] - df["exon_start"]
    df["exon_len"] = df["exon_len"].abs()
    if pairs:
        df["upstream_alu_len"]  = df["upstream_alu_end"] - df["upstream_alu_start"]
        df["upstream_alu_len"]  = df["upstream_alu_len"].abs()
        df["downstream_alu_len"]  = df["downstream_alu_end"] - df["downstream_alu_start"]
        df["downstream_alu_len"]  = df["downstream_alu_len"].abs()
        df["dist1"] = df["downstream_alu_start"] - df["upstream_alu_end"]
        df["dist2"] = df["upstream_alu_start"] - df["downstream_alu_end"]
        df["intra_alu_len"] = df[["dist1", "dist2"]].abs().min(axis=1)
        return df[["upstream_alu_len","downstream_alu_len","exon_len","intra_alu_len","dist1","dist2"]]
    elif not pairs:
        df["alu_len"] = df["alu_end"] - df["alu_start"]
        df["alu_len"] = df["alu_len"].abs()
        return df[["alu_len","exon_len"]]


def generate_pairgrid_plot(df, cols, nrows_sample, window_size, title, inversion):
    #pairgrid_pair_cols = ["upstream_alu_len", "downstream_alu_len", "exon_len", 
    #                      "upstream_dist", "downstream_dist", "label"]
    pairgrid_pair_cols = cols
    
    column_labels = {
        'upstream_alu_len': "Upstream Alu\nLength",
        'downstream_alu_len': "Downstream Alu\nLength",
        'exon_len': "Exon Length",
        'upstream_dist': "Upstream Alu\nDistance",
        'downstream_dist': "Downstream Alu\nDistance",
        'label': "label"
    }

    label_col = {
        "skippable": "Skippable",
        "constitutive": "Constitutive",
        "skippable_inverted": "Skippable inverted",
        "constitutive_inverted": "Constitutive inverted",
        "skippable_noninverted": "Skippable non-inverted",
        "constitutive_noninverted": "Constitutive non-inverted"
    }
    
    df_renamed = df.rename(columns=column_labels)
    pairgrid_pair_cols_renamed = [column_labels[col] for col in pairgrid_pair_cols]
    df_renamed['label'] = df_renamed['label'].map(label_col)

    blue="#0077BB"
    light_blue="#86B4CF"
    orange="#EE7733"
    light_orange="#EEC5AE"
    if not inversion:
        hue_order = ["Skippable", "Constitutive"]
    elif inversion == "inverted":
        hue_order = ["Skippable inverted", "Constitutive inverted"]
    elif inversion == "non-inverted":
        hue_order = ["Skippable non-inverted", "Constitutive non-inverted"]

    g = sns.PairGrid(
        df_renamed[pairgrid_pair_cols_renamed].sample(n=10000, random_state=42), 
        hue="label", 
        hue_order=hue_order[::-1], # reversed
        palette=[orange, blue]
    )
    g.map_diag(sns.histplot, rasterized=True) # XXX
    g.map_upper(sns.scatterplot, s=10, alpha=.2, rasterized=True)
    g.map_lower(sns.kdeplot)

    #for ax in g.axes.flatten():
    #    if hasattr(ax, 'collections'):
    #        for col in ax.collections:
    #            col.set_rasterized(True)
    
    # def plot_diag(x, **kws):
    #     label = kws.pop('label', None)
    #     sns.histplot(x, **kws, zorder=1 if label == hue_order[0] else 2)

    # def plot_upper(x, y, **kws):
    #     label = kws.pop('label', None)
    #     sns.scatterplot(x=x, y=y, **kws, s=10, alpha=.2, zorder=1 if label == hue_order[0] else 2)

    # def plot_lower(x, y, **kws):
    #     label = kws.pop('label', None)
    #     sns.kdeplot(x=x, y=y, **kws, zorder=1 if label == hue_order[0] else 2)


    # g.map_diag(plot_diag)
    # g.map_upper(plot_upper)
    # g.map_lower(plot_lower)
    
    g.fig.set_size_inches(10, 10)  # Set the desired width and height in inches

    #g.map_offdiag(sns.scatterplot) 
    #g.map_diag(sns.kdeplot)
 
    g.add_legend(title='')
    g._legend.set_bbox_to_anchor((0.5, 0.98))
    g._legend.set_loc('center')
    g._legend.set_title(None)
    g._legend._legend_box.align = "center"
    
    for text in g._legend.get_texts():
        text.set_fontsize(14)
    g._legend._ncol = 2
    
    for ax in g.axes.flatten():
        ax.set_xlabel(ax.get_xlabel(), fontsize=14)
        ax.set_ylabel(ax.get_ylabel(), fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
    g.fig.suptitle(f"{title} ({window_size} windows; {nrows_sample:,.0f} row sampling with seed=42)", y=1.05, fontsize=16) #1.02
    plt.tight_layout()


def generate_singlealu_pairgrid_plot(df, cols, nrows_sample, window_size, title):
    #pairgrid_pair_cols = ["upstream_alu_len", "downstream_alu_len", "exon_len", 
    #                      "upstream_dist", "downstream_dist", "label"]
    pairgrid_pair_cols = cols
    
    # map labels
    column_labels = {
        'alu_len': "Alu Length",
        'exon_len': "Exon Length",
        'dist': "Distance",
        'label': "label"
    }

    label_col = {
        "skippable (hexevent)": "Skippable (HEXEvent)",
        "skippable (exonskipdb)": "Skippable (ExonSkipDB)",
        "constitutive (hexevent)": "Constitutive (HEXEvent)",
    }
    
    blue="#0077BB"
    light_blue="#86B4CF"
    teal = "#00bb99"
    royal_blue="#331fcc"
    orange="#EE7733"
    light_orange="#EEC5AE"

    hue_order = ["Skippable (ExonSkipDB)", 
                 "Constitutive (HEXEvent)", 
                 "Skippable (HEXEvent)"]
    color_order = [royal_blue, orange, teal]
        
    df_renamed = df.rename(columns=column_labels)
    pairgrid_pair_cols_renamed = [column_labels[col] for col in pairgrid_pair_cols]
    df_renamed['label'] = df_renamed['label'].map(label_col)
    #display(df_renamed)
    
    g = sns.PairGrid(
        df_renamed[pairgrid_pair_cols_renamed].sample(n=10000, random_state=42), 
        hue="label", 
        hue_order=hue_order, # reversed [::-1]
        palette=color_order,
        height=2
    )
    g.map_diag(sns.histplot, rasterized=True)
    g.map_upper(sns.scatterplot, s=10, alpha=.2, rasterized=True)
    g.map_lower(sns.kdeplot)
    g.fig.set_size_inches(10, 10)
    
    g.add_legend(title='')
    g._legend.set_bbox_to_anchor((0.4, 0.96)) #0.98
    g._legend.set_loc('center')
    g._legend.set_title(None)
    g._legend._legend_box.align = "center"
    
    for text in g._legend.get_texts():
        text.set_fontsize(14)
    g._legend._ncol = 3
    
    for ax in g.axes.flatten():
        ax.set_xlabel(ax.get_xlabel(), fontsize=14)
        ax.set_ylabel(ax.get_ylabel(), fontsize=14)
        ax.tick_params(axis='both', which='major', labelsize=12)
    g.fig.suptitle(f"{title} ({window_size} windows; {nrows_sample:,.0f} row sampling with seed=42)", y=1.1, fontsize=16) #1.02


def jitter(values, j):
    """
    Ref: https://stackoverflow.com/questions/64553046
    """
    return values + np.random.normal(j,0.1,values.shape)

def run_ks_2sample_test_2sided(df1, df2, col):
    return ks_2samp(df1[col], df2[col], alternative="two-sided", method="auto")

def run_ks_test(df1, df2_inverted, df3_noninverted, cv=0, n_sample=0):
    """
    cv=10
    n_sample=1000
    """
    pairgrid_pair_cols = ["upstream_alu_len", "downstream_alu_len", "exon_len", 
                          "upstream_dist", "downstream_dist", "intra_alu_len"]
    fixed_categories = ["skippable", "constitutive", "inverted", "noninverted", "all"]
    
    if cv > 0:
        ks_dct = {}
        for col in pairgrid_pair_cols:
            ks_dct[col] = {}
            for category in fixed_categories:
                ks_dct[col][category] = {"pval": [], "stat": []}
        
        for i in range(cv):
            df1_skippable = df1[df1["label"] == "skippable"].sample(n=n_sample, random_state=i)
            df1_constitutive = df1[df1["label"] == "constitutive"].sample(n=n_sample, random_state=i)
            df2_inverted_skippable = df2_inverted[df2_inverted["label"] == "skippable_inverted"].sample(n=n_sample, random_state=i)
            df2_inverted_constitutive = df2_inverted[df2_inverted["label"] == "constitutive_inverted"].sample(n=n_sample, random_state=i)
            df3_noninverted_skippable = df3_noninverted[df3_noninverted["label"] == "skippable_noninverted"].sample(n=n_sample, random_state=i) # new
            df3_noninverted_constitutive = df3_noninverted[df3_noninverted["label"] == "constitutive_noninverted"].sample(n=n_sample, random_state=i) # new
    
            for col in pairgrid_pair_cols:
                ks_res_skip = (
                    run_ks_2sample_test_2sided(df1=df1_skippable,
                                               df2=df2_inverted_skippable,
                                               col=col)
                )
                ks_res_cons = (
                    run_ks_2sample_test_2sided(df1=df1_constitutive,
                                               df2=df2_inverted_constitutive,
                                               col=col)
                )
                ks_res_inv = (
                    run_ks_2sample_test_2sided(df1=df2_inverted_skippable,
                                               df2=df2_inverted_constitutive,
                                               col=col)
                )
                ks_res_noninv = (
                    run_ks_2sample_test_2sided(df1=df3_noninverted_skippable,
                                               df2=df3_noninverted_constitutive,
                                               col=col)
                )
                ks_res_all = (
                    run_ks_2sample_test_2sided(df1=df1_skippable,
                                               df2=df1_constitutive,
                                               col=col)
                )
                ks_dct[col]["skippable"]["pval"].append(ks_res_skip.pvalue)
                ks_dct[col]["skippable"]["stat"].append(ks_res_skip.statistic)
                ks_dct[col]["constitutive"]["pval"].append(ks_res_cons.pvalue)
                ks_dct[col]["constitutive"]["stat"].append(ks_res_cons.statistic)
                ks_dct[col]["inverted"]["pval"].append(ks_res_inv.pvalue)
                ks_dct[col]["inverted"]["stat"].append(ks_res_inv.statistic)
                ks_dct[col]["noninverted"]["pval"].append(ks_res_noninv.pvalue)
                ks_dct[col]["noninverted"]["stat"].append(ks_res_noninv.statistic)
                ks_dct[col]["all"]["pval"].append(ks_res_all.pvalue)
                ks_dct[col]["all"]["stat"].append(ks_res_all.statistic)
                # ks_res_skip = ks_2samp(df1_skippable[col],
                #                  df2_inverted_skippable[col],
                #                  alternative="two-sided", method="auto")

                # ks_res_cons = ks_2samp(df1_constitutive[col],
                #                         df2_inverted_constitutive[col],
                #                         alternative="two-sided", method="auto")

                # ks_res_inv = ks_2samp(df2_inverted_skippable[col],
                #                     df2_inverted_constitutive[col],
                #                     alternative="two-sided", method="auto")

                # ks_res_all = ks_2samp(df1_skippable[col],
                #                     df1_constitutive[col],
                #                     alternative="two-sided", method="auto")
        return ks_dct

    df1_skippable = df1[df1["label"] == "skippable"] #.sample(n=1000, random_state=1)
    df1_constitutive = df1[df1["label"] == "constitutive"] #.sample(n=1000, random_state=1)
    df2_inverted_skippable = df2_inverted[df2_inverted["label"] == "skippable_inverted"] #.sample(n=1000, random_state=1)
    df2_inverted_constitutive = df2_inverted[df2_inverted["label"] == "constitutive_inverted"] #.sample(n=1000, random_state=1)
    df3_noninverted_skippable = df3_noninverted[df3_noninverted["label"] == "skippable_noninverted"] #.sample(n=1000, random_state=1)
    df3_noninverted_constitutive = df3_noninverted[df3_noninverted["label"] == "constitutive_noninverted"] #.sample(n=1000, random_state=1)

    print(f"len skippable: {len(df1_skippable)}")
    print(f"len constitutive: {len(df1_constitutive)}")
    print(f"len skippable_inverted: {len(df2_inverted_skippable)}")
    print(f"len constitutive_inverted: {len(df2_inverted_constitutive)}")
    print(f"len skippable_noninverted: {len(df3_noninverted_skippable)}")
    print(f"len constitutive_noninverted: {len(df3_noninverted_constitutive)}")
    
    ks_dct = {}
    for col in pairgrid_pair_cols:
        print(f"\n{col}: skippable vs inverted_skippable:")
        ks_res_skip = ks_2samp(df1_skippable[col],
                                 df2_inverted_skippable[col],
                                 alternative="two-sided", method="auto")
        print("skip all:", df1_skippable[col].describe())
        print("skip inv:", df2_inverted_skippable[col].describe())
        print(ks_res_skip)
        print(f"{col}: constitutive vs inverted_constitutive:")
        ks_res_cons = ks_2samp(df1_constitutive[col],
                                 df2_inverted_constitutive[col],
                                 alternative="two-sided", method="auto")
        print("const all:", df1_constitutive[col].describe())
        print("const inv:", df2_inverted_constitutive[col].describe())
        print(ks_res_cons)
        print(f"{col}: inverted_skippable vs inverted_constitutive:")
        ks_res_inv = ks_2samp(df2_inverted_skippable[col],
                              df2_inverted_constitutive[col],
                              alternative="two-sided", method="auto")
        print("inv skip:", df2_inverted_skippable[col].describe())
        print("inv const:", df2_inverted_constitutive[col].describe())
        print(ks_res_inv)
        print(f"{col}: noninverted_skippable vs noninverted_constitutive:")
        ks_res_noninv = ks_2samp(df3_noninverted_skippable[col],
                              df3_noninverted_constitutive[col],
                              alternative="two-sided", method="auto")
        print("noninv skip:", df3_noninverted_skippable[col].describe())
        print("noninv const:", df3_noninverted_constitutive[col].describe())
        print(ks_res_noninv)
        print(f"{col}: skippable vs constitutive:")
        ks_res_all = ks_2samp(df1_skippable[col],
                              df1_constitutive[col],
                              alternative="two-sided", method="auto")
        print("all skip:", df1_skippable[col].describe())
        print("all const:", df1_constitutive[col].describe())
        print(ks_res_all)
        
        ks_dct[col] = {"skippable": {"pval": ks_res_skip.pvalue, 
                                     "stat": ks_res_skip.statistic, 
                                     "stat_sign": ks_res_skip.statistic_sign}, 
                       "constitutive": {"pval": ks_res_cons.pvalue, 
                                        "stat": ks_res_cons.statistic, 
                                        "stat_sign": ks_res_cons.statistic_sign}, 
                       "inverted": {"pval": ks_res_inv.pvalue, 
                                    "stat": ks_res_inv.statistic, 
                                    "stat_sign": ks_res_inv.statistic_sign}, 
                       "noninverted": {"pval": ks_res_noninv.pvalue, 
                                    "stat": ks_res_noninv.statistic, 
                                    "stat_sign": ks_res_noninv.statistic_sign},
                       "all": {"pval": ks_res_all.pvalue, 
                               "stat": ks_res_all.statistic, 
                               "stat_sign": ks_res_all.statistic_sign}}
        print(ks_dct[col])

        if col == "exon_len":  #downstream_dist
            print(list(df1_skippable[col])[0:10])
            print(list(df2_inverted_skippable[col])[0:10])
            #df_skippable = pd.DataFrame.from_dict({"downstream_distance": sorted(list(df1_skippable[col])),
            #                                       "downstream_distance_inverted": sorted(list(df2_inverted_skippable[col]))},
            #                                      orient="index")
            
            df_skippable = pd.DataFrame.from_dict({"exon_len": sorted(list(df1_skippable[col])),
                                                   "exon_len_inv": sorted(list(df2_inverted_skippable[col]))},
                                                  orient="index")
            df_skippable = df_skippable.transpose()
            df_skippable["type"] = "skippable"
            df_constitutive = pd.DataFrame.from_dict({"exon_len": sorted(list(df1_constitutive[col])),
                                                      "exon_len_inv": sorted(list(df2_inverted_constitutive[col]))}, 
                                                     orient="index")
            df_constitutive = df_constitutive.transpose()
            df_constitutive["type"] = "constitutive"
            
            df_merged = pd.concat([df_skippable, df_constitutive], ignore_index=True)

            plt.subplots_adjust(top=0.9)
            sns.ecdfplot(df1_constitutive[col], stat="proportion", color="blue", label="all")
            sns.ecdfplot(df2_inverted_constitutive[col], stat="proportion", color="grey", label="inverted")
            plt.suptitle('Constitutive CDF', fontsize = 16)
            plt.legend()
            plt.show()
            plt.clf()
            
            plt.subplots_adjust(top=0.9)
            sns.ecdfplot(df1_skippable[col], stat="proportion", color="grey", label="all")
            sns.ecdfplot(df2_inverted_skippable[col], stat="proportion", color="purple", label="inverted")
            plt.suptitle('Skippable CDF', fontsize = 16)
            plt.legend()
            plt.show()

    return ks_dct


def plot_exon_length(alu_dct, window, exon_distance=False, zoom=False, show_title=True, inv=None, ax=None, show_legend=True):
    skippable_pairs_df = add_label_col(
        alu_dct["skippable"], "skippable")
    constitutive_pairs_df = add_label_col(
        alu_dct["constitutive"], "constitutive")
    skippable_inverted_pairs_df = add_label_col(
        alu_dct["skippable_inverted"], "skippable_inverted")
    constitutive_inverted_pairs_df = add_label_col(
        alu_dct["constitutive_inverted"], "constitutive_inverted")
    skippable_noninverted_pairs_df = add_label_col(
        alu_dct["skippable_noninverted"], "skippable_noninverted")
    constitutive_noninverted_pairs_df = add_label_col(
        alu_dct["constitutive_noninverted"], "constitutive_noninverted")

    concat_pairs_df = pd.concat([skippable_pairs_df, constitutive_pairs_df])
    concat_inverted_pairs_df = pd.concat([skippable_inverted_pairs_df, constitutive_inverted_pairs_df])
    concat_noninverted_pairs_df = pd.concat([skippable_noninverted_pairs_df, constitutive_noninverted_pairs_df])
    
    concat_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len", 
                     "intra_alu_len","alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_pairs_df, pairs=True))
    concat_inverted_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len", 
                              "intra_alu_len", "alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_inverted_pairs_df, pairs=True))
    concat_noninverted_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len",
                                 "intra_alu_len", "alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_noninverted_pairs_df, pairs=True))
    
        
    #subset_df = df[df['type'].isin(['skippable', 'const'])]
    # generate_violin_plot(df=concat_pairs_df[concat_pairs_df["label"].isin(["skippable", "constitutive"])], 
    #                      window_size=window, title="All", fill_colour="lightgrey")
    # generate_violin_plot(df=concat_inverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_inverted", "constitutive_inverted"])], 
    #                      window_size=window, title="Inverted", fill_colour="steelblue")
    # generate_violin_plot(df=concat_noninverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_noninverted", "constitutive_noninverted"])], 
    #                      window_size=window, title="Non-inverted", fill_colour="orange")
    if not exon_distance:
        generate_violin_plot_multi(
            df_all=concat_pairs_df[concat_pairs_df["label"].isin(["skippable", "constitutive"])],
            df_inv=concat_inverted_pairs_df,
            df_noninv=concat_noninverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_noninverted", "constitutive_noninverted"])], 
            window_size=window, 
            title="",
            show_legend=show_legend,
            ax=ax)
    
    if exon_distance:
        #generate_histogram(df=concat_pairs_df[concat_pairs_df["label"].isin(["skippable", "constitutive"])],
        #                   window_size=window, title="All")
        if inv == "inverted" and show_title == True:
            generate_histogram_multi(df=concat_inverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_inverted", "constitutive_inverted"])], 
                                     window_size=window, title="Inverted", 
                                     legend=False, ylabel=True, ax=ax, zoom=zoom, show_title=show_title)
        elif inv == "inverted" and show_title == False:
            generate_histogram_multi(df=concat_inverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_inverted", "constitutive_inverted"])], 
                                     window_size=window, title="Inverted", 
                                     legend=True, ylabel=True, ax=ax, 
                                     zoom=zoom, show_title=show_title, 
                                     legend_pos="single")
        elif inv == "non-inverted" and show_title==True:
            generate_histogram_multi(df=concat_noninverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_noninverted", "constitutive_noninverted"])], 
                                     window_size=window, title="Non-inverted", 
                                     legend=True, ylabel=False, ax=ax, zoom=zoom, show_title=show_title)
        elif inv == "non-inverted" and show_title==False:
            generate_histogram_multi(df=concat_noninverted_pairs_df, #[concat_pairs_df["label"].isin(["skippable_noninverted", "constitutive_noninverted"])], 
                                     window_size=window, title="Non-inverted", 
                                     legend=True, ylabel=True, ax=ax, 
                                     zoom=zoom, show_title=show_title, 
                                     legend_pos="single")


def generate_violin_plot_multi(df_all, df_inv, df_noninv, window_size, title, show_legend=True, ax=None):
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    blue="#0077BB"
    light_blue="#86B4CF"
    orange="#EE7733"
    light_orange="#EEC5AE"

    # df_all["hue_label"] = "All"
    # df_inv["hue_label"] = "Inverted"
    # df_noninv["hue_label"] = "Non-inverted"
    ##hue_order = ["Inverted", "Non-inverted", "All"]
    df_all["hue_label"] = "Inv and\nNon-inv"
    df_inv["hue_label"] = "Inv"
    df_noninv["hue_label"] = "Non-inv"
    
    df_plot = pd.concat([df_all, df_inv, df_noninv], axis=0, ignore_index=True)
    # Define conditions and corresponding values
    conditions = [
        df_plot['label'].isin(["skippable", "skippable_inverted", "skippable_noninverted"]),
        df_plot['label'].isin(["constitutive", "constitutive_inverted", "constitutive_noninverted"])
    ]
    values = ["Skippable (S)", "Constitutive (C)"]
    df_plot['group'] = np.select(conditions, values, default=np.nan)

    if not ax:
        fig, ax = plt.subplots(figsize=(4,6))

    fontsize_min = 5
    #plt.figure(figsize=(4, 6))
    sns.boxplot(x='hue_label', y='exon_len', data=df_plot, showcaps=False, 
                #boxprops={'facecolor': fill_colour}, 
                hue="group", #hue_order and palette
                order=["Inv", "Non-inv", "Inv and\nNon-inv"],
                palette=[blue, orange],
                #linecolor=["#0077BB", "#EE7733"],
                showfliers=False, 
                whiskerprops={'linewidth': 0.5},
                linewidth=0.25, 
                ax=ax)
    
    ax.set_ylim(-20, 350)
    ax.set_yticks(np.arange(0, 351, 50))

    ax.set_xlabel("") # Window
    ax.set_ylabel("Exon length (bp)", fontsize=fontsize_min, labelpad=axis_label_pad)

    ax.tick_params(axis='x', rotation=90, labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    
    if show_legend:
        # ax.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.25), 
        #           ncol=1, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
        ax.legend(title='', loc='lower center', bbox_to_anchor=(0.5, -1.25),
                  ncol=1, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
    else:
        ax.get_legend().remove()

    for spine in ax.spines.values():
        spine.set_linewidth(0.25)

    if not ax:
        plt.show()
    #plt.title(f'{title}: {window_size}')
    #plt.xlabel('Category')
    #plt.ylabel('Exon length (bp)')
    #plt.xticks(rotation=45, ha='center')
    #plt.ylim(-20,400)

    #plt.legend(title='', loc='upper right', bbox_to_anchor=(1, 1))
    #plt.legend(title='', loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.legend(title='', ncol=2, loc='upper center', bbox_to_anchor=(0.5, 1.1))
    #plt.xlabel("")
    #plt.show()


def generate_violin_plot(df, window_size, title, fill_colour):
    print(df.head())
    plt.figure(figsize=(3, 6))
    #sns.violinplot(x='label', y='exon_len', data=df)
    # sns.stripplot(x='label', y='exon_len', data=df, jitter=True, color='black', alpha=0.5)
    sns.boxplot(x='label', y='exon_len', data=df, showcaps=False, boxprops={'facecolor': fill_colour}, 
                 showfliers=False, whiskerprops={'linewidth': 2})
    #sns.kdeplot(data=df, x="exon_len", hue="label", fill=True, color="blue")
    plt.title(f'{title}: {window_size}')
    plt.xlabel('Category')
    plt.ylabel('Exon Length (bp)')
    plt.xticks(rotation=45, ha='center')
    #plt.yscale('log') 
    #plt.xscale("log")
    plt.show()


def generate_histogram_multi(df, window_size, title, legend, ylabel, ax=None, zoom=False, show_title=True, legend_pos=None):
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    blue="#0077BB"
    orange="#EE7733"

    fontsize_min = 5

    if zoom:
        n_bins=500
        ax.set_xlim(-500,500)
        ax.set_ylim(0,2800)
        ax.set_xticks([-500, -400, -300, -200, -100, 0, 100, 200, 300, 400, 500])
        ax.set_xticklabels(ax.get_xticks(), fontsize=fontsize_min)
        ax.set_yticks([0, 500, 1000, 1500, 2000, 2500])
        ax.set_yticklabels(ax.get_yticks(), fontsize=fontsize_min)
    else:
        n_bins=100
        ax.set_xlim(-5000,5000)
        ax.set_ylim(0,10500)
        ax.set_xticks([-5000, -2500, 0, 2500, 5000])
        ax.set_xticklabels(ax.get_xticks(), fontsize=fontsize_min)
        ax.set_yticks([0, 2500, 5000, 7500, 10000])
        ax.set_yticklabels(ax.get_yticks(), fontsize=fontsize_min)
    ax.tick_params(axis='both', which='both', 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)

    df_melted = pd.melt(df, id_vars=['label'], value_vars=['upstream_dist', 'downstream_dist'], 
                        value_name='distance')
    #print(df_melted.columns)
    
    label_mapping = {"skippable_inverted": "Skippable",
                     "skippable_noninverted": "Skippable",
                     "constitutive_inverted": "Constitutive", 
                     "constitutive_noninverted": "Constitutive"}
    df_melted["label"] = df_melted["label"].map(label_mapping)
    
    if not ax:
        fig, ax = plt.subplots(figsize=(4,6))
    
    #sns.set_style('ticks')
    plt.rcParams['font.family'] = 'DejaVu Sans'

    sns.histplot(data=df_melted, x='distance', hue='label', kde=False, 
                 hue_order=["Skippable", "Constitutive"],
                 alpha=0.7, multiple="stack", ax=ax, 
                 palette=[blue, orange], legend=legend, bins=n_bins)

    ax.set_xlabel("Distance to exon boundary (bp)", fontsize=fontsize_min, labelpad=axis_label_pad) #12
    #ax.tick_params(axis='x', rotation=45, labelsize=12)
    
    if ylabel:
        ax.set_ylabel("Frequency", fontsize=fontsize_min, labelpad=axis_label_pad)
    else:
        ax.set_ylabel("")

    if show_title:
        ax.set_title(f"{title}", fontsize=fontsize_min)
    else:
        ax.set_title("")
    ax.get_yaxis().set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))

    # ax.set_xticklabels(ax.get_xticklabels(), fontsize=fontsize_min)
    # ax.set_yticklabels(ax.get_yticklabels(), fontsize=fontsize_min)
    # ax.tick_params(axis='both', which='both', length=1, pad=2)

    if legend and not legend_pos:
        #sns.move_legend(ax, loc="lower center", bbox_to_anchor=(0.5, -0.5), 
        #            ncol=2, title='', fancybox=True, shadow=False, borderaxespad=0.)
        sns.move_legend(ax, loc="center right", bbox_to_anchor=(0.985, 0.9), #(1.26, 0.5)
                        ncol=1, title='', fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
    elif legend and legend_pos=="single":
        #sns.move_legend(ax, loc="center right", bbox_to_anchor=(0.97, 0.85),
        #                ncol=1, title='', fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
        sns.move_legend(ax, loc="upper center", bbox_to_anchor=(0.5, 1.15),  # Adjust y-position as needed
                        ncol=2, title='', fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)

    for spine in ax.spines.values():
        spine.set_linewidth(0.25)


def generate_histogram(df, window_size, title):
    print(df.head())
    ['upstream_dist', 'downstream_dist']
    df_melted = pd.melt(df, id_vars=['label'], value_vars=['upstream_dist', 'downstream_dist'], 
                        value_name='distance')
    print(df_melted.columns)
    plt.figure(figsize=(8, 3))
    ax = sns.histplot(data=df_melted, x='distance', hue='label', kde=False, alpha=0.7, multiple="stack")

    plt.xlabel('Distance to exon (bp)')
    plt.ylabel('Frequency')

    sns.move_legend(ax, loc="lower center", bbox_to_anchor=(0.5, -0.5), ncol=2, title='Exon group')
    #plt.legend(handles=handles, labels=labels, title='Exon group', loc='upper center', bbox_to_anchor=(0.5, -0.1), shadow=True, ncol=2)
    plt.title(f"{title}: {window_size}")
    plt.show()
    

def run_pairplot_window(alu_dct, window, ks_test=False, cv=0, n_sample=0, inverted=False):
    """
    window: '5kb'
    """
    pairgrid_pair_cols = ["upstream_alu_len", "downstream_alu_len", "exon_len", 
                          "upstream_dist", "downstream_dist", "label"]

    skippable_pairs_df = (
        add_label_col(alu_dct["skippable"], "skippable")
        )
    constitutive_pairs_df = (
        add_label_col(alu_dct["constitutive"], "constitutive")
        )
    skippable_inverted_pairs_df = (
        add_label_col(alu_dct["skippable_inverted"], "skippable_inverted")
        )
    constitutive_inverted_pairs_df = (
        add_label_col(alu_dct["constitutive_inverted"], "constitutive_inverted")
        )
    # new
    skippable_noninverted_pairs_df = (
        add_label_col(alu_dct["skippable_noninverted"], "skippable_noninverted")
        )
    constitutive_noninverted_pairs_df = (
        add_label_col(alu_dct["constitutive_noninverted"], "constitutive_noninverted")
        )

    concat_pairs_df = pd.concat([skippable_pairs_df, constitutive_pairs_df])
    concat_inverted_pairs_df = pd.concat([skippable_inverted_pairs_df, constitutive_inverted_pairs_df])
    # new
    concat_noninverted_pairs_df = pd.concat([skippable_noninverted_pairs_df, constitutive_noninverted_pairs_df])

    # compute lengths of Alus and exons
    #concat_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len"]] = (
    #   add_alu_len_exon_len_cols(concat_pairs_df, pairs=True))
    #concat_inverted_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len"]] = (
    #   add_alu_len_exon_len_cols(concat_inverted_pairs_df, pairs=True))
    
    concat_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len", 
                     "intra_alu_len","alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_pairs_df, pairs=True))
    concat_inverted_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len", 
                              "intra_alu_len", "alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_inverted_pairs_df, pairs=True))
    # new
    concat_noninverted_pairs_df[["upstream_alu_len", "downstream_alu_len", "exon_len", 
                              "intra_alu_len", "alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]] = (
       add_alu_len_exon_len_cols(concat_noninverted_pairs_df, pairs=True))
    
    cols_display = ["exon_chr", "exon_start", "exon_end", 
                    "upstream_alu_start", "upstream_alu_end", "downstream_alu_start", "downstream_alu_end",
                    "intra_alu_len", "alu1_end_to_alu2_start", "alu2_end_to_alu1_start"]
    #display(concat_pairs_df[cols_display].head())
    #display(concat_inverted_pairs_df[cols_display].head())
    #display(concat_noninverted_pairs_df[cols_display].head())

    # run KS test
    if ks_test:
        ks_res_dct = run_ks_test(df1=concat_pairs_df, 
                                 df2_inverted=concat_inverted_pairs_df,
                                 df3_noninverted=concat_noninverted_pairs_df,
                                 cv=cv,
                                 n_sample=n_sample)
        return ks_res_dct

    # count number of Alus in each category
    #display(concat_pairs_df["label"].value_counts())
    #display(concat_inverted_pairs_df["label"].value_counts())
    #display(concat_noninverted_pairs_df["label"].value_counts())

    # plot
    if not inverted:
        generate_pairgrid_plot(df=concat_pairs_df,
                               cols=pairgrid_pair_cols,
                               nrows_sample=10000,
                               window_size=window,
                               title="Pairs", 
                               inversion=None)
    elif inverted:
        generate_pairgrid_plot(df=concat_inverted_pairs_df,
                               cols=pairgrid_pair_cols,
                               nrows_sample=10000,
                               window_size=window,
                               title="Inverted pairs", 
                               inversion="inverted")
    # re-run and output 1 plot per cell to save
    # re-run and save plot to file


def plot_ks_result(ks_dct, window, cv=0, n_sample=0, ax=None, show_legend=True, ax_legend=None):
    """
    Expects a dictionary with "pval" and "stat" entries (resulting from a 2 sample ks test).
    """
    # Create dataframe from nested dictionary of dictionaries
    # Ref: https://stackoverflow.com/questions/13575090
    ks_df = pd.DataFrame.from_dict({(i,j): ks_dct[i][j] 
                                    for i in ks_dct.keys()
                                    for j in ks_dct[i].keys()},
                       orient='index').reset_index().rename(columns={"level_0": "feature", "level_1": "fixed_category"})
    if cv>0:
        #num_pvalues = ks_df['pval'].apply(len).max()
        pval_columns = [f'pval{i+1}' for i in range(cv)]
        negpval_columns = [f'neglog10pval{i+1}' for i in range(cv)]
        stat_columns = [f'stat{i+1}' for i in range(cv)]
        #ks_df[pval_columns] = pd.DataFrame(ks_df['pval'].tolist(), index=ks_df.index) # 240423 added multiple testing correction
        pvals = ks_df['pval'].tolist()
        pvals_corrected = []
        for pval_list in pvals:
            _, bonferroni_corrected_p_values, _, _ = multipletests(pval_list, method='bonferroni')
            pvals_corrected.append(bonferroni_corrected_p_values)
        ks_df[pval_columns] = pd.DataFrame(pvals_corrected, index=ks_df.index)
        ks_df[stat_columns] = pd.DataFrame(ks_df['stat'].tolist(), index=ks_df.index)
        #display(ks_df)
        for pval_col in pval_columns:
            prob_temp = np.where(ks_df[pval_col] > 1.0e-10, ks_df[pval_col], 1.0e-10)
            ks_df[f"neglog10{pval_col}"] = -np.log10(prob_temp)
        #return ks_df
        id_vars = ["feature", "fixed_category", "pval", "stat"] + stat_columns + pval_columns
        ks_df_melt = pd.melt(ks_df, id_vars=id_vars, value_vars=negpval_columns, value_name="neglog10pvals")
        #display(ks_df_melt.groupby(["feature", "fixed_category"])["neglog10pvals"].nunique())

        # category_mapping = {'skippable': 'Skippable:\nInverted vs Non-inverted', 
        #                     'constitutive': 'Constitutive:\nInverted vs Non-inverted', 
        #                     'inverted': 'Inverted:\nSkippable vs Constitutive', 
        #                     'noninverted': 'Non-inverted:\nSkippable vs Constitutive', 
        #                     'all': 'All:\nSkippable vs Constitutive'}

        category_mapping = {'skippable': 'Skip:\nInv vs Non-inv', 
                            'constitutive': 'Const:\nInv vs Non-inv', 
                            'inverted': 'Inv:\nSkip vs Const', 
                            'noninverted': 'Non-inv:\nSkip vs Const', 
                            'all': 'Inv and Non-inv:\nSkip vs Const'}
        
        def replace_with_subscript(text):
            lines = text.split('\n')
            processed_lines = []
            for line in lines:
                parts = line.split(' vs ')
                new_parts = []
                for part in parts:
                    subparts = part.split('_')
                    subscripted_text = f'{subparts[0]}'
                    for subpart in subparts[1:]:
                        subscripted_text += f'$_{{\mathrm{{{subpart}}}}}$'
                    new_parts.append(subscripted_text)
                processed_lines.append(' vs '.join(new_parts))
            return '\n'.join(processed_lines)
        
        # def replace_with_subscript(text):
        #     parts = text.split(' vs ')
        #     new_parts = []
        #     for part in parts:
        #         subparts = part.split('_')
        #         subscripted_text = f'{subparts[0]}'
        #         for subpart in subparts[1:]:
        #             subscripted_text += f'$_{{\mathrm{{{subpart}}}}}$'
        #         new_parts.append(subscripted_text)
        #     return ' vs '.join(new_parts)
        
        category_mapping_subscripts = {"skippable": replace_with_subscript("S_inv vs S_non-inv"), 
                                       "constitutive": replace_with_subscript("C_inv vs C_non-inv"),
                                       "inverted": replace_with_subscript("Inv_S vs Inv_C"),
                                       "noninverted": replace_with_subscript("Non-inv_S vs \n Non-inv_C"),
                                       "all": replace_with_subscript("S vs C")}
        
        feature_mapping = {"upstream_alu_len": "Upstream Alu length", 
                           "downstream_alu_len": "Downstream Alu length",
                           "exon_len": "Exon length", 
                           "upstream_dist": "Upstream Alu distance to exon",
                           "downstream_dist": "Downstream Alu distance to exon", 
                           "intra_alu_len": "Intra-Alu distance"}
        feature_order = ["Upstream Alu length", "Downstream Alu length",
                         "Upstream Alu distance to exon", "Downstream Alu distance to exon",
                         "Intra-Alu distance", "Exon length"]
        ks_df_melt["fixed_category"] = ks_df_melt["fixed_category"].map(category_mapping_subscripts)
        ks_df_melt["feature"] = ks_df_melt["feature"].map(feature_mapping)
        
        # Ref: https://stackoverflow.com/questions/69315871
        if not ax:
            fig, ax1 = plt.subplots(figsize=(12, 6))
        
        colorblind_palette = sns.color_palette('colorblind')
        selected_colors = colorblind_palette[2:8]
        selected_colors = selected_colors[1:] + [selected_colors[0]]

        fontsize_min = 5
        #sns.set_theme() #font_scale=1.4
        #sns.set_style('ticks')
        plt.rcParams['font.family'] = 'DejaVu Sans'
        
        tick_params_width = 0.25
        tick_params_length = 1
        tick_params_pad = 2
        axis_label_pad = 2

        sns.barplot(data=ks_df_melt, x='fixed_category', y='neglog10pvals', hue="feature", 
                    palette=selected_colors, hue_order=feature_order, ax=ax, edgecolor='none', 
                    err_kws={'linewidth': 0.5}) #viridis
        sns.stripplot(data=ks_df_melt, x='fixed_category', y='neglog10pvals', hue="feature", 
                      palette=selected_colors, dodge=True, alpha=0.6, legend=False, 
                      hue_order=feature_order, edgecolor="black", linewidth=0.25, ax=ax, size=2) #edgecolor="black"
        if ax:
            ax.axhline(-np.log(0.05), ls='--', linewidth=0.5)
            ax.tick_params(axis='x', rotation=90, labelsize=fontsize_min, 
                           pad=tick_params_pad, length=tick_params_length, width=tick_params_width) #14
            ax.tick_params(axis='y', labelsize=fontsize_min, 
                           pad=tick_params_pad, length=tick_params_length, width=tick_params_width) #14
            ax.set_xlabel("") # Window
            ax.set_ylabel("-log10(p_adj)", fontsize=fontsize_min, labelpad=axis_label_pad) #16
            handles, labels = ax.get_legend_handles_labels()
            if show_legend and ax_legend:
                # ax_legend.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.25), #(0.5,1.12)
                #           ncol=2, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min) #14
                ax_legend.legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.1),
                                 ncol=2, fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)
                ax_legend.axis('off')
                ax.get_legend().remove()
            else:
                ax.get_legend().remove()
            for spine in ax.spines.values():
                spine.set_linewidth(0.25)
        
        if not ax:
            ax1.axhline(-np.log(0.05), ls='--')
            ax1.legend(title='', loc='upper center', bbox_to_anchor=(0.5, 1.25), #(0.5,1.12)
                      ncol=2, fancybox=True, shadow=False, borderaxespad=0.)
            plt.xticks(rotation=45, ha='right')
            plt.xlabel('Category')
            plt.ylabel('-log10 p-value')
            plt.title(f'{window} flanking window ({cv}-fold CV randomly choosing {n_sample} samples): KS test p-values')
            plt.tight_layout()
            plt.show()

        return 
    
    # Convert p-values to -log10(p-values) and account for dividing by 0
    # Ref: https://stackoverflow.com/questions/21610198
    prob_temp = np.where(ks_df["pval"] > 1.0e-10, ks_df["pval"], 1.0e-10)
    #ks_df_temp["neglog10pval"] = np.where(prob_temp > 1.0e-10, -np.log10(prob_temp), 10)
    ks_df["neglog10pval"] = -np.log10(prob_temp)
    display(ks_df)

    plt.figure(figsize=(12, 6))
    ax = sns.barplot(data=ks_df, x='fixed_category', y='neglog10pval', hue="feature", palette='viridis', errorbar=None)
    ax.axhline(-np.log(0.05), ls='--')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Category')
    plt.ylabel('-log10 p-value')
    plt.title(f'{window} flanking window: KS test p-values')
    plt.tight_layout()
    plt.show()