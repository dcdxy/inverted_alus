#!/usr/bin/env python
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib_venn import venn2,venn3
import matplotlib.ticker as ticker
import matplotlib as mpl
from scipy.stats import ks_2samp, ttest_ind

from generate_exon_alu_df import groupby_exons
from mfe import add_alu_cols

hue_order = ["skipped", "constitutive"]
hue_order2 = ["inverted", "non-inverted"]
palette = {
    'inverted': 'tab:blue',
    'skipped': 'tab:red',
    'non-inverted': 'tab:orange',
    'constitutive': 'tab:green'
}
palette_random = {'inverted': 'tab:olive', 'non-inverted': 'tab:grey'}


def generate_kdeplot_4sets_deprecated(df_skip, df_const):
    """ 
    df_skip: skipped_exonskipdb
    df_const: constitutive_exonskipdb
    """
    fig, ax = plt.subplots(figsize=(10,6))
    sns.kdeplot(data=df_skip[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
               hue_order=hue_order2, legend=False)
    sns.kdeplot(data=df_const[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
                hue_order=hue_order2, linestyle='--', legend=False)
    line = Line2D([0,1],[0,1],linestyle='-', color='tab:blue') #inverted, skipped
    line2 = Line2D([0,1],[0,1],linestyle='-', color='tab:orange') #non-inverted, skipped
    line3 = Line2D([0,1],[0,1],linestyle='--', color='tab:blue') #inverted, constitutive
    line4 = Line2D([0,1],[0,1],linestyle='--', color='tab:orange') # non-inverted, constitutive
    ax.legend([line, line3, line2, line4],["inverted, skipped", "inverted, constitutive",
                                           "non-inverted, skipped", "non-inverted, constitutive"])
    plt.show()


def generate_kdeplot_5sets(df_inv, df_noninv, df_randinv, df_randnoninv, title, ax=None): # 240426: df_skip, df_const, df_randskip, df_randconst
    """ 
    df_skip: skipped_exonskipdb
    df_const: constitutive_exonskipdb
    """
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    fontsize_min = 5
    plt.rcParams['font.family'] = 'DejaVu Sans'

    palette = {
    #'inverted': 'tab:blue',
    'skippable': 'tab:blue',
    #'non-inverted': 'tab:orange',
    'constitutive': 'tab:orange'
    }
    
    blue="#0077BB"
    light_blue="#86B4CF"
    orange="#EE7733"
    light_orange="#EEC5AE"
    
    palette_random = {'inverted': '#CCCCCC', 'non-inverted': '#CCCCCC', 
                      "skippable": '#CCCCCC', "constitutive": '#CCCCCC'}  # 'tab:grey'
    
    hue_order2 = ["inverted", "non-inverted"]
    hue_order3 = ["skippable", "constitutive"]
    
    mapping_inv = {"skippable_inverted": "skippable", "constitutive_inverted": "constitutive"}
    df_inv["exon_set"] = df_inv["type"].map(mapping_inv)
    
    mapping_noninv = {"skippable_noninverted": "skippable", "constitutive_noninverted": "constitutive"}
    df_noninv["exon_set"] = df_noninv["type"].map(mapping_noninv)
    
    df_randinv["exon_set"] = df_randinv["type"].map(mapping_inv)
    df_randnoninv["exon_set"] = df_randnoninv["type"].map(mapping_noninv)

    sns.set_theme() #font_scale=1.5
    sns.set_style('ticks')

    #fig, ax = plt.subplots(figsize=(10,6))
    
    # sns.kdeplot(data=df_inv[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
    #            hue_order=hue_order2, style="inversion", legend=False, palette=palette, label=["Skippable"]) #label=["Skippable", "b"]
    # sns.kdeplot(data=df_noninv[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
    #             hue_order=hue_order2, linestyle='--', legend=False, palette=palette, label="Constitutive")
    # sns.kdeplot(data=df_randskip[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
    #            hue_order=hue_order2, linestyle='-', legend=False, palette=palette_random, label="Random (Skippable)")
    # sns.kdeplot(data=df_randconst[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
    #            hue_order=hue_order2, linestyle="-", legend=False, palette=palette_random, label="Random (Constitutive)")
    
    # XXX split into individual distributions
    # sns.kdeplot(data=df_inv[['seq_corrected_mfe', 'exon_set']], x="seq_corrected_mfe", hue="exon_set",
    #            hue_order=hue_order3, legend=True, palette=palette, label="Inverted") #label=["Skippable", "b"]
    # sns.kdeplot(data=df_noninv[['seq_corrected_mfe', 'exon_set']], x="seq_corrected_mfe", hue="exon_set",
    #             hue_order=hue_order3, linestyle='--', legend=True, palette=palette, label="Non-inverted")
    # sns.kdeplot(data=df_randinv[['seq_corrected_mfe', 'exon_set']], x="seq_corrected_mfe", hue="exon_set",
    #            hue_order=hue_order3, linestyle='-', legend=True, palette=palette_random, label="Random (Skippable)")
    # sns.kdeplot(data=df_randnoninv[['seq_corrected_mfe', 'exon_set']], x="seq_corrected_mfe", hue="exon_set",
    #            hue_order=hue_order3, linestyle="-", legend=True, palette=palette_random, label="Random (Constitutive)")
    
    sns.kdeplot(data=df_inv[df_inv["type"] == "skippable_inverted"][['seq_corrected_mfe', "type"]], 
                x="seq_corrected_mfe", 
                #color="tab:blue",
                color=blue,
                legend=True, 
                #hue="exon_set",
                #hue_order=hue_order3,
                #palette=palette,  
                label="Skippable inverted", ax=ax, linewidth=0.5) #label=["Skippable", "b"] XXX ax is new
    sns.kdeplot(data=df_inv[df_inv["type"] == "constitutive_inverted"][['seq_corrected_mfe', "type"]], 
                x="seq_corrected_mfe", 
                #color="tab:orange",
                color=orange,
                legend=True, 
                #hue="exon_set",
                #hue_order=hue_order3,
                #palette=palette, 
                label="Constitutive inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(data=df_noninv[df_noninv["type"] == "skippable_noninverted"][['seq_corrected_mfe', 'type']], 
                x="seq_corrected_mfe",
                #color="tab:blue",
                color=light_blue,
                linestyle='--', 
                legend=True,
                #hue="exon_set",
                #hue_order=hue_order3, 
                #palette=palette, 
                label="Skippable non-inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(data=df_noninv[df_noninv["type"] == "constitutive_noninverted"][['seq_corrected_mfe', 'type']], 
                x="seq_corrected_mfe",
                #color="tab:orange",
                color=light_orange,
                linestyle='--', 
                legend=True,
                #hue="exon_set",
                #hue_order=hue_order3, 
                #palette=palette, 
                label="Skippable non-inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_randinv[df_randinv["type"] == "skippable_inverted"][['seq_corrected_mfe', 'type']],
        x="seq_corrected_mfe",
        linestyle='-', 
        legend=True,
        color="#CCCCCC",
        #hue="exon_set",
        #hue_order=hue_order3,
        #palette=palette_random,
        label="Random (Skippable inverted)", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_randinv[df_randinv["type"] == "constitutive_inverted"][['seq_corrected_mfe', 'type']],
        x="seq_corrected_mfe",
        linestyle='-', 
        legend=True,
        color="#CCCCCC",
        #hue="exon_set",
        #hue_order=hue_order3,
        #palette=palette_random,
        label="Random (Constitutive inverted)", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_randnoninv[df_randnoninv["type"] == "skippable_noninverted"][['seq_corrected_mfe', 'type']],
        x="seq_corrected_mfe",
        linestyle='-', 
        legend=True,
        color="#CCCCCC",
        #hue="exon_set",
        #hue_order=hue_order3,
        #palette=palette_random,
        label="Random (Skippable non-inverted)", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_randnoninv[df_randnoninv["type"] == "constitutive_noninverted"][['seq_corrected_mfe', 'type']],
        x="seq_corrected_mfe",
        linestyle='-', 
        legend=True,
        color="#CCCCCC",
        #hue="exon_set",
        #hue_order=hue_order3,
        #palette=palette_random,
        label="Random (Constitutive non-inverted)", ax=ax, linewidth=0.5)
    
    line1 = Line2D([0,1],[0,1],linestyle='-', color=blue) #skippable inverted
    line2 = Line2D([0,1],[0,1],linestyle='--', color=light_blue, dashes=[1,1]) #skippable non-inverted
    line3 = Line2D([0,1],[0,1],linestyle='-', color=orange) #constitutive inverted
    line4 = Line2D([0,1],[0,1],linestyle='--', color=light_orange, dashes=[1,1]) # constitutive non-inverted
    line8 = Line2D([0,1],[0,1],linestyle='-', color='#CCCCCC') # random constitutive non-inverted
    
    #line1 = Line2D([0,1],[0,1],linestyle='-', color='tab:blue') #inverted, skipped
    #line2 = Line2D([0,1],[0,1],linestyle='-', color='tab:orange') #non-inverted, skipped
    #line3 = Line2D([0,1],[0,1],linestyle='--', color='tab:blue') #inverted, constitutive
    #line4 = Line2D([0,1],[0,1],linestyle='--', color='tab:orange') # non-inverted, constitutive
    
    #line5 = Line2D([0,1],[0,1],linestyle=':', color='tab:olive') #inverted, randskip
    #line6 = Line2D([0,1],[0,1],linestyle=':', color='tab:grey') # non-inverted, randskip
    #line7 = Line2D([0,1],[0,1],linestyle='-.', color='tab:olive') #inverted, randconst
    #line8 = Line2D([0,1],[0,1],linestyle='-', color='#CCCCCC') # non-inverted, randconst

    legend_lines  = [line1, line2, line3, line4, line8] 
    legend_labels = ["Skippable inverted", "Skippable non-inverted",
                     "Constitutive inverted", "Constitutive non-inverted",
                     "Random"]
    #[*zip(legend_lines, legend_labels)]
    #ax.legend(legend_lines, legend_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #g.fig.legend(handles, labels, loc='upper center', ncol=3, fontsize=14, bbox_to_anchor=(0.5, 1.05))
    ax.legend(legend_lines, legend_labels, loc="upper center", bbox_to_anchor=(0.5, 1.16), 
              ncol=3, borderaxespad=0., fontsize=fontsize_min, 
              columnspacing=0.5, borderpad=0.5, handletextpad=0.4) #handlelength=0.6

    #ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #plt.title(title)
    #plt.xlabel("Length-adjusted MFE")
    ax.set_xlabel("Length-adjusted MFE (kcal/mol)", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax.set_ylabel("Density", fontsize=fontsize_min, labelpad=axis_label_pad)

    ax.tick_params(axis='x', labelsize=fontsize_min, pad=tick_params_pad, 
                   length=tick_params_length, width=tick_params_width) #14
    ax.tick_params(axis='y', labelsize=fontsize_min, pad=tick_params_pad, 
                   length=tick_params_length, width=tick_params_width) #14
    
    ax.xaxis.grid(True, linewidth=0.25)
    #ax.xaxis.grid(True, which='both', linestyle='--', linewidth=0.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    #plt.show()


def generate_kdeplot_4sets(df_inv, df_noninv, ax, show_legend=True):
    """ 
    df_skip: skipped_exonskipdb
    df_const: constitutive_exonskipdb
    """
    import matplotlib.ticker as ticker

    fontsize_min = 5
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    plt.rcParams['font.family'] = 'DejaVu Sans'

    blue="#0077BB"
    light_blue="#86B4CF"
    orange="#EE7733"
    light_orange="#EEC5AE"
    
    mapping_inv = {"skippable_inverted": "skippable", "constitutive_inverted": "constitutive"}
    df_inv["exon_set"] = df_inv["type"].map(mapping_inv)
    
    mapping_noninv = {"skippable_noninverted": "skippable", "constitutive_noninverted": "constitutive"}
    df_noninv["exon_set"] = df_noninv["type"].map(mapping_noninv)
    
    sns.set_theme() #font_scale=1.5
    sns.set_style('ticks')
    #fig, ax = plt.subplots(figsize=(10,6))
    
    sns.kdeplot(
        data=df_inv[df_inv["type"] == "skippable_inverted"][
            ['seq_corrected_mfe', "type"]], 
        x="seq_corrected_mfe",
        color=blue,
        legend=True,
        label="Skippable inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_inv[df_inv["type"] == "constitutive_inverted"][
            ['seq_corrected_mfe', "type"]],
        x="seq_corrected_mfe",
        color=orange,
        legend=True, 
        label="Constitutive inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_noninv[df_noninv["type"] == "skippable_noninverted"][
            ['seq_corrected_mfe', 'type']], 
        x="seq_corrected_mfe",
        color=light_blue,
        linestyle='--', 
        legend=True,
        label="Skippable non-inverted", ax=ax, linewidth=0.5)
    sns.kdeplot(
        data=df_noninv[df_noninv["type"] == "constitutive_noninverted"][
            ['seq_corrected_mfe', 'type']], 
        x="seq_corrected_mfe",
        color=light_orange,
        linestyle='--', 
        legend=True, 
        label="Skippable non-inverted", ax=ax, linewidth=0.5)
    
    line1 = Line2D([0,1],[0,1],linestyle='-', color=blue) #skippable inverted
    line2 = Line2D([0,1],[0,1],linestyle='--', color=light_blue, dashes=[1,1]) #skippable non-inverted
    line3 = Line2D([0,1],[0,1],linestyle='-', color=orange) #constitutive inverted
    line4 = Line2D([0,1],[0,1],linestyle='--', color=light_orange, dashes=[1,1]) # constitutive non-inverted
    
    legend_lines  = [line1, line2, line3, line4] 
    legend_labels = ["Skippable inverted", "Skippable non-inverted",
                     "Constitutive inverted", "Constitutive non-inverted"]
    if show_legend:
        #ax.legend(legend_lines, legend_labels, loc="upper center",
        #          bbox_to_anchor=(0.5, 1.16), ncol=2, borderaxespad=0., fontsize=14)
        ax.legend(legend_lines, legend_labels, loc="upper center", bbox_to_anchor=(0.5, 1.3), 
                  ncol=2, borderaxespad=0., fontsize=fontsize_min, 
                  columnspacing=0.5, borderpad=0.5, handletextpad=0.4)
    else:
        legend = ax.get_legend()
        if legend:
            legend.remove()

    ax.set_xlabel("Length-adjusted MFE (kcal/mol)", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax.set_ylabel("Density", fontsize=fontsize_min, labelpad=axis_label_pad)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(10))

    ax.tick_params(axis='both', labelsize=fontsize_min, pad=tick_params_pad, 
                   length=tick_params_length, width=tick_params_width)
    ax.xaxis.grid(True, linewidth=0.25)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)   
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
    #ax.xaxis.grid(True, which='both', linestyle='--', linewidth=0.5)


def generate_kdeplot_6sets(df_skip, df_const, df_randskip, df_randconst, title):
    """ 
    df_skip: skipped_exonskipdb
    df_const: constitutive_exonskipdb
    """
    sns.set(font_scale=1.5)
    sns.set_style('ticks')
    fig, ax = plt.subplots(figsize=(10,6))
    sns.kdeplot(data=df_skip[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
               hue_order=hue_order2, legend=False, palette=palette, label=["Skippable", "b"])
    sns.kdeplot(data=df_const[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
                hue_order=hue_order2, linestyle='--', legend=False, palette=palette, label="Constitutive")
    sns.kdeplot(data=df_randskip[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
               hue_order=hue_order2, linestyle=':', legend=False, palette=palette_random, label="Random (Skippable)")
    sns.kdeplot(data=df_randconst[['seq_corrected_mfe', 'inversion']], x="seq_corrected_mfe", hue="inversion",
               hue_order=hue_order2, linestyle="-.", legend=False, palette=palette_random, label="Random (Constitutive)")
    line1 = Line2D([0,1],[0,1],linestyle='-', color='tab:blue') #inverted, skipped
    line2 = Line2D([0,1],[0,1],linestyle='-', color='tab:orange') #non-inverted, skipped
    line3 = Line2D([0,1],[0,1],linestyle='--', color='tab:blue') #inverted, constitutive
    line4 = Line2D([0,1],[0,1],linestyle='--', color='tab:orange') # non-inverted, constitutive
    line5 = Line2D([0,1],[0,1],linestyle=':', color='tab:olive') #inverted, randskip
    line6 = Line2D([0,1],[0,1],linestyle=':', color='tab:grey') # non-inverted, randskip
    line7 = Line2D([0,1],[0,1],linestyle='-.', color='tab:olive') #inverted, randconst
    line8 = Line2D([0,1],[0,1],linestyle='-.', color='tab:grey') # non-inverted, randconst

    legend_lines  = [line1, line2, line3, line4, line5, line6, line7, line8] 
    legend_labels = ["skipped inverted", "skipped non-inverted",
                     "constitutive inverted", "constitutive non-inverted",
                     "random (skipped inverted)", "random (skipped non-inverted)",
                     "random (constitutive inverted)", "random (constitutive non-inverted)"]
    #[*zip(legend_lines, legend_labels)]
    ax.legend(legend_lines, legend_labels, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    #g.fig.legend(handles, labels, loc='upper center', ncol=3, fontsize=14, bbox_to_anchor=(0.5, 1.05))

    #ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title(title)
    plt.xlabel("Length-adjusted MFE")
    plt.show()


def generate_violinplot_2sets(df):
    """
    Input
    df: concatenated skippable and constitutive df (e.g. df_exonskipdb)
    """
    plt.figure(figsize=(16,6))
    ax = sns.violinplot(data=df[['seq_corrected_mfe', 'type', 'inversion']],
                        y="seq_corrected_mfe",
                        x="type",
                        hue="inversion",
                        hue_order=hue_order2,
                        inner='quartile',
                        split=True,
                        scale="count")


def generate_violinplot_4sets(df):
    """
    Input
    df: concatenated skippable and constitutive df (e.g. df_exonskipdb)
    """
    plt.figure(figsize=(16,6))
    ax = sns.violinplot(data=df[['seq_corrected_mfe', 'type', 'inversion']],
                        y="seq_corrected_mfe",
                        x="type",
                        hue="inversion",
                        hue_order=hue_order2,
                        inner='box') #df[['mfe','type']]


def generate_multiset_df_violinplots(df_skip, df_const, df_randskip, df_randconst, case):
    """
    Case 1: df_skipped vs df_constitutive vs df_rand
    Case 2: df_skipped vs df_constitutive vs df_randskip vs df_randconst
    Case 3: df_skipped vs df_randskip
    """
    df_skip      = df_skip.copy()
    df_const     = df_const.copy()
    df_randskip  = df_randskip.copy()
    df_randconst = df_randconst.copy()
    df_skip["group"] = "skippable"
    df_const["group"] = "constitutive"
    if case == 1:
        df_randskip["group"]  = "random"
        df_randconst["group"] = "random"
        violinplot_df = (pd.concat([df_skip, df_const, df_randskip, df_randconst]).reset_index())
    elif case == 2:
        df_randskip["group"]  = "random (skippable)"
        df_randconst["group"] = "random (constitutive)"
        violinplot_df = (pd.concat([df_skip, df_const, df_randskip, df_randconst]).reset_index())
    elif case == 3:
        df_randskip["group"]  = "random (skippable)"
        violinplot_df = (pd.concat([df_skip, df_randskip]).reset_index())
    return add_alu_cols(violinplot_df)


def generate_violinplot_multisets(df):
    """
    Input
    df: concatenated skippable and constitutive df (e.g. df_exonskipdb)
    """
    sns.set_style("whitegrid")
    plt.figure(figsize=(16,6))
    ax = sns.violinplot(data=df[['seq_corrected_mfe', 'group', 'inversion']],
                        y="seq_corrected_mfe",
                        x="group",
                        hue="inversion",
                        hue_order=hue_order2,
                        inner='box') #df[['mfe','type']]
    ax.set(xlabel=None)
    plt.ylabel("Length-adjusted MFE")
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


def generate_mfe_raw_plot(df, legend=True, legend_pos="on", ax=None):
    """
    Input
    df: concatenated skippable and constitutive df (e.g. df_exonskipdb)
    """
    fontsize_min = 5
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    plt.rcParams['font.family'] = 'DejaVu Sans'
    mpl.rcParams['savefig.dpi'] = 300

    # size of ticks
    # size of grid lines

    sns.set_theme() #font_scale=1.2
    sns.set_style("whitegrid")
    
    
    def format_with_commas(x, pos):
        # ChatGPT answer
        return '{:,.0f}'.format(x)
    
    #sns.histplot(data=df[['mfe', 'exon_type', 'seq_length']],  # 'type'
    #             x="seq_length",
    #             y="mfe",
    #             hue="exon_type",
    #             hue_order=hue_order, stat="count")
    if not ax:
        fig, ax = plt.subplots(figsize=(6,6))
    
    df["seq_length_kb"] = df["seq_length"]/1000
    hue_order = ["skippable", "constitutive"]
    sns.scatterplot(data=df, x="seq_length_kb", y="mfe", hue="exon_type", hue_order=hue_order, s=2, alpha=0.5, ax=ax)

    # Overlay regression lines for each hue category
    for hue_category in hue_order: #df["exon_type"].unique()
        subset_data = df[df["exon_type"] == hue_category]
        sns.regplot(data=subset_data, x="seq_length_kb", y="mfe", scatter=False, ci=None, ax=ax, line_kws={"linewidth": 0.5})
    
    if legend:
        if legend_pos == "upper":
            ax.legend(title="Exon type", loc='upper center', bbox_to_anchor=(0.5, 1.25), 
                      ncol=1, fancybox=True, shadow=False, borderaxespad=0, 
                      fontsize=fontsize_min, title_fontsize=fontsize_min)
            legend = ax.get_legend()
            for text in legend.get_texts():
                text.set_text(text.get_text().title())
        elif legend_pos == "on":
            ax.legend(title="Exon type", loc='upper center', bbox_to_anchor=(0.7, 0.96), 
                      ncol=1, fancybox=True, shadow=False, borderaxespad=0., framealpha=1,
                      fontsize=fontsize_min, title_fontsize=fontsize_min)
            legend = ax.get_legend()
            for text in legend.get_texts():
                text.set_text(text.get_text().title())
    else:
        ax.get_legend().remove()
    
    ax.set_ylim(None,0)
    ax.set_xlim(0,18) #None

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))

    ax.set_xlabel("Sequence length (kbp)", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax.set_ylabel("MFE (kcal/mol)", fontsize=fontsize_min, labelpad=axis_label_pad)
    
    ax.tick_params(axis='x', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    
    for collection in ax.collections:
        collection.set_rasterized(True)
    
    ax.grid(True, linewidth=0.25)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    
    if not ax:
        plt.show()


def generate_mfe_corrected_plot(df, legend=True, legend_pos="on", ax=None):
    """
    Input
    df: concatenated skippable and constitutive df (e.g. df_exonskipdb)
    """
    fontsize_min = 5
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    plt.rcParams['font.family'] = 'DejaVu Sans'
    mpl.rcParams['savefig.dpi'] = 300

    sns.set_theme() #font_scale=1.2
    sns.set_style("whitegrid")
    
    def format_with_commas(x, pos):
        # ChatGPT answer
        return '{:,.0f}'.format(x)
    
    # sns.histplot(data=df[['seq_corrected_mfe', 'exon_type', 'seq_length']],  # 'type'
    #              x="seq_length",
    #              y="seq_corrected_mfe",
    #              hue="exon_type",
    #              hue_order=hue_order)
    if not ax:
        fig, ax = plt.subplots(figsize=(6,6))
    df["seq_length_kb"] = df["seq_length"]/1000
    hue_order = ["skippable", "constitutive"]
    sns.scatterplot(data=df, x="seq_length_kb", y="seq_corrected_mfe", hue="exon_type", hue_order=hue_order, s=2, alpha=0.5, ax=ax)

    # Overlay regression lines for each hue category
    for hue_category in hue_order: #df["exon_type"].unique()
        subset_data = df[df["exon_type"] == hue_category]
        sns.regplot(data=subset_data, x="seq_length_kb", y="seq_corrected_mfe", scatter=False, ci=None, ax=ax, line_kws={"linewidth": 0.5})
    
    if legend:
        if legend_pos == "upper":
            ax.legend(title="Exon type", loc='upper center', bbox_to_anchor=(0.5, 1.25), 
                      ncol=1, fancybox=True, shadow=False, borderaxespad=0., 
                      fontsize=fontsize_min, title_fontsize=fontsize_min)
            legend = ax.get_legend()
            for text in legend.get_texts():
                text.set_text(text.get_text().title())
        elif legend_pos == "on":
            ax.legend(title="Exon type", loc='upper center', bbox_to_anchor=(0.7, 0.96), 
                      ncol=1, fancybox=True, shadow=False, borderaxespad=0., framealpha=1,
                      fontsize=fontsize_min, title_fontsize=fontsize_min)
            legend = ax.get_legend()
            for text in legend.get_texts():
                text.set_text(text.get_text().title())
    else:
        ax.get_legend().remove()
        
    ax.set_ylim(None,0)
    ax.set_xlim(0,18) #None

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))

    ax.set_xlabel("Sequence length (kbp)", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax.set_ylabel("Length-adjusted MFE (kcal/mol)", fontsize=fontsize_min, labelpad=axis_label_pad)
    
    ax.tick_params(axis='x', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    
    for collection in ax.collections:
        if isinstance(collection, mpl.collections.PathCollection):  # e.g., rasterize only specific types
            collection.set_rasterized(True)
        #collection.set_rasterized(True)
    
    ax.grid(True, linewidth=0.25)
    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
    
    if not ax:
        plt.show()


def relplot_1panel(data_df, tissue_type):
    plt.figure(figsize=(10,6))
    sns.scatterplot(data=data_df,
                    x=tissue_type,
                    y="seq_corrected_mfe",
                    hue="type")
    plt.xlabel(f"{tissue_type.capitalize()}")
    plt.ylabel("Length-adjusted MFE")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


def relplot_2panel(data_df, psi_col, psi_type, col_to_wrap, col_hue, kind="scatter"):
    ax=sns.relplot(data=data_df,
                   x=psi_col,
                   hue=col_hue,
                   col=col_to_wrap,
                   y="seq_corrected_mfe",
                   kind=kind,
                   style="hexevent_type",
                   alpha=0.5,
                   palette=palette)
    ax.set_xlabels(f"{psi_type.capitalize()} PSI across tissues \n(ASCOT database)", clear_inner=False)
    ax.set_ylabels("Length-adjusted MFE")
    #plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()


def relplot_multipanel(melt_df):
    ax=sns.relplot(data=melt_df,
                   y="seq_corrected_mfe",
                   x="value",
                   hue="inversion",
                   kind="scatter",
                   col="variable",
                   style="type",
                   alpha=0.5,
                   col_wrap=4,
                   palette=palette) #col_wrap
    ax.set_xlabels("PSI across tissues\n(ASCOT database)", clear_inner=False)
    ax.set_ylabels("Length-adjusted MFE")
    plt.show()


def relplot_2panel_2x(data_df, psi_type):
    """
    """
    if psi_type == "mean":
        psi_col = "tissue_psi_mean"
    elif psi_type == "max":
        psi_col = "tissue_psi_max"
    elif psi_type == "min":
        psi_col = "tissue_psi_min"
    elif psi_type == "median":
        psi_col = "tissue_psi_median"
    
    relplot_2panel(data_df, psi_col, psi_type, col_to_wrap="type", col_hue="inversion")
    relplot_2panel(data_df, psi_col, psi_type, col_to_wrap="inversion", col_hue="type")


def generate_venn_diagrams(df_skip, df_const, option=3):
    """
    *********DEPRECATED -- use find_shared_alu_alu_pairs() instead***********
    df_skip: mfe_skippable_exonskipdb
    df_const: mfe_constitutive_exonskipdb
    option == 1 or 2 (or 3 new for non-mfe dfs)
    """
    if option == 1:
        set1 = set(df_skip[df_skip.columns[0:5]]
                   .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        set2 = set(df_const[df_const.columns[0:5]]
                    .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        venn2([set1, set2], ('Skippable', 'Constitutive'))
        plt.title("Alu-Alu pairs")
        plt.show()
    elif option == 2:
        set1b = set(df_skip[df_skip.columns[0:3]]
                    .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        set2b = set(df_const[df_const.columns[0:3]]
                    .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        #list(set1b)[0]
        venn2([set1b, set2b], ('Skippable', 'Constitutive'))
        plt.title("Exons")
        plt.show()
    elif option == 3:
        #.columns[1:6]
        keep_cols = ["chr", "start", "stop", "alu1", "alu2"]
        set1b = set(df_skip[keep_cols] 
                    .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        set2b = set(df_const[keep_cols]
                    .apply(lambda x: '-'.join(x.astype(str)),axis=1).values)
        venn2([set1b, set2b], ('Skippable', 'Constitutive'))
        plt.title("Alu-Alu pairs")
        plt.show()


def find_shared_alu_alu_pairs(exons_df1: pd.DataFrame, exons_df2: pd.DataFrame,
                              label1: str, label2: str):
    """
    Code from option 3 of generate_venn_diagrams()
    """
    keep_cols = ["chr", "start", "stop", "alu1", "alu2"]
    set1b = set(exons_df1[keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set2b = set(exons_df2[keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    venn2([set1b, set2b], (label1, label2))
    plt.title("Alu-Alu pairs")
    plt.show()


def find_shared_alu_alu_pairs3(exons_df1: pd.DataFrame, exons_df2: pd.DataFrame, exons_df3: pd.DataFrame,
                              label1: str, label2: str, label3: str):
    """
    Code from option 3 of generate_venn_diagrams()
    """
    keep_cols = ["chr", "start", "stop", "alu1", "alu2"]
    set1b = set(exons_df1[keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set2b = set(exons_df2[keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set3b = set(exons_df3[keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    venn3([set1b, set2b, set3b], (label1, label2, label3))
    plt.title("Alu-Alu pairs")
    plt.show()


def find_shared_alu_alu_pairs_any(exons_df_list, label_list, colour_list, file_type, ax=None):
    """
    file_type: mfe or bedops
    
    Plot either Venn2 or Venn3 depending on size of lists
    # hexevent_df_dct["skipped"], hexevent_df_dct["constitutive"], esdb_df_dct["constitutive"]
    # bedops_skippable_pairs, bedops_constitutive_pairs, bedops_exonskipdb_skippable_pairs
    """
    fontsize_min = 5
    #plt.rcParams.update({'font.size': 20})
    plt.rcParams.update({'font.size': fontsize_min}) #20
    #plt.rcParams['figure.figsize'] = [20, 10]
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    if file_type == "mfe":
        keep_cols = ["chr", "start", "stop", "alu1", "alu2"]
    elif file_type == "bedops":
        keep_cols = ["exon_chr", "upstream_alu_start", "downstream_alu_end",
                     "upstream_alu_subfamily", "downstream_alu_subfamily"]
    if len(exons_df_list) == 2:
        set1 = set(exons_df_list[0][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df_list[1][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        v = venn2([set1, set2], (label_list[0], label_list[1]))
        v.get_patch_by_id('10').set_color(colour_list[0])
        v.get_patch_by_id('01').set_color(colour_list[1])
    elif len(exons_df_list) == 3:
        set1 = set(exons_df_list[0][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df_list[1][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set3 = set(exons_df_list[2][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        print(f"The union of sets {label_list[0]} and {label_list[1]} is of size: {len(set1.union(set2)):,}")
        print(f"Set {label_list[2]} is of size: {len(set3):,}")
        v = venn3([set1, set2, set3], set_labels=(label_list[0], label_list[1], label_list[2]), ax=ax)  #added ax=axes[0][0]
        v.get_patch_by_id('100').set_color(colour_list[0])
        v.get_patch_by_id('010').set_color(colour_list[1])
        if v.get_patch_by_id('001'):
            v.get_patch_by_id('001').set_color(colour_list[2])
        if v.get_patch_by_id('110'):
            v.get_patch_by_id('110').set_color(colour_list[0])
        for text in v.set_labels:
            text.set_fontsize(fontsize_min) #20
        # for x in range(len(v.subset_labels)):
        #     if v.subset_labels[x] is not None:
        #         v.subset_labels[x].set_fontsize(fontsize_min) #20

        for x in range(len(v.subset_labels)):
            if v.subset_labels[x] is not None:
                # Get the label text, apply formatting with comma separator, and set it back
                label_text = v.subset_labels[x].get_text()
                if label_text.isdigit():  # Check if the label is numeric
                    formatted_text = f"{int(label_text):,}"  # Add commas to thousands
                    v.subset_labels[x].set_text(formatted_text)
                # Set font size after formatting
                v.subset_labels[x].set_fontsize(fontsize_min)
    plt.subplots_adjust(bottom=0, top=1, left=0, right=1)      
    #plt.title("Alu-Alu pairs")
    #plt.show()


def find_shared_alu_alu_pairs_bedops(exons_df_list, label_list, colour_list):
    """
    Plot either Venn2 or Venn3 depending on size of lists
    # bedops_skippable_pairs, bedops_constitutive_pairs, bedops_exonskipdb_skippable_pairs
    """
    plt.rcParams['figure.figsize'] = [10, 5]
    keep_cols = ["exon_chr", "upstream_alu_start", "downstream_alu_end", 
                 "upstream_alu_subfamily", "downstream_alu_subfamily"]
    if len(exons_df_list) == 2:
        set1 = set(exons_df_list[0][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df_list[1][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        v = venn2([set1, set2], (label_list[0], label_list[1]))
        v.get_patch_by_id('10').set_color(colour_list[0])
        v.get_patch_by_id('01').set_color(colour_list[1])
    elif len(exons_df_list) == 3:
        set1 = set(exons_df_list[0][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df_list[1][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set3 = set(exons_df_list[2][keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        v = venn3([set1, set2, set3], (label_list[0], label_list[1], label_list[2]))
        v.get_patch_by_id('100').set_color(colour_list[0])
        v.get_patch_by_id('010').set_color(colour_list[1])
        print(colour_list[2])
        v.get_patch_by_id('001').set_color(colour_list[2])
    plt.title("Alu-Alu pairs")
    plt.show()


def find_shared_exons(exons_df1: pd.DataFrame, exons_df2: pd.DataFrame, label1: str, label2: str):
    """
    Code from option 3 of generate_venn_diagrams()
    """
    exon_cols = ["exon_chr", "exon_start", "exon_end", "exon_gene"]
    exon_keep_cols = ["exon_chr", "exon_start", "exon_end"]
    
    #bedops_skippable.groupby(by=exon_cols, group_keys=True).count().reset_index()[exon_cols]
    #bedops_exonskipdb_skippable.groupby(by=exon_cols, group_keys=True).count().reset_index()[exon_cols]

    #bedops_skippable_exons = (bedops_skippable
    #                          .groupby(by=exon_cols, group_keys=True)
    #                          .count()
    #                          .reset_index())
    #bedops_skippable_exons[["exon_start", "exon_end"]] = (bedops_skippable_exons[["exon_start", "exon_end"]]
    #                                                      .astype(int))
    exons_df1_groupby = groupby_exons(exons_df1)
    exons_df2_groupby = groupby_exons(exons_df2)

    set1b = set(exons_df1_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set2b = set(exons_df2_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    #print(set1b.intersection(set2b))
    
    venn2([set1b, set2b], (label1, label2))
    plt.title("Exons")
    plt.show()
    

def find_shared_exons3(exons_df1: pd.DataFrame, exons_df2: pd.DataFrame, exons_df3: pd.DataFrame,
                       label1: str, label2: str, label3: str):
    """
    Code from option 3 of generate_venn_diagrams()
    """
    exon_cols = ["exon_chr", "exon_start", "exon_end", "exon_gene"]
    exon_keep_cols = ["exon_chr", "exon_start", "exon_end"]
    
    exons_df1_groupby = groupby_exons(exons_df1)
    exons_df2_groupby = groupby_exons(exons_df2)
    exons_df3_groupby = groupby_exons(exons_df3)

    set1b = set(exons_df1_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set2b = set(exons_df2_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    set3b = set(exons_df3_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
    
    venn3([set1b, set2b, set3b], (label1, label2, label3))
    plt.title("Shared exons")
    plt.show()


def find_shared_exons_any(exons_df_list, label_list, colour_list, file_type, verbose=False, ax=None):
    """
    file_type: mfe or bedops
    
    Plot either Venn2 or Venn3 depending on size of lists
    # hexevent_df_dct["skipped"], hexevent_df_dct["constitutive"], esdb_df_dct["constitutive"]
    # bedops_skippable_pairs, bedops_constitutive_pairs, bedops_exonskipdb_skippable_pairs
    """
    fontsize_min = 5
    #plt.rcParams.update({'font.size': 20})
    plt.rcParams.update({'font.size': fontsize_min})
    #plt.rcParams['figure.figsize'] = [20, 10]
    plt.rcParams['font.family'] = 'DejaVu Sans'
    
    if file_type == "mfe":
        exon_keep_cols = ["exon_chr", "exon_start", "exon_end"] #XXX doesn't work bc need to intersection with original
    elif file_type == "bedops":
        exon_keep_cols = ["exon_chr", "exon_start", "exon_end"]
    if len(exons_df_list) == 2:
        exons_df1_groupby = groupby_exons(exons_df_list[0])
        exons_df2_groupby = groupby_exons(exons_df_list[1])
        set1 = set(exons_df1_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df2_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        v = venn2([set1, set2], (label_list[0], label_list[1]))
        v.get_patch_by_id('10').set_color(colour_list[0])
        v.get_patch_by_id('01').set_color(colour_list[1])
        if verbose:
            print(set1.intersection(set2))
    elif len(exons_df_list) == 3:
        exons_df1_groupby = groupby_exons(exons_df_list[0])
        exons_df2_groupby = groupby_exons(exons_df_list[1])
        exons_df3_groupby = groupby_exons(exons_df_list[2])
        set1 = set(exons_df1_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set2 = set(exons_df2_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        set3 = set(exons_df3_groupby[exon_keep_cols].apply(lambda x: '-'.join(x.astype(str)), axis=1).values)
        print(f"The union of sets {label_list[0]} and {label_list[1]} is of size: {len(set1.union(set2)):,}")
        print(f"Set {label_list[2]} is of size: {len(set3):,}")
        v = venn3([set1, set2, set3], (label_list[0], label_list[1], label_list[2]), ax=ax) #added ax=axes[0][0]
        v.get_patch_by_id('100').set_color(colour_list[0])
        v.get_patch_by_id('010').set_color(colour_list[1])
        
        for text in v.set_labels:
            text.set_fontsize(fontsize_min)
        
        # for x in range(len(v.subset_labels)):
        #     if v.subset_labels[x] is not None:
        #         v.subset_labels[x].set_fontsize(fontsize_min)

        for x in range(len(v.subset_labels)):
            if v.subset_labels[x] is not None:
                # Get the label text, apply formatting with comma separator, and set it back
                label_text = v.subset_labels[x].get_text()
                if label_text.isdigit():  # Check if the label is numeric
                    formatted_text = f"{int(label_text):,}"  # Add commas to thousands
                    v.subset_labels[x].set_text(formatted_text)
                # Set font size after formatting
                v.subset_labels[x].set_fontsize(fontsize_min)
        
        if v.get_patch_by_id('001'):
            v.get_patch_by_id('001').set_color(colour_list[2])
        if v.get_patch_by_id('110'):
            v.get_patch_by_id('110').set_color(colour_list[0])
    plt.subplots_adjust(bottom=0, top=1, left=0, right=1)        
    #plt.title("Exons")
    #plt.show()


def remove_shared_exons(bedops_df, exon_set, exon_type, exon_db, verbose=False):
    """
    bedops_df: bedops_skippable, bedops_exonskipdb_skippable, or bedops_constitutive
    exon_set: skippable_set, or constitutive_set
    """
    new_df = bedops_df.copy()
    new_df = new_df.dropna(subset=["exon_start", "exon_end"], how='any') # drop row if start/end undefined
    new_df[["exon_start", "exon_end"]] = new_df[["exon_start", "exon_end"]].astype(int)
    new_df["exon_key"] = (new_df[["exon_chr", "exon_start", "exon_end"]].apply(
        lambda x: '-'.join(x.astype(str)), axis=1).values)
    new_df = new_df[new_df['exon_key'].isin(list(exon_set))]
    new_df["exon_type"] = exon_type
    new_df["exon_db"] = exon_db
    if verbose:
        print(f"\nExon database: {exon_db}\nExon type: {exon_type}")
        print(f"Original number of entries: {len(bedops_df):,}")
        print(f"New number of entries: {len(new_df):,}")
    return new_df


def plot_redi_rna_editing_burden(df, ax, show_table=True):
    tick_params_width = 0.25
    tick_params_length = 1
    tick_params_pad = 2
    axis_label_pad = 2

    fontsize_min = 5
    plt.rcParams['font.family'] = 'DejaVu Sans'

    exon_def_mapping = {".": "Other", "WS": "WS", "SS":"SS", "WW":"WW", "SW":"SW"}
    category_order = ["WS", "SS", "SW", "WW", "Other"]
    custom_palette = ["#2ca02c", "#9467bd"]
    
    filtered_redi = df
    filtered_redi['MES_category_map'] = filtered_redi['MES_category'].map(exon_def_mapping)

    if show_table:
        ttest_upstream_redi_results = []
        grouped_data = filtered_redi.groupby("MES_category_map")

        for name, group in grouped_data:
            inv1 = group[group["psi_category"] == "PSI < 50%"]["rediportal_count"]
            inv2 = group[group["psi_category"] == "PSI > 50%"]["rediportal_count"]

            t_stat, p_value = ttest_ind(inv1, inv2)
            ttest_upstream_redi_results.append({"MES_category_map": name, "t_stat": t_stat, "p_value": p_value, "p_adj": p_value*5})

        ttest_upstream_redi_results_df = pd.DataFrame(ttest_upstream_redi_results)
        display(ttest_upstream_redi_results_df)
    
    sns.set_style("ticks")
    sns.violinplot(
        data=filtered_redi, x='MES_category_map', hue="psi_category", y="rediportal_count", 
        inner="box", order=category_order, palette=custom_palette, linewidth=0.5, ax=ax
    )

    ax.set_ylabel("RNA editing burden", fontsize=fontsize_min, labelpad=axis_label_pad)
    ax.set_xlabel("Exon definition category", fontsize=fontsize_min, labelpad=axis_label_pad)

    ax.tick_params(axis='x', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    ax.tick_params(axis='y', labelsize=fontsize_min, 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)

    ax.set_ylim(-99, 500)
    ax.tick_params(axis='both', which='both', 
                   length=tick_params_length, width=tick_params_width, pad=tick_params_pad)
    
    ax.yaxis.grid(True, linewidth=0.25)
    ax.xaxis.grid(False)

    ax.legend(title='', fancybox=True, shadow=False, borderaxespad=0., fontsize=fontsize_min)

    for spine in ax.spines.values():
        spine.set_linewidth(0.25)
