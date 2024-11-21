#!/usr/bin/env python
import os
import pandas as pd
from pandas import Series
from numpy import var, mean, sqrt
from statsmodels.stats.weightstats import ztest
from scipy.stats.contingency import odds_ratio
from scipy import stats


def cohend(d1: Series, d2: Series) -> float:
    """Copied from: https://stackoverflow.com/questions/21532471"""
    # calculate the size of samples
    n1, n2 = len(d1), len(d2)
    # calculate the variance of the samples
    s1, s2 = var(d1, ddof=1), var(d2, ddof=1)
    # calculate the pooled standard deviation
    s = sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
    # calculate the means of the samples
    u1, u2 = mean(d1), mean(d2)
    # return the effect size
    return (u1 - u2) / s


def run_stats_tests_2groups(df_skip, df_const):
    """
    df_skip: skipped
    df_const: constitutive
    """
    # Test hypothesis that df_skip mean is less than df_const mean (alt='less')
    t_result = stats.ttest_ind(df_skip.seq_corrected_mfe, df_const.seq_corrected_mfe,
                               equal_var=False, alternative='less')
    z_result = ztest(df_skip.seq_corrected_mfe, df_const.seq_corrected_mfe, alternative="smaller")
    
    # Effect size
    d_result = cohend(d1=df_skip.seq_corrected_mfe, d2=df_const.seq_corrected_mfe)
    
    #print(f"t-test: {t_result}")
    #print(f"z-test: {z_result}")
    #print(f"cohen's d: {d_result}")
    return {"t": t_result, "z": z_result, "d": d_result}


def perform_mfe_stats_tests(df1, df2):
    """
    Generally, df1 = skipped, df2 = constitutive
    1-tailed, alternative = (less or smaller) for t-test and z-test
    """
    stats_dct = {}
    
    # Test hypothesis that df_skip mean is less than df_const mean (alt='less')
    stats_dct["ttest"] = stats.ttest_ind(df1.seq_corrected_mfe, df2.seq_corrected_mfe,
                                         equal_var=False, alternative='less')
    stats_dct["ztest"] = ztest(df1.seq_corrected_mfe, df2.seq_corrected_mfe, alternative="smaller")
    
    # Effect size
    stats_dct["cohend"] = cohend(d1=df1.seq_corrected_mfe, d2=df2.seq_corrected_mfe)
    return stats_dct


def calculate_odds_ratio(row):
    """
    This function assumes that there is a column called "table" which contains
    2x2 contingency tables.
    Then, it calculates OR on a per-row basis, and returns the OR and CI boundaries.
    """
    # Extract the contingency table from the row
    table = row['table']
    print(table)
    # Calculate odds ratio using scipy.stats.contingency.odds_ratio
    odds_ratio_res = odds_ratio(table)
    odds_ratio_val = odds_ratio_res.statistic
    print(odds_ratio_val)
    odds_ratio_ci = odds_ratio_res.confidence_interval(confidence_level=0.95)
    print(odds_ratio_ci.low, odds_ratio_ci.high)
    
    return odds_ratio_val, odds_ratio_ci.low, odds_ratio_ci.high