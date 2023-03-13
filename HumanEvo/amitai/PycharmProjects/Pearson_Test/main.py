from pydoc import help
from scipy.stats.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def main(origin, output):
    """
    Read a csv and create a plot of p values over threshold using the pearson correlation test
    :param origin: origin of data
    :param output: where to output data
    :return: nothing
    """
    # help(pearsonr)
    # origin = "Z://Matlab_Data//training_data//chr_1_training.csv"
    # df = df[:10000]

    df = pd.read_csv(origin, names=['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520'])
    df = handle_NaNs(df)
    df = best_Var(df)
    # print(df)
    meth_mat = df.values
    time_vec = [8.1, 3.8, 1, 7, 7.7, 6.1, 5.1] # Time in thousands of years
    q = 0.05
    unsorted_p_vals = do_pearson(meth_mat, time_vec)
    sorted_p_vals = sorted(unsorted_p_vals, key=lambda tup: tup[0], reverse=False)
    # print("sorted p vals")
    print(sorted_p_vals)
    N = len(sorted_p_vals)
    i = np.arange(1, N + 1)
    p_vals_only = [elem[0] for elem in sorted_p_vals]
    rejected, pval_corr_arr = fdrcorrection(p_vals_only, q, method='i')

    plt.title("FDR - Pearson Correlation")
    plt.plot(i, p_vals_only, 'b.', label='$p(i)$')
    plt.xlabel('Sorted tests ')
    plt.ylabel('p value')
    plt.plot(i, q * i / N, 'r', label='$q i / N $')
    plt.legend()
    plt.savefig(output)
    plt.clf()
    # print("pval corrected arrs\n", pval_corr_arr)
    #
    # adjusted_qs = calc_q_adjusted(pval_corr_arr)
    #
    # print("adjusted_qs")
    # print(adjusted_qs)


def best_Var(df):
    """
    Reindex the data by the best variance
    :param df: dataframe
    :return: newly ordered dataframe
    """
    temp = df.reindex(df.var(axis=1).sort_values().index, axis=0)[::-1]
    return temp


def calc_q_adjusted(q_corr):
    """
    calculate the q-adujusted value (min([n * p / i), q_corr of all next indices in list]
    :param q_corr: list of q_corr  values
    :return: q_adjusted list
    """

    q_adjusted = []

    for i in range(len(q_corr)):
        q_adjusted.append(min(q_corr[i:]))

    return q_adjusted


def handle_NaNs(df):
    """
    handle the NANS in the data, fill them with mean by row
    :param df: dataframe
    :return: df with NANS filled with mean
    """
    temp = df.apply(lambda x: x.fillna(x.mean(skipna=True)), axis=1)
    return temp.dropna(axis=0)


def do_pearson(meth_math, time_vec):
    """
    do a pearson test for every row in the data
    :param meth_math: methylation matrix
    :param time_vec: time vector labels
    :return: pearson p values for every row
    """
    p_vals_arr = []
    for i, row in enumerate(meth_math):
        #     print(row)
        p_vals_arr.append((pearsonr(row, time_vec)[1], i))
    print("pearson p vals")
    print(p_vals_arr)
    return p_vals_arr


if __name__ == '__main__':
    """
    do the pearson correlation test for every chrosmosome in the data
    """
    origin_arr = ["Z://Matlab_Data//chr_1_united.csv", "Z://Matlab_Data//chr_2_united.csv",
                  "Z://Matlab_Data//chr_3_united.csv", "Z://Matlab_Data//chr_4_united.csv",
                  "Z://Matlab_Data//chr_5_united.csv", "Z://Matlab_Data//chr_6_united.csv",
                  "Z://Matlab_Data//chr_7_united.csv", "Z://Matlab_Data//chr_8_united.csv",
                  "Z://Matlab_Data//chr_9_united.csv", "Z://Matlab_Data//chr_10_united.csv",
                  "Z://Matlab_Data//chr_11_united.csv", "Z://Matlab_Data//chr_12_united.csv",
                  "Z://Matlab_Data//chr_13_united.csv", "Z://Matlab_Data//chr_14_united.csv",
                  "Z://Matlab_Data//chr_15_united.csv","Z://Matlab_Data//chr_16_united.csv",
                  "Z://Matlab_Data//chr_17_united.csv", "Z://Matlab_Data//chr_18_united.csv",
                  "Z://Matlab_Data//chr_19_united.csv", "Z://Matlab_Data//chr_20_united.csv",
                  "Z://Matlab_Data//chr_21_united.csv", "Z://Matlab_Data//chr_22_united.csv"]
    output_arr = []
    for i in range(22):
        output_arr.append("Z://FDR//Graphs//chr_"+ str(i+1) +  "_FDR_all.png" )



    for i in range(len(origin_arr)):
        main(origin_arr[i], output_arr[i])

    print("done")
