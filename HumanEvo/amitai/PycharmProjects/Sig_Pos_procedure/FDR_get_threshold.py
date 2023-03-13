"""
This file receives two outputs of the procedure for the simulated data and the regular data
and returns the threshold for FDR<= 0.05 for the data
"""


import pandas as pd
import sys

real_results ="interval_procedure_output.csv"
simu_results = "interval_procedure_output_simu.csv"


col_names = ['chromosome',  'position', 'pvalue', '-log 10 p value', 'pvalue', '-log 10 p value', 'min_log_10_p_val']
ALPHA = sys.argv[1]


def main():
    df1 = pd.read_csv(real_results, names=col_names, dtype=float)
    df2 = pd.read_csv(simu_results, names=col_names, dtype=float)
    num_rows_1 = df1.shape[0]
    thresh = 0
    for (i, row1)in df1[::-1].iterrows():
        significant_sim = df2['min_log_10_p_val'] > row1['min_log_10_p_val']
        a = df2[significant_sim].shape[0]
        b = num_rows_1 - i
        if((a/b) < ALPHA):
            thresh = row1['min_log_10_p_val']
            break

    return thresh





if __name__ == '__main__':
    thresh = main()
    print("The threshold for FDR = ", ALPHA, " is ", thresh)