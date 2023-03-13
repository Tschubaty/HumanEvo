import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.datasets import load_digits
from scipy import stats


names = ['I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775',
         'I4315', 'I1053', 'I7421', 'I3255', 'I5835', 'I2514', 'I5748', 'I1633', 'I5950', 'I5838', 'I5742', 'I5743',
          'I2520', 'I2935', 'I2978', 'I2980', 'I3133', 'I4634', 'I1965', 'I1631', 'I1632', 'I1961',
         'I1496', 'I5077', 'I1962', 'I4438', 'I2134', 'I1507', 'I4878', 'I4432', 'I4873',  'I4596',
         'I5233', 'I0708', 'I1960', 'I4914', 'I4875', 'I4877', 'I2139',  'I1734', 'I5236', 'I5235']
time_vec = np.asarray(
    [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7,
      5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7.3, 7.3,  7.6, 7.7, 7.8, 7.9, 7.9, 8, 8,
     8.1, 8.1, 8.1, 8.4, 8.5, 8.9,  9.2, 10, 10.8])




def main():
    origin_arr = [
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_1_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_2_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_3_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_4_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_5_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_6_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_7_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_8_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_9_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_10_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_11_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_12_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_13_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_14_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_15_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_16_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_17_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_18_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_19_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_20_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_21_simu_united.csv",
        # "X://Matlab_Data//Matlab_data_28.2.19//chr_22_simu_united.csv",



        "chr_1_simu_united.csv",
        "chr_2_simu_united.csv",
        "chr_3_simu_united.csv",
        "chr_4_simu_united.csv",
        "chr_5_simu_united.csv",
        "chr_6_simu_united.csv",
        "chr_7_simu_united.csv",
        "chr_8_simu_united.csv",
        "chr_9_simu_united.csv",
        "chr_10_simu_united.csv",
        "chr_11_simu_united.csv",
        "chr_12_simu_united.csv",
        "chr_13_simu_united.csv",
        "chr_14_simu_united.csv",
        "chr_15_simu_united.csv",
        "chr_16_simu_united.csv",
        "chr_17_simu_united.csv",
        "chr_18_simu_united.csv",
        "chr_19_simu_united.csv",
        "chr_20_simu_united.csv",
        "chr_21_simu_united.csv",
        "chr_22_simu_united.csv",

    ]

    # output_add ="X://Matlab_Data//Matlab_data_28.2.19//simu//united//meth_unified_simu.csv"
    output_add ="meth_unified_simu.csv"
    df_arr = []
    print("here")

    for i in range(len(origin_arr)):
        df_arr.append(pd.read_csv(origin_arr[i], names=names))
        new_col = [i+1]*df_arr[i].shape[0]
        df_arr[i].insert(loc=0, column='chr', value=new_col)


    print(len(df_arr))
    df = pd.concat(df_arr)
    df.to_csv(output_add, sep=',', index=False, header=None)


main()




