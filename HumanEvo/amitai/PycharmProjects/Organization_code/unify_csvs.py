
import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.datasets import load_digits
from scipy import stats


names = ['I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775',
         'I4315', 'I1053', 'I7421', 'I3255', 'I5835', 'I2514', 'I5748', 'I1633', 'I5950', 'I5838', 'I5742', 'I5743',
         'I2105', 'I2520', 'I2935', 'I2978', 'I2980', 'I3133', 'I4634', 'I1965', 'I1631', 'I1632', 'I1961', 'LBK',
         'I1496', 'I5077', 'I1962', 'I4438', 'I15', 'I2134', 'I1507', 'I4878', 'I4432', 'I4873', 'Los', 'LaB', 'I4596',
         'I5233', 'I0708', 'I1960', 'I4914', 'I4875', 'I4877', 'I2139', 'SF1', 'I1734', 'I5236', 'I5235', 'M45', 'Ust']
time_vec = np.asarray(
    [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7,
     4.9, 5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7, 7.3, 7.3, 7.5, 7.6, 7.7, 7.8, 7.9, 7.9, 8, 8, 8, 8,
     8.1, 8.1, 8.1, 8.4, 8.5, 8.9, 9, 9.2, 10, 10.8, 12, 45])




def main():
    origin_arr = [
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_1_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_2_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_3_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_4_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_5_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_6_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_7_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_8_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_9_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_10_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_11_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_12_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_13_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_14_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_15_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_16_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_17_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_18_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_19_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_20_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_21_united.csv",
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_22_united.csv",
    ]

    output_add ="Z://Matlab_Data//Matlab_data_28.2.19//Unified_data//meth_unified.csv"
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




