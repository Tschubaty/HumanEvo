import pandas as pd
from matplotlib import cm

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.preprocessing import Imputer
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from collections import OrderedDict

import seaborn
from matplotlib.colors import ListedColormap

NAN = "nan"
names = ['I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775',
         'I4315', 'I1053', 'I7421', 'I3255', 'I5835', 'I2514', 'I5748', 'I1633', 'I5950', 'I5838', 'I5742', 'I5743',
         'I2105', 'I2520', 'I2935', 'I2978', 'I2980', 'I3133', 'I4634', 'I1965', 'I1631', 'I1632', 'I1961', 'LBK',
         'I1496', 'I5077', 'I1962', 'I4438', 'I15', 'I2134', 'I1507', 'I4878', 'I4432', 'I4873', 'Los', 'LaB', 'I4596',
         'I5233', 'I0708', 'I1960', 'I4914', 'I4875', 'I4877', 'I2139', 'SF1', 'I1734', 'I5236', 'I5235', 'M45', 'Ust']
time_vec = np.asarray(
    [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7,
     4.9, 5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7, 7.3, 7.3, 7.5, 7.6, 7.7, 7.8, 7.9, 7.9, 8, 8, 8, 8,
     8.1, 8.1, 8.1, 8.4, 8.5, 8.9, 9, 9.2, 10, 10.8, 12, 45])

cmap = colors


def legend_generator():
    """
    generate the legend of the final graph (concatenate the name with the years
    :return: array of concatenated names + time_vec elements
    """
    arr = []
    for i in range(len(names)):
        arr.append(names[i] + " - " + str(time_vec[i]) + "K years")
    return arr


def main():
    """
    This is the main function, reads the csv and does principle component analysis on it
    :return:
    """
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
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_22_united.csv"
    ]
    # origin_arr = ["Z://Matlab_Data//Matlab_data_28.2.19//working_example//chr_example_united.csv"]

    output_arr = []
    for i in range(len(origin_arr)):
        output_arr.append(
            "Z://PCA//outputs//chr_" + str(
                i + 1) + "_PCA.png")

    for i in range(len(origin_arr)):
        print("chr " + str(i))
        PCA_x(origin_arr[i], output_arr[i], i+1)


def PCA_x(origin, output, chr):
    """
    Perform pca on the csv file in the origin and output a graph into the file
    :param origin: the origin of the csv file
    :param output: address to output the file
    :param chr: chromosome being examined
    :return: nothing
    """
    df = pd.read_csv(origin, names=names)


    df = handle_NaNs(df)

    x = df.values


    pca = PCA(n_components=2)

    principalComponents = pca.fit_transform(np.transpose(x))
    exp_variance = pca.explained_variance_ratio_
    print("explained var", exp_variance)

    principalDf = pd.DataFrame(data=principalComponents
                               , columns=['principal component 1', 'principal component 2'])


    headers = pd.DataFrame({'humans': names})

    finalDf = pd.concat([principalDf, headers], axis=1)

    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(1, 1, 1)

    # plt.ylim(-1, 1)
    # plt.xlim(-1, 1)
    ax.set_xlabel('Dim1 ('+ str(exp_variance[0]*100)[:8]+'%'+')', fontsize=15)
    ax.set_ylabel('Dim2 ('+str(exp_variance[1]*100)[:8]+'%'+')', fontsize=15)
    ax.set_title('Two Component PCA Chr ' + str(chr), fontsize=20)
    targets = legend_generator()
    x = cm.get_cmap('tab20b', 8)
    for i, target in enumerate(targets):
        plt.scatter(principalDf['principal component 1'].iloc[i]
                    , principalDf['principal component 2'].iloc[i]
                    , s=30, c=[x(i / len(names))])

        text = ax.annotate(names[i], (principalDf.iloc[i, 0], principalDf.iloc[i, 1]))
        text.set_fontsize(9)

    ax.legend(targets, loc=2, prop={'size': 6}, )
    ax.grid()
    plt.savefig(output)
    plt.clf()


def handle_NaNs(df):
    """
    handle the NANS in the data, fill them with zero
    :param df: dataframe
    :return: df with NANS filled with mean
    """
    temp = df.apply(lambda x: x.fillna(x.mean()), axis=1)  ## possible - to switch to fill with mean value
    return temp.dropna(axis=0)


main()
