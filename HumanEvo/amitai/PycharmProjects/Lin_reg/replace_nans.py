
import pandas as pd

def handle_NaNs(origin, output):
    """
    handle the NANS in the data, fill them with zero
    :param df: dataframe
    :return: df with NANS filled with mean
    """
    names = ['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520']
    df = pd.read_csv(origin, names=names)
    print(df)
    temp = df.apply(lambda x: x.fillna(x.mean()), axis=1)  ## possible - to switch to fill with mean value
    temp = df.dropna(axis=0)


    temp.to_csv(output)
    print("done")




def main():
    """
    Execute the program
    """
    origin_arr = ["Z://Matlab_Data//training_data//meth_window_example_1.csv"]
    # origin_arr = ["Z://Matlab_Data//Matlab_data_january//chr_1_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_2_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_3_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_4_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_5_united.csv"
    #               ,"Z://Matlab_Data//Matlab_data_january//chr_5_united.csv" , "Z://Matlab_Data//Matlab_data_january//chr_6_united.csv","Z://Matlab_Data//Matlab_data_january//chr_7_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_8_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_9_united.csv"
    #               , "Z://Matlab_Data//Matlab_data_january//chr_10_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_11_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_12_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_13_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_14_united.csv",
    #               "Z://Matlab_Data//Matlab_data_january//chr_15_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_16_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_17_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_18_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_19_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_20_united.csv"
    #               , "Z://Matlab_Data//Matlab_data_january//chr_21_united.csv", "Z://Matlab_Data//Matlab_data_january//chr_22_united.csv"]


    output_arr = []
    for i in range(len(origin_arr)):
        output_arr.append(
            "Z://Matlab_Data//Matlab_data_january//nans_removed//chr_" + str(
                i + 1) + "_nans_removed.csv")

    for i in range(len(origin_arr)):
        handle_NaNs(origin_arr[i], output_arr[i])


main()