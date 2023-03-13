import matplotlib.pyplot as plt
from sklearn.linear_model import SGDClassifier
import numpy as np
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
import pandas as pd
import os
import math

TEST_NUM = 2


def main(origin, output):
    # origin = "Z://Matlab_Data//training_data//chr_1_training.csv"

    names = ['p_0708', 'p_1053', 'p_1116', 'p_1496', 'p_1507', 'p_1631', 'p_2520']
    df = pd.read_csv(origin, names=names)
    df = handle_NaNs(df)
    time_vec = [8.1,3.8, 1, 7, 7.7, 6.1, 5.1] # Time in thousands of years

    # increase by 10 to create integer values
    for elem in time_vec:
        elem*=10
    df = best_Var(df)
    data = df.values

    # Split the data into train and test

    X_train = np.transpose(np.matrix(data[:, :-TEST_NUM]))
    Y_train = np.transpose(np.matrix(time_vec[:-TEST_NUM]))
    X_test = np.transpose(np.matrix(data[:, -TEST_NUM:]))
    Y_test = np.transpose(np.matrix(time_vec[-TEST_NUM:]))




    # # Perform linear regression prediction
    regr = linear_model.LinearRegression()
    print(X_train)
    print(Y_train)
    regr.fit(X_train, Y_train)
    weights = np.transpose(regr.coef_)

    # Linear regression Coefficients:




    # gradient_descent on training
    clf = SGDClassifier(loss="hinge", alpha=0.01, max_iter=200, fit_intercept=True)
    Y_temp = np.array ( Y_train, dtype=int)
    weights = clf.fit(X_train, Y_temp)


    prediction = clf.predict(X_test)

    # calculate distance using MSE
    print("MSE distance = ", math.sqrt(mean_squared_error(prediction, Y_test)))

    print("X test", X_test)
    Y_test = np.array(Y_test.T)
    print("Y test", type(Y_test))
    print("prediction", type(prediction))
    # plt.scatter(X_test, Y_test, color='black')
    x_axis = [i for i in range(10)]
    plt.title("Linear Regression Classification \n" + os.path.basename(origin)[:-4])
    plt.ylabel("Thousand Years")
    plt.xlabel("Name of Sample")

    plt.scatter(names[-TEST_NUM:], Y_test, color='red', linewidth=2, label="Actual Years")

    plt.plot(names[-TEST_NUM:], prediction, color='blue', linewidth=1, label="Years prediction")
    plt.legend()
    plt.savefig(output)
    plt.clf()


def handle_NaNs(df):
    """
    handle the NANS in the data, fill them with zero
    :param df: dataframe
    :return: df with NANS filled with mean
    """
    temp = df.apply(lambda x: x.fillna(0), axis=1)
    return temp.dropna(axis=0)


def best_Var(df):
    """
    Reindex the data by the best variance
    :param df: dataframe
    :return: newly ordered dataframe
    """
    temp = df.reindex(df.var(axis=1).sort_values().index, axis=0)[::-1]
    return temp


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
                 "Z://Matlab_Data//chr_15_united.csv", "Z://Matlab_Data//chr_16_united.csv",
                 "Z://Matlab_Data//chr_17_united.csv", "Z://Matlab_Data//chr_18_united.csv",
                 "Z://Matlab_Data//chr_19_united.csv", "Z://Matlab_Data//chr_20_united.csv",
                 "Z://Matlab_Data//chr_21_united.csv", "Z://Matlab_Data//chr_22_united.csv"]

    output_arr = []
    for i in range(len(origin_arr)):
       output_arr.append("Z://Lin_reg//Graphs_grad_desc//chr_" + str(i + 1) + "_SGD_training_" + str(TEST_NUM )+ ".png")

    for i in range(len(origin_arr)):
       main(origin_arr[i], output_arr[i])
# main("", "Z://Examples//Lin_reg_examples//Lin_reg_and_grad_descent.png")
print("done")

