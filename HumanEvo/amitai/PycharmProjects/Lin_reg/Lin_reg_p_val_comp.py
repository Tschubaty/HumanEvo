import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.datasets import load_digits
from scipy import stats

resolution = 10000
HIGHEST_VAR_NUM = 100000

names = ['I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775',
         'I4315', 'I1053', 'I7421', 'I3255', 'I5835', 'I2514', 'I5748', 'I1633', 'I5950', 'I5838', 'I5742', 'I5743',
         'I2105', 'I2520', 'I2935', 'I2978', 'I2980', 'I3133', 'I4634', 'I1965', 'I1631', 'I1632', 'I1961', 'LBK',
         'I1496', 'I5077', 'I1962', 'I4438', 'I15', 'I2134', 'I1507', 'I4878', 'I4432', 'I4873', 'Los', 'LaB', 'I4596',
         'I5233', 'I0708', 'I1960', 'I4914', 'I4875', 'I4877', 'I2139', 'SF1', 'I1734', 'I5236', 'I5235', 'M45', 'Ust']
time_vec = np.asarray(
    [1, 1, 1.2, 1.6, 1.8, 2.2, 2.5, 2.8, 3, 3.1, 3.3, 3.4, 3.5, 3.8, 3.8, 4, 4.2, 4.3, 4.4, 4.5, 4.5, 4.6, 4.7, 4.7,
     4.9, 5.1, 5.1, 5.1, 5.2, 5.4, 5.4, 5.8, 6.1, 6.1, 6.1, 7, 7, 7, 7.3, 7.3, 7.5, 7.6, 7.7, 7.8, 7.9, 7.9, 8, 8, 8, 8,
     8.1, 8.1, 8.1, 8.4, 8.5, 8.9, 9, 9.2, 10, 10.8, 12, 45])


class LinearRegression(linear_model.LinearRegression):
    """
    LinearRegression class after sklearn's, but calculate t-statistics
    and p-values for model coefficients (betas).
    Additional attributes available after .fit()
    are `t` and `p` which are of the shape (y.shape[1], X.shape[1])
    which is (n_features, n_coefs)
    This class sets the intercept to 0 by default, since usually we include it
    in X.
    """

    def init(self, *args, **kwargs):
        if not "fit_intercept" in kwargs:
            kwargs['fit_intercept'] = False
        super(LinearRegression, self) \
            .init(*args, **kwargs)

    def fit(self, X, y, n_jobs=1):
        self = super(LinearRegression, self).fit(X, y, n_jobs)

        sse = np.sum((self.predict(X) - y) ** 2, axis=0) / float(X.shape[0] - X.shape[1])
        se = np.array([
            np.sqrt(np.diagonal(sse[i] * np.linalg.inv(np.dot(X.T, X))))
            for i in range(sse.shape[0])
        ])

        self.t = self.coef_ / se
        self.p = 2 * (1 - stats.t.cdf(np.abs(self.t), y.shape[0] - X.shape[1]))
        return self


def createarr():
    """
    This method is a random array generator with numbers between 0,1,
    immitating p values with significant findings (size between 1-2 million)
    :return:
    """
    arr = []
    for i in range(1000000):
        arr.append(random.random())
    for j in range(300000):
        arr.append(random.uniform(0.8, 1))
    for i in range(100000):
        arr.append(random.uniform(0.2, 0.3))
    for i in range(300000):
        arr.append(random.random())
    return arr


def visualize(signif_vals_arr, N, output):
    """
    This is a method that visualizes p values array with N as the number
    of tessts run
    :param signif_vals_arr: array of values between 0,2
    :param N: number of tests
    :param output: output address
    :return: output a graph where the x axis is the mean of the signif vals
    arr per position (since it is a large size we reduce it to a given resolution)
    into the address at output
    """
    print(N)
    x_ax = np.arange(0, N, resolution)
    y_ax = np.zeros((int(N / resolution))+1)
    non_zero_vals_arr = []
    j = 0
    for i in range(0, N, resolution):
        l = signif_vals_arr[i: i + resolution]
        non_zero_vals = np.count_nonzero(l)
        non_zero_vals_arr.append(non_zero_vals)
        if (non_zero_vals == 0):
            x = 'nan'
        else:
            x = l[l.nonzero()].mean()
        y_ax[j] = x
        j += 1
    print(y_ax)
    plt.title("Linear Regression")
    plt.xlabel("Methylation position")
    plt.ylabel("p value")
    plt.scatter(x_ax, y_ax)
    plt.savefig(output)
    plt.clf()


    plt.title("Non Zero Vals")
    plt.xlabel("Methylation position")
    plt.ylabel("number of highest variance sites")
    plt.scatter(x_ax, non_zero_vals_arr)
    plt.savefig("Z://Lin_reg//Graphs_lin_reg_4.3//chr_" + output[-13] + "_num_highest_variance.png")
    plt.clf()



def lin_reg(origin, chr_num):
    """
    Perform simple linear regression
    :param origin: address of csv file to be read
    :param chr_num: the number of chromosome being examined
    :return: array of p_values per position
    """
    df = pd.read_csv(origin, names=names)
    orig_size = df.shape[0]
    print("size of chr ",  orig_size)
    df = handle_NaNs(df)
    print("NAN's handles")

    df = best_Var(df)
    print("best variance taken")

    df = df.reset_index().set_index('index', drop=False)
    print(df)
    data = df.values

    # Split the data into train and test

    X_train = np.asarray(np.matrix(data))
    Y_train = time_vec

    # # basic linear regression using sklearn

    p_value_arr = []

    y = Y_train.reshape(-1, 1)
    for i, meth_position in enumerate(X_train):
        x = meth_position.reshape(-1, 1)

        # this is the actual linear regression
        if (np.count_nonzero(x[1:]) == 0):
            pass
        else:
            regr = LinearRegression().fit(x[1:], y)
            intercept, r_value, p_value, slope, chr_index= regr.intercept_, regr.score(x[1:], y) \
                , regr.p, regr.coef_, x[0]
            # print("intercept: %f    p value: %f  " % (intercept, p_value))
            #
            # print("r-squared: %f" % r_value ** 2)

            # insert values into the array as tuples in the following order:
            # (chr num, index of meth position, p, r**2, )
            p_value_arr.append((chr_num, int(x[0][0]), p_value[0][0], r_value ** 2))

            # # plot the p value
            # plt.title(p_value[0][0])
            # plt.plot(x, y, 'o', label='original data')
            # plt.plot(x, intercept + slope * x, 'r', label='fitted line')
            # plt.legend()
            # plt.show()
    return p_value_arr, orig_size


def handle_NaNs(df):
    """
    handle the NANS in the data, fill them with zero
    :param df: dataframe
    :return: df with NANS filled with mean
    """
    print(df)
    temp = df.apply(lambda x: x.fillna(x.mean()), axis=1)  ## possible - to switch to fill with mean value
    return temp.dropna(axis=0)


def best_Var(df):
    """
    Reindex the data by the best variance
    :param df: dataframe
    :return: newly ordered dataframe
    """
    temp = df.reindex(df.var(axis=1).sort_values().index, axis=0)[::-1]
    return temp[:HIGHEST_VAR_NUM]


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
        "Z://Matlab_Data//Matlab_data_28.2.19//chr_22_united.csv"
    ]

    output_arr = []
    for i in range(len(origin_arr)):
        output_arr.append(
            "Z://Lin_reg//Graphs_lin_reg_4.3//chr_" + str(
                i + 1) + "_lin_reg.png")
    for i in range(len(origin_arr)):
        unsorted_p_vals , N= lin_reg(origin_arr[i], i + 1)

        # for FDR
        # sorted_p_vals = sorted(unsorted_p_vals, key=lambda tup: tup[2], reverse=False)
        # print(sorted_p_vals)
        # i_arr = np.arange(1, N + 1)
        # thresh = [(i * SIGNIFICANT_THRESH / N) for i in i_arr]
        # rev_thresh = thresh[::-1]
        # signif_vals_arr = []
        # for j, single_test in enumerate(sorted_p_vals[::-1]):
        #     if (single_test[2] < rev_thresh[i]):
        #         print("found significant", single_test[2], rev_thresh[i])
        #
        #         signif_vals_arr = sorted_p_vals[:N - i]
        #         break
        # print("significant vals arr", signif_vals_arr)
        # visualize(sorted_p_vals, len(signif_vals_arr) , output_arr[i])
        print("Size = ", N)
        p_vals = np.zeros(N)

        for elem in unsorted_p_vals:

            p_vals[elem[1]] = elem[2]

        print("p valas arr size", p_vals.shape)
        print("p vals arr, ", p_vals)
        visualize(p_vals, N, output_arr[i])


main()
