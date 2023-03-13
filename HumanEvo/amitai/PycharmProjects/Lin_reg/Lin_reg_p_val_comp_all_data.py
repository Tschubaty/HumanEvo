import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model
from sklearn.datasets import load_digits
from scipy import stats

resolution = 10000
HIGHEST_VAR_NUM = 1000000
MOST_SIGNIFICANT_P_VAL = 10000


names = ['chr', 'I1116', 'I5319', 'I4532', 'I2862', 'I4529', 'I3758', 'I5725', 'I2861', 'I3957', 'I3565', 'I4792', 'I4775',
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



def lin_reg(origin):
    """
    Perform simple linear regression
    :param origin: address of csv file to be read
    :param chr_num: the number of chromosome being examined
    :return: array of p_values per position
    """
    df = pd.read_csv(origin, names=names)

    print("num of methylation positions",  df.shape[0])
    df = handle_NaNs(df)
    print("NAN's handles")

    df = best_Var(df)
    print("best variance taken")

    df = df.reset_index().set_index('index', drop=False)
    print(df)
    data = df.values


    X_train = np.asarray(np.matrix(data))
    Y_train = time_vec

    # # basic linear regression using sklearn

    p_value_arr = []

    y = Y_train.reshape(-1, 1)
    for i, meth_position in enumerate(X_train):
        x = meth_position[2:].reshape(-1, 1)

        # this is the actual linear regression
        if (np.count_nonzero(x[2:]) == 0):
            pass
        else:
            regr = LinearRegression().fit(x, y)
            intercept, r_value, p_value, slope, chr_index= regr.intercept_, regr.score(x, y) \
                , regr.p, regr.coef_, meth_position[1]

            # insert values into the array as tuples in the following order:
            # (chr num, index of meth position, p, r**2 )
            chr_num = meth_position[0]
            p_value_arr.append((int(chr_num), int(chr_index), p_value[0][0], r_value ** 2))

    #sort by most significant p value
    ret_val = sorted(p_value_arr, key=lambda x: x[2])
    # reduce data to most significant
    ret_val = ret_val[:MOST_SIGNIFICANT_P_VAL]
    # resort by chromosomes
    ret_val = sorted(ret_val, key=lambda x: x[0])

    return ret_val


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
    origin = "Z://Matlab_Data//Matlab_data_28.2.19//Unified_data//meth_unified.csv"
    output = "Z://Lin_reg//Graphs_lin_reg_4.3//all_data_most_sig_p//all_chr_lin_reg_output.csv"

    tuples = lin_reg(origin)

    print("finished linear regression")
    df = pd.DataFrame(tuples, columns=['position', 'chr', 'p', 'r2'])
    chr = df['chr']
    df.drop(labels=['chr'], axis=1, inplace=True)
    df.insert(0, 'chr', chr)

    df.to_csv(output, sep=',', index=False, header=None)


main()