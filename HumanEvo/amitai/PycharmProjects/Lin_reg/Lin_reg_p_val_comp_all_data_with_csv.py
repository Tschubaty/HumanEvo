import csv
import numpy as np
import pandas as pd
from sklearn import linear_model
from scipy import stats
import sys

resolution = 10000
HIGHEST_VAR_NUM = 10000
MOST_SIGNIFICANT_P_VAL = 10000
BIG_NUM = 99999


# for the LUSTRE
sys.path.insert(0, r'Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]



names = m.names

time_vec = m.time_vec

# class LinearRegression(linear_model.LinearRegression):
#     """
#     LinearRegression class after sklearn's, but calculate t-statistics
#     and p-values for model coefficients (betas).
#     Additional attributes available after .fit()
#     are `t` and `p` which are of the shape (y.shape[1], X.shape[1])
#     which is (n_features, n_coefs)
#     This class sets the intercept to 0 by default, since usually we include it
#     in X.
#     """
#
#     def init(self, *args, **kwargs):
#         if not "fit_intercept" in kwargs:
#             kwargs['fit_intercept'] = False
#         super(LinearRegression, self) \
#             .init(*args, **kwargs)
#
#     def fit(self, X, y, n_jobs=1):
#         self = super(LinearRegression, self).fit(X, y, n_jobs)
#
#         sse = np.sum((self.predict(X) - y) ** 2, axis=0) / float(X.shape[0] - X.shape[1])
#         se = np.array([
#             np.sqrt(np.diagonal(sse[i] * np.linalg.inv(np.dot(X.T, X))))
#             for i in range(sse.shape[0])
#         ])
#
#         self.t = self.coef_ / se
#         self.p = 2 * (1 - stats.t.cdf(np.abs(self.t), y.shape[0] - X.shape[1]))
#         return self

def lin_reg(origin, output):
    """
    Perform simple linear regression
    :param origin: address of csv file to be read
    :param chr_num: the number of chromosome being examined
    :return: array of p_values per position
    """
    # origin =

    with open(origin) as csvfile:
        df = csv.reader(csvfile, delimiter=',')
        analysis = []
        added = 0
        smallest_var = BIG_NUM
        position_counter = np.full((22), 1)
        # y = time_vec.reshape(-1,1)
        y = time_vec
        for i, row in enumerate(df):
            nans = False
            chr = int(row[0])
            for i, elem in enumerate(row):
                if elem == '':
                    nans = True
            if(not nans):
                x = np.array(row[1:])
                x = x.astype(np.float64)

                # # this is for the method of calculation using sklearn
                # x = x.reshape(-1,1)
                # regr = LinearRegression().fit(x, y)
                # p_value, chr_index , pos= regr.p[0][0], chr, position_counter[chr-1]
                p_value , chr_index , pos= stats.linregress(x,y)[3], chr, position_counter[chr-1]
                analysis.append([chr_index , pos, p_value])
            position_counter[int(chr) -1] +=1

    write_csv(analysis, output)


def write_csv(arr, output):
    df = pd.DataFrame(arr)
    print("df")
    print(df)

    df.to_csv(output, sep=',', index=False, header=None)



def main():
    # working examples
    origin = "X://Matlab_Data//Matlab_data_28.2.19//working_example//meth_unified.csv"
    output ="X://Matlab_Data//Matlab_data_28.2.19//working_example//all_chr_lin_reg_output.csv"

    # origin = "X://Matlab_Data//Matlab_data_28.2.19//Unified_data//meth_unified.csv"
    # output = "X://Matlab_Data//Matlab_data_28.2.19//all_data_most_sig_p//all_chr_lin_reg_output.csv"

    lin_reg(origin, output)



main()

