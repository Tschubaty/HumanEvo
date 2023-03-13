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
INTERVAL = 1 # interval in thousands of years
START_METH = 1
END_METH = 0 # this is for not including the last sample




# this is for getting the arrays and dictionaries
sys.path.insert(0, r'X://PycharmProjects//Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]


#
# # for the LUSTRE
# sys.path.insert(0, r'Organization_code')
# m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
# del sys.path[0]

# interval names
#['1-2K_5_samples', '2-3K_3_samples', '3-4K_7_samples', '4-5K_10_samples', '5-6K_7_samples', '6-7K_3_samples', '7-8K_11_samples', '8-9K_10_samples', '9-10K_2_samples', '10-11K_2_samples', '12-13K_1_samples', '45-46K_1_samples']
#interval mean years
#[1.32, 2.5, 3.414285714285714, 4.48, 5.299999999999999, 6.099999999999999, 7.454545454545456, 8.21, 9.1, 10.4, 12.0, 45.0]
# interval num samples
#[5, 3, 7, 10, 7, 3, 11, 10, 2, 2, 1, 1]



interv_names = []
interv_time_vec = []
interv_num_samples = []



def create_intervals(times):
    cur_interv = 1
    num_samples_cur_interval = 0
    for i, time in enumerate(times):
        if((time< cur_interv + INTERVAL) and (time>= cur_interv)):

            num_samples_cur_interval +=1
        else:
            if(num_samples_cur_interval != 0):
                interv_names.append(str(cur_interv) + '-' + str(cur_interv + INTERVAL) + 'K_'+str(num_samples_cur_interval) + '_samples')
                interv_time_vec.append(np.mean(times[i - num_samples_cur_interval: i]))
                interv_num_samples.append(num_samples_cur_interval)
                cur_interv = int(np.min(times[i:]))
                num_samples_cur_interval = 1

        if (i == len(times) - 1):
            interv_names.append(str(cur_interv) + '-' + str(cur_interv + INTERVAL) + 'K_' + str(num_samples_cur_interval) + '_samples')
            interv_time_vec.append(np.mean(time))
            interv_num_samples.append(1)







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
        y = interv_time_vec
        for i, row in enumerate(df):
            nans = False
            chr = int(row[0])
            for i, elem in enumerate(row):
                if elem == '':
                    nans = True
            if(not nans):
                full = np.array(row[START_METH:]).astype(np.float64)
                #divide the row into intervals
                cur_loc = 0
                x = []
                for elem in interv_num_samples:
                    x.append(np.mean(full[cur_loc:cur_loc + elem]))
                    cur_loc+=elem
                x = np.array(x)
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
    df.to_csv(output, sep=',', index=False, header=None)



def main():
    # working examples
    origin = "meth_unified.csv"
    output ="1K_intervals_all.csv"

    lin_reg(origin, output)


create_intervals(m.time_vec)
main()
