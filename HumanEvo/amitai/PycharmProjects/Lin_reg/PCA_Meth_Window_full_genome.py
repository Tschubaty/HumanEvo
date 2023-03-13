import sys

# # this is for getting the arrays and dictionaries
# sys.path.insert(0, r'X://PycharmProjects//Organization_code')
# m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
# del sys.path[0]

# for the LUSTRE
sys.path.insert(0, r'Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]



import pandas as pd
from matplotlib import cm

from sklearn import preprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np


GRAPH_ORIGIN = [0, 0]

names = ['chr'] + m.names


def legend_generator():
    """
    generate the legend of the final graph (concatenate the name with the years
    :return: array of concatenated names + time_vec elements
    """
    arr = []
    for i in range(len(m.names)):
        arr.append(m.names[i] + " - " + str(m.time_vec[i]) + "K years")
    return arr


def construct_coverage_vec():
    """
    This is a method for constructing a coverage vector from the coverage dictionary
    :return: a coverage vector in the correct order
    """
    arr = []
    for elem in m.names:
        arr.append(m.cov_dict[elem])
    return arr


def projection(a, b):
    return ((np.dot(a, b)) / np.linalg.norm(a))


def time_and_coverage(x, pc_1, pc_2):
    """
    Find the vector of the principal component portion of the projection of the methylation on the time and coverage
    :param x: the methylation vector
    :return:
    """
    time = np.array(m.time_vec)
    coverage = construct_coverage_vec()
    time_proj = []
    coverage_proj = []
    for vec in x:
        vec_proj_time = projection(time, vec)
        time_proj.append(vec_proj_time)
        vec_proj_coverage = projection(coverage, vec)
        coverage_proj.append(vec_proj_coverage)
    return projection(pc_1, time_proj), projection(pc_2, time_proj), projection(pc_1, coverage_proj), projection(pc_2, coverage_proj)


def main():
    """
    This is the main function, reads the csv and does principle component analysis on it
    :return:
    """

    # origin_arr = ["X://Matlab_Data//training_data//chr_1_training.csv"] # for training
    # origin_arr = ["X://Matlab_Data//Matlab_data_28.2.19//Unified_data//meth_unified.csv"]
    origin_arr = ["meth_unified.csv"]  # for LUSTRE



    output_arr = []
    for i in range(len(origin_arr)):
        # output_arr.append(
        #     "X://PCA//outputs//chr_wise_two_comp_PCA_correlation//meth_unified.png")
        output_arr.append("meth_unified.csv")

    for i in range(len(origin_arr)):
        PCA_x(origin_arr[i], output_arr[i])


def PCA_x(origin, output):
    """
    Perform pca on the csv file in the origin and output a graph into the file
    :param origin: the origin of the csv file
    :param output: address to output the file
    :param chr: chromosome being examined
    :return: nothing
    """
    df = pd.read_csv(origin, names=names)


    df = df.drop(['chr'], axis=1) #drop the chromosome number
    df = handle_NaNs(df)

    x = df.values

    x_scaled = preprocessing.scale(x, with_mean=True, with_std=True)  # added 21.3

    pca = PCA(n_components=2)

    principalComponents = pca.fit_transform(np.transpose(x_scaled))  # added 21.3
    component_vectors = pca.components_

    exp_variance = pca.explained_variance_ratio_

    principalDf = pd.DataFrame(data=principalComponents
                               , columns=['principal component 1', 'principal component 2'])

    headers = pd.DataFrame({'humans': m.names})

    finalDf = pd.concat([principalDf, headers], axis=1)

    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(1, 1, 1)

    # plt.ylim(-1, 1)
    # plt.xlim(-1, 1)
    ax.set_xlabel('Dim1 (' + str(exp_variance[0] * 100)[:8] + '%' + ')', fontsize=15)
    ax.set_ylabel('Dim2 (' + str(exp_variance[1] * 100)[:8] + '%' + ')', fontsize=15)
    ax.set_title('Two Component PCA - All Data', fontsize=20)
    targets = legend_generator()
    cmap = cm.get_cmap('tab20b', 8)
    for i, target in enumerate(targets):
        plt.scatter(principalDf['principal component 1'].iloc[i]
                    , principalDf['principal component 2'].iloc[i]
                    , s=30, c=[cmap(i / len(m.names))])

        text = ax.annotate(m.names[i], (principalDf.iloc[i, 0], principalDf.iloc[i, 1]))
        text.set_fontsize(9)

    # get the boundaries of the graph
    lat1, lat2, long1, long2 = principalDf['principal component 1'].min(), principalDf['principal component 1'].max(), principalDf[
        'principal component 2'].min(), principalDf['principal component 2'].max()
    xrange = np.arange(lat1, lat2, ((lat2 - lat1) / 10))
    ax.legend(targets, loc=2, prop={'size': 6}, )

    t1, t2, c1, c2= time_and_coverage(x, component_vectors[0], component_vectors[1])
    x1, y1 = 0,0
    x2, y2 = t1,t2
    line_eqn = lambda x: ((y2 - y1) / (x2 - x1)) * (x - x1) + y1
    # this is for drawing the time and coverage lines

    lines = []
    from matplotlib.legend import Legend
    lines.append(ax.plot(xrange, [line_eqn(x) for x in xrange], color='m', linestyle=':', linewidth = 1))
    text1 = plt.annotate('Time projection', xy=(xrange[-1], line_eqn(xrange[-1])), color='m')
    text1.set_fontsize(9)

    x2, y2 = c1, c2
    line_eqn = lambda x: ((y2 - y1) / (x2 - x1)) * (x - x1) + y1
    lines.append(ax.plot(xrange, [line_eqn(x) for x in xrange], color='c', linestyle=':', linewidth = 1))
    text2 = plt.annotate('Coverage projection', xy=(xrange[-1], line_eqn(xrange[-1])), color='c')
    text2.set_fontsize(9)
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
