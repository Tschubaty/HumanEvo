import sys
import csv

# this is for getting the arrays and dictionaries
sys.path.insert(0, r'X://PycharmProjects//Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]




from pandas import DataFrame
import pandas as pd
from scipy.stats import uniform
from scipy.stats import randint
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import csv
from matplotlib import colors

#origin = "X://Matlab_Data//Matlab_data_28.2.19//all_data_most_sig_p//all_chr_lin_reg_output.csv"
origin = "X://Matlab_Data//Matlab_data_28.2.19//working_example//all_chr_lin_reg_output_we.csv"
CHR_NUM = 22
MAX_MLOG = 25
NUM_MOST_SIGNIFICANT = 100
cmap = colors
x = cm.get_cmap('tab20b', 22)
MAX_P_VAL = 1

print('kaka')
max = m.num_cpgs_per_chr

# -log_10(pvalue)
df = pd.read_csv(origin, names=['chromosome', 'position', 'pvalue'])
ndf = pd.DataFrame({'chromosome':[], 'position':[], 'pvalue':[]})

#
# for i in range(len(max)):
#     for j in range(max[i]):
#         ndf.loc[sum(max[:i]) + j ] = [int(i+1), int(j+1), MAX_P_VAL]
#
# print('here')
# ndf[['chromosome']] = ndf[['chromosome']].astype(int)
# ndf[['position']] = ndf[['position']].astype(int)
#
#
# for i, row in enumerate(df):
#     ndf.loc[((ndf.chromosome == df.chromosome.iloc[i]) & (ndf.position == df.position.iloc[i])), 'pvalue'] = df.iloc[i].pvalue
#
# df = ndf
#
#
#





print(df.iloc[:3])
# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['-log 10 p value'] = -np.log10(df['pvalue'])
# df = df.sort_values('chromosome')

df['ind'] = range(1, len(df)+1)
df_grouped = df.groupby(('chromosome'))



fig = plt.figure()
ax = fig.add_subplot(111)
colors = [x(i/CHR_NUM) for i in range(CHR_NUM)]
print(colors)
x_labels = []
x_labels_pos = []


most_significant = df.nlargest(NUM_MOST_SIGNIFICANT, '-log 10 p value')
print(most_significant)
most_significant.drop("ind", inplace=True, axis=1)



for num, (name, group) in enumerate(df_grouped):

    group.plot(kind='scatter', x='ind', y='-log 10 p value',color=colors[num%len(colors)], ax=ax, s=1)
    x_labels.append(name)
    x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, MAX_MLOG])
ax.set_xlabel('Chromosome Number')
plt.savefig('X://Graphs//Manhattan_graph_7.3//fig.png', dpi=4000)
plt.show()



most_significant.to_csv('mahnattan_output.csv', sep=',', index=False, header=False, encoding='utf-8' )