
import sys

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

origin = "X://Matlab_Data//Matlab_data_28.2.19//all_data_most_sig_p//all_chr_lin_reg_output.csv"
CHR_NUM = 22
MAX_MLOG = 25
SIG_VAL = 0.05
cmap = colors
x = cm.get_cmap('tab20b', 22)

max = m.num_cpgs_per_chr



# -log_10(pvalue)
df = pd.read_csv(origin, names=['chromosome', 'position', 'pvalue'])
N = df.shape[0]
# How to plot gene vs. -log10(pvalue) and colour it by chromosome?
df['ind'] = range(1, len(df) + 1)
df['-log 10 p value'] = -np.log10(df.pvalue)
df.sort_values('pvalue', inplace=True)
df['sorted index'] = range(len(df))
corr_q_Vals = [((N*df['pvalue'][i])/(i+1)) for i in range(0,len(df))]
adj_q_val = []
cur_min = 1
for q in corr_q_Vals[::-1]:
    if(q<cur_min):
        cur_min = q
    adj_q_val.append(cur_min)

adj_q_val = adj_q_val[::-1]
df['adj q value'] = adj_q_val
df['-log 10 adj q value'] = -np.log10(adj_q_val)

for index, elem in df.iloc[::-1].iterrows():
    # print("pval", elem['pvalue'])
    # print("bound", (elem['sorted index'] * SIG_VAL) / len(df))
    if((elem['pvalue']) < (elem['sorted index'] * SIG_VAL) / len(df)):
        print("success")
        df = df.iloc[:index]
        break

print(df)

print("num significant vals", len(df))
df_grouped = df.groupby(('chromosome'))







fig = plt.figure()
ax = fig.add_subplot(111)
colors = [x(i/CHR_NUM) for i in range(22)]

x_labels = []
x_labels_pos = []
prev_loc = 0
for num, (name, group) in enumerate(df_grouped):
    group.plot(kind='scatter', x='ind', y='-log 10 p value',color=colors[num%len(colors)], ax=ax)
    x_labels.append(name)
    x_labels_pos.append((len(group) - (len(group)- prev_loc)/2))
    prev_loc = x_labels[-1]
ax.set_xticks(x_labels_pos)
ax.set_xticklabels(x_labels)
ax.set_xlim([0, len(df)])
ax.set_ylim([0, MAX_MLOG])
ax.set_xlabel('Chromosome Number')

plt.show()