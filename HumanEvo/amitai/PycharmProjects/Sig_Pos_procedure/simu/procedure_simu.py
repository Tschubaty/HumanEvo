import intervals_without_ust_simu
import intervals_all_simu



import sys
import pandas as pd
import numpy as np



# this is for getting the arrays and dictionaries
sys.path.insert(0, r'X://PycharmProjects//Organization_code')
m = __import__('Arrays_and_dictionaries', globals(), locals(), [])
del sys.path[0]




def procedure(origin_1, origin_2,  output): #todo if is ust is included in simu add origin_2
    df1 = pd.read_csv(origin_1, names=['chromosome', 'position', 'pvalue'], dtype=float)
    df2 = pd.read_csv(origin_2, names=['chromosome', 'position', 'pvalue'], dtype=float)
    df1['-log 10 p value'] = -np.log10(df1['pvalue'])
    df2['-log 10 p value'] = -np.log10(df2['pvalue'])
    analysis = []
    for (i, row1), (j,row2) in zip(df1.iterrows(), df2.iterrows()):
        # if(row1['-log 10 p value']>THRESHOLD):
        min_log_10_p_val = min( row1['-log 10 p value'], row2['-log 10 p value'])
        analysis.append([row1['chromosome'],  row1['position'], row1['pvalue'], row1['-log 10 p value'], row2['pvalue'], row2['-log 10 p value'], min_log_10_p_val  ])

    analysis.sort(key=lambda i:i[-1]) # sort by minimal threshold
    write_csv(analysis, output)



def write_csv(arr, output):
    df = pd.DataFrame(arr)
    df.to_csv(output, sep=',', index=False, header=None)




def main():
    # working examples
    origin1 = "1K_intervals_all_simu.csv"
    origin2 ="1K_intervals_without_ust_simu.csv"
    output ="interval_procedure_output_simu.csv"
    procedure(origin1, origin2, output)



main()