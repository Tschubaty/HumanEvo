import numpy as np
import pandas as pd

UCSC_col_names = ['chr', 'f_pos', 'l_pos', 'name', 'full_name', 'strand']
interval_proc_names = ['chr', 'CpG_pos', 'p_all', '-log10_p_all', 'p_n_ust', '-log10_p_n_ust', '-log10_p_min']
genomic_locations_names = ['chr_1', 'chr_2', 'chr_3', 'chr_4', 'chr_5', 'chr_6', 'chr_7', 'chr_8', 'chr_9', 'chr_10', 'chr_11', 'chr_12', 'chr_13', 'chr_14', 'chr_15', 'chr_16', 'chr_17', 'chr_18', 'chr_19', 'chr_20', 'chr_21', 'chr_22']
genomic_locations_src = 'genomic_locations.csv'
interval_proc_src = 'interval_procedure_output_thresh.csv'




def get_positions(gene):
    if(gene['strand'] == '+'):
        first = gene['f_pos'] - 5000
        last = gene['l_pos']
    else:
        first = gene['f_pos']
        last = gene['l_pos'] + 5000
    return first,last



def getGenes(src, sig_positions):
    output_arr = []
    not_in_gene_arr = []
    # gene locs is a tuple with ('chr#', genomic position on chromosome)
    # src is the source of the UCSC gene name xlsx file
    df = pd.read_excel(src, names = UCSC_col_names)
    df_grouped = df.groupby('chr')
    found = False
    for pos in sig_positions:
        chr_genes = df_grouped.get_group(pos[0])
        for i, gene in chr_genes.iterrows():
            first, last = get_positions(gene)
            if((pos[1] > first) and (pos[1] < last)):
                output_arr.append([pos[0], pos[1], gene['name'], gene['full_name'], gene['strand']])
                found = True
        if(not found):
            not_in_gene_arr.append([pos[0], pos[1]])
        found = False
    return output_arr, not_in_gene_arr


def get_genomic_pos(CpG_pos_file, gen_pos_file):
    # CpG_pos_file is the file of the output of the interval procedure
    # gen_pos_file is a file matching a CpG position to a general genomic position
    df1 = pd.read_csv(CpG_pos_file, names = interval_proc_names)
    df2 =  pd.read_csv(gen_pos_file, names = genomic_locations_names)
    gen_pos_tuples_arr = []
    for i, row in df1.iterrows():
        chr = int(row['chr'])
        pos = int(row['CpG_pos'])
        gen_pos = df2.iloc[pos - 1, chr - 1]
        gen_pos_tuples_arr.append((chr, gen_pos))

    return gen_pos_tuples_arr


pos_arr = get_genomic_pos(interval_proc_src,genomic_locations_src )
in_gene, not_in_gene = getGenes('UCSC_Genes.xlsx', pos_arr)


pd.DataFrame(in_gene).to_csv("sig_vals_intervals.csv", header=None, index=False)
