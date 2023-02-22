gi_genes = gintervals;
gi_genes = gi_genes.readbed('UCSC_Genes.txt','strandcolumn',6,...
    'inamecolumn',4,'mdname','UCSCName','mdcol',5);
gi_genes = gi_genes.removeduplicates;
%save gi_genes gi_genes