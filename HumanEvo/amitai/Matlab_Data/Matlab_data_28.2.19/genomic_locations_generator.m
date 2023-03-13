format long

load('C:\Users\amitai\Google Drive\p_CpGs\gc_CpGs.mat');

chr1 = gc_CpGs.coordinates{1,1};
chr2 = gc_CpGs.coordinates{1,2};
chr3 = gc_CpGs.coordinates{1,3};
chr4 = gc_CpGs.coordinates{1,4};
chr5 = gc_CpGs.coordinates{1,5};
chr6 = gc_CpGs.coordinates{1,6};
chr7 = gc_CpGs.coordinates{1,7};
chr8 = gc_CpGs.coordinates{1,8};
chr9 = gc_CpGs.coordinates{1,9};
chr10 = gc_CpGs.coordinates{1,10};
chr11 = gc_CpGs.coordinates{1,11};
chr12 = gc_CpGs.coordinates{1,12};
chr13 = gc_CpGs.coordinates{1,13};
chr14 = gc_CpGs.coordinates{1,14};
chr15 = gc_CpGs.coordinates{1,15};
chr16 = gc_CpGs.coordinates{1,16};
chr17 = gc_CpGs.coordinates{1,17};
chr18 = gc_CpGs.coordinates{1,18};
chr19 = gc_CpGs.coordinates{1,19};
chr20 = gc_CpGs.coordinates{1,20};
chr21 = gc_CpGs.coordinates{1,21};
chr22 = gc_CpGs.coordinates{1,22};


fn = strcat('genomic_locations.csv') ;

T = padcat(chr1, chr2 ,chr3 ,chr4 ,chr5 ,chr6 ,chr7 ,chr8 ,chr9 ,chr10 ,chr11 ,chr12 ,chr13 ,chr14 ,chr15 ,chr16 ,chr17 ,chr18 ,chr19 ,chr20 ,chr21 ,chr22);


dlmwrite(fn, T, 'precision', '%i');




