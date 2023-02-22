gc = load('gc_CpGs.mat');
gc = gc.gc_CpGs;

chr_names = gc.get('chr_names');
chr_names(end) = [];

samples = {'1116', '1496', '1507','2520', '4873', '4875', '4877', '4878', '4914', '5077', '5233', ...
    '5235', '5236', '5725'};
g_names = {'farmer','farmer', 'HG','farmer', 'HG',...    
    'HG','HG','HG','HG','farmer','HG','HG','HG','farmer'};
sex = {'M','M','M','M','F','F','F','M','M','M','F','M','M','M'};
Location = {'Serbia','Hungary','Hungary', 'Bulgaria',  'Serbia', 'Serbia', 'Serbia','Serbia','Serbia', 'Croatia','Serbia',...
    'Serbia','Serbia',  'Croatia'};

no_samples = length(samples);
samplist = cell(1,no_samples);
for samp = 1:no_samples
    SAMPLE = samples{samp};
    fprintf(1,'loading sample I%s ...',SAMPLE);
    load(sprintf('f_%s',SAMPLE),'fas');
    samplist{samp} = fas;
    fprintf(1,' done\n');
end


dm = DMRs;
dm = dm.groupDMRs(samplist,g_names,gc,'chromosomes',chr_names,...
    'min_finite',[0.8 0.8],'delta',0.3,'min_CpGs',12);
dm = dm.annotate('genes',gi_genes,'cgis',gi_CGIs);