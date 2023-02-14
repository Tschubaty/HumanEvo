% hard-coded parameters
% max_coverage = 100;
% max_tOct = 0.25;
% max_aOga = max_tOct;

% load coordinates
gc = load('gc_CpGs.mat');
gc = gc.gc_CpGs;

% load reference
% load RefBone5

% initialize
% as = amSample;

% chromosome names
chr_names = gc.get('chr_names');
chr_names(end) = [];

% sample list
samples = {'0708', '1053', '1116', '1496', '1507', '1631','1632', '1633', ...
    '1734', '1960', '1961', '1962', '1965', '2134', '2139', '2514',...
    '2520', '2861', '2862', '2935', '2978', '2980', '3133', '3255', '3565',...
    '3758', '3957', '4315', '4432', '4438', '4529', '4532', '4596', '4634',...
    '4775', '4792', '4873', '4875', '4877', '4878', '4914', '5077', '5233', ...
    '5235', '5236', '5319', '5725', '5742', '5743', '5748', '5835', '5838',...
    '5950', '7421' ,'Alt', 'Vin', 'Den', 'Ust', 'Los', 'LBK', 'LaB', 'M45',...
    'I15','Bal','SF1','Ust'};
g_names = {'farmer','farmer','farmer','farmer', 'HG','farmer','farmer', 'farmer',...
    'farmer','HG','HG','HG','HG','HG','HG','farmer','farmer','farmer',...
    'farmer','farmer','farmer','farmer','farmer','farmer','farmer','farmer','farmer',...
    'farmer','HG','HG','farmer','farmer','HG','farmer','farmer','farmer', 'HG',...    
    'HG','HG','HG','HG','farmer','HG','HG','HG','farmer','farmer','farmer',...
    'farmer','farmer','farmer','farmer','farmer','farmer','HG','HG','HG','HG','HG',...
    'farmer', 'HG', 'HG','farmer', 'HG', 'HG','HG'};
sex = {'M','M','M','M','M','F','M','F','M','F','M','M','M','M','F','M','M',...
    'F','M','M','M','F','M','F','F','F','M','M','M','M','M','M','M','F','M','M',...
    'F','F','F','M','M','M','F','M','M','M','M','F','F','M','M','M','M','M',...
    'F','F','F','M','M','M','M','F','M','M','F','M'};
Location = {'Turkey','Russia','Serbia','Hungary','Hungary', 'Armenia', 'Armenia',...
    'Armenia', 'Ukraine', 'Siberia', 'Siberia', 'Siberia','Siberia', ...
    'Siberia', 'Siberia', 'Iran', 'Bulgaria', 'Scotland', 'Scotland', 'Scotland', ...
    'Scotland', 'Scotland', 'Scotland', 'England', 'California', 'Spain', 'Russia', 'Uzbekistan', ...
    'Latvia', 'Latvia', 'Turkey','Turkey','Latvia', 'Turkmenistan', 'Kazakhstan', ...
    'Kazakhstan', 'Serbia', 'Serbia', 'Serbia','Serbia','Serbia', 'Croatia','Serbia',...
    'Serbia','Serbia', 'Alaska', 'Croatia', 'Turkey', 'Turkey' ,'Netherlands','Germany', ...
    'Spain', 'Ethiopia','Uzbekistan','Siberia','Croatia','Siberia','Siberia','Luxembourg',...
    'Germany','Spain','China','Turkey','South Africa','Sweden','Siberia'};

no_samples = length(samples);
samplist = cell(1,no_samples);
for samp = 1:no_samples
    SAMPLE = samples{samp};
    fprintf(1,'loading sample I%s ...',SAMPLE);
    load(sprintf('f_%s',SAMPLE),'fas');
    samplist{samp} = fas;
    fprintf(1,' done\n');
end
%}

%
% add bone reference
%{
mms = mms.mergepositions;
samplist{end+1} = mms;
no_samples = no_samples + 1;
g_names = [g_names {'Farmer'}];
clear fas mms
%}

%{
% plot a low-variance region (chr=9)
width = 10000;
idx = 1025790;
gidx = gc.coordinates{chr}(idx);
reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
plotregion(reg,gc,samplist);
%}

%{
% plot a sex-DMR (chr=9)
width = 700;
idx = 763708;
gidx = gc.coordinates{chr}(idx);
reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
plotregion(reg,gc,samplist);
for ii = 1:(no_samples-1)
    text(idx+22,ii+0.5,sex{ii});
end
text(idx+22,12.5,'F');
vals_male = nanmean(meth_male(:,idx-10:idx+10),2);
vals_female = nanmean(meth_female(:,idx-10:idx+10),2);
xx = [ones(1,8) 2*ones(1,4)];
figure, plot(xx,[vals_male; vals_female],'o','markerfacecolor',[0.8 0.8 0.8],...
    'markeredgecolor','k','markersize',8);
set(gca,'xlim',[0 3],'xticklabel',[]);
ylabel('Methylation');
text(1,0.8,'Female','horizontal','center');
text(2,0.6,'Male','horizontal','center');
title(sprintf('chr9:%s-%s',numcommas(108342986),numcommas(108343991)));
%}

%{
% plot a lifestyle-DMR (chr=9)
width = 3000;
idx = 416659;
gidx = gc.coordinates{chr}(idx);
reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
plotregion(reg,gc,samplist);
for ii = 1:(no_samples-1)
    text(idx+0,ii+0.5,g_names{ii});
end
vals_farmer = nanmean(meth_farmer(:,idx-13:idx+13),2);
vals_hg = nanmean(meth_hg(:,idx-13:idx+13),2);
xx = [ones(1,7) 2*ones(1,4)];
figure, plot(xx,[vals_hg; vals_farmer],'o','markerfacecolor',[0.8 0.8 0.8],...
    'markeredgecolor','k','markersize',8);
set(gca,'xlim',[0 3],'xticklabel',[]);
ylabel('Methylation');
text(1,0.4,'HGs','horizontal','center');
text(2,0.6,'Farmers','horizontal','center');
title(sprintf('chr9:%s-%s',numcommas(69893518),numcommas(69899553)));
%}

%{
% plot a ethnic-DMR (chr=9)
width = 3000;
idx = 212869;
gidx = gc.coordinates{chr}(idx);
reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
plotregion(reg,gc,samplist);
for ii = 1:(no_samples-1)
    text(idx+0,ii+0.5,Location{ii});
end
vals_farmer = nanmean(meth_farmer(:,idx-13:idx+13),2);
vals_hg = nanmean(meth_hg(:,idx-13:idx+13),2);
xx = [ones(1,7) 2*ones(1,4)];
figure, plot(xx,[vals_hg; vals_farmer],'o','markerfacecolor',[0.8 0.8 0.8],...
    'markeredgecolor','k','markersize',8);
set(gca,'xlim',[0 3],'xticklabel',[]);
ylabel('Methylation');
text(1,0.4,'HGs','horizontal','center');
text(2,0.6,'Farmers','horizontal','center');
title(sprintf('chr9:%s-%s',numcommas(69893518),numcommas(69899553)));
%}

% text(idx+22,12.5,'F');
% vals_male = nanmean(meth_male(:,idx-10:idx+10),2);
% vals_female = nanmean(meth_female(:,idx-10:idx+10),2);
% xx = [ones(1,8) 2*ones(1,4)];
% figure, plot(xx,[vals_male; vals_female],'o','markerfacecolor',[0.8 0.8 0.8],...
%     'markeredgecolor','k','markersize',8);
% set(gca,'xlim',[0 3],'xticklabel',[]);
% ylabel('Methylation');
% text(1,0.8,'Female','horizontal','center');
% text(2,0.6,'Male','horizontal','center');
% title(sprintf('chr9:%s-%s',numcommas(108342986),numcommas(108343991)));

%{
% create methylation data matrix
meth = zeros(no_samples,length(samplist{1}.getNo_Ts(chr)));
for samp = 1:no_samples
    meth(samp,:) = samplist{samp}.getmethylation(chr);
end
% PCA
NO_VARS = 100;
dat = meth(1:end-1,:)';
idx = find(isnan(sum(dat,2)));
dat(idx,:) = [];
no_variables = size(dat,1);
to_keep = randperm(no_variables,NO_VARS);
dat = dat(to_keep,:);
vsm = vsmatrix(dat);
vsm = vsm.set('col_sampleset',samples);
% compute variance
vr = nanvar(meth(1:no_samples-1,:),0,1);
% compute average variance in a window
%tmplt = 1/W*ones(1,W);
%vr = nanconv(vr,tmplt,'same');

% choose a high-variance region
%width = 10000;
%idx = 956830;
%gidx = gc.coordinates{chr}(idx);
%reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
reg = 'chr3:134723533-134724663';
plotregion(reg,gc,samplist);
% choose sex difference
idx_male = strmatch('M',sex);
idx_female = [strmatch('F',sex); no_samples];
%meth_male = mean(meth(idx_male,:),1);
%meth_female = mean(meth(idx_female,:),1);
meth_male = meth(idx_male,:);
meth_female = meth(idx_female,:);
no_positions = size(meth,2);
h = nan(1,no_positions);
for ii = 1:10:no_positions
    h(ii) = ttest2(meth_male(:,ii),meth_female(:,ii));
end
meth_diff = meth_male - meth_female;
W = 500; tmplt = ones(1,W);
meth_diff = nanconv(meth_diff,tmplt,'same');
width = 3000;
idx = 763708;
%idx = 411070;
gidx = gc.coordinates{chr}(idx);
reg = sprintf('chr%d:%d-%d',chr,gidx-width,gidx+width);
plotregion(reg,gc,samplist);
for ii = 1:(no_samples-1)
    text(idx+40,ii+0.5,sex{ii});
end
% HG-farmers
idx_farmer=strmatch('Farmer',g_names);
idx_hg=strmatch('HG',g_names);
idx_croatia=strmatch('Serbia',Location);
%}

% detectDMRs
dm = DMRs;
dm = dm.groupDMRs(samplist,g_names,gc,'chromosomes',chr_names,...
    'min_finite',[2 3],'delta',0.35);