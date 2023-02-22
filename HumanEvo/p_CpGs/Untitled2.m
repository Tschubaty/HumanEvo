% 
% 
% gc = load('gc_CpGs.mat');
% gc = gc.gc_CpGs;
% 
% chr_names = gc.get('chr_names');
% chr_names(end) = [];
% 
load 'RefBone5'
% 
% max_coverage = 100;
% max_tOct = 0.25;
% 
% 
% % 
% % samples = {'0708','1053','1116','1496','1507','1631','1632','1633','1734',...
% %     '1960','1961','1962','1965','2134','2139','2514','2520','2861',...
% %     '2862','2935','2978','2980','3133','3255','3565','3758','3957','4315',...
% %     '4432','4438','4529','4532','4596','4634','4775','4792','4873','4875',...
% %     '4877','4878','4914','5077','5233','5235','5236','5319','5725','5742',...
% %     '5743','5748','5835','5838','5950','7421'};
% % g_names = {'farmer','farmer','farmer','farmer','HG','farmer','farmer',...
% %     'farmer','farmer','HG','HG','HG','HG','HG','HG','farmer','farmer','farmer',...
% %     'farmer','farmer','farmer','farmer','farmer','farmer','farmer','farmer',...
% %     'farmer','farmer','HG','HG','farmer','farmer','HG','farmer','farmer','farmer',...
% %     'HG','HG','HG','HG','HG','farmer','HG','HG','HG','farmer','farmer','farmer',...
% %     'farmer','farmer','farmer','farmer','farmer','farmer'};
% 
% samples = {'1116', '1496', '1507','2520', '4873', '4875', '4877', '4878', '4914', '5077', '5233', ...
%     '5235', '5236', '5725'};
samples = {'1116', '1496'};
% g_names = {'farmer','farmer', 'HG','farmer', 'HG',...    
%     'HG','HG','HG','HG','farmer','HG','HG','HG','farmer'};
% 
no_samples = length(samples);
% samplist = cell(1,no_samples);
% for samp = 1:no_samples
%     SAMPLE = samples{samp};
%     fprintf(1,'loading sample I%s ...',SAMPLE);
%     load(sprintf('simu_%s',SAMPLE),'simu_fas');
%     samplist{samp} = simu_fas;
%     fprintf(1,' done\n');
% end
% 
% for delta = 0.1
%     dms = DMRs;
%     dms = dms.groupDMRs(samplist,g_names,gc,'chromosomes',chr_names,...
%         'min_finite',[0.8 0.8],'delta',delta);
%     save (['DMR_54HGXfarm',num2str(delta*100),'delta'], 'dms')
% end
% clear ('samplist')
% 
% samplist = cell(1,no_samples);
% for samp = 1:no_samples
%     SAMPLE = samples{samp};
%     fprintf(1,'loading sample I%s ...',SAMPLE);
%     load(sprintf('f_%s',SAMPLE),'fas');
%     samplist{samp} = fas;
%     fprintf(1,' done\n');
% end
% 
% for delta = 0.1
%     dm_unfilt = DMRs;
%     dm_unfilt = dm_unfilt.groupDMRs(samplist,g_names,gc,'chromosomes',chr_names,...
%         'min_finite',[0.8 0.8],'delta',delta);
%     save (['DMR_54HGXfarm',num2str(delta*100),'delta'], 'dm_unfilt','-append')
% end

bins = linspace(0,1,101);

for sampl = 1:no_samples
    SAMPLE = samples{sampl};
    fprintf(1,'loading sample I%s ...',SAMPLE);
    load(sprintf('simu_log_%s',SAMPLE),'simu_fas');
    load(sprintf('f_log_%s',SAMPLE),'fas');
    W = simu_fas.methylation.win_size;
    methmms=mms.getmethylation;
    
    for chr = 1:22
        methmms{chr} = nanmerge(methmms{chr},'ave');
        mms_local{chr} = nansmooth(methmms{chr},W(chr),'same');
    end
    
    tmodern = mms_local; modern = [];
    tsimu = simu_fas.getmethylation; simu = [];
    tsamp = fas.getmethylation; samp = [];
    coverage = [];
    
    for chr = 1:22
%         summer_simusamp (chr)= nansum ((simu{chr}-samp{chr}).^2);
%         summer_simumodern (chr)= nansum ((simu{chr}-modern{chr}).^2);
%         summer_sampmodern  (chr)= nansum ((samp{chr}-modern{chr}).^2);
        modern = [modern;tmodern{chr}./100];
        simu = [simu;tsimu{chr}];
        samp = [samp;tsamp{chr}];
        coverage = [coverage fas.getNo_Cs(chr) + fas.getNo_Ts(chr)];
    end
    
%     simusamp (sampl) = sqrt(sum(summer_simusamp));
%     simumodern (sampl) = sqrt(sum(summer_simumodern));
%     sampmodern (sampl) = sqrt(sum(summer_sampmodern));

%     simusamp(sampl) = nansum( (simu - samp).^2 ./ (simu + samp) );
%     simumodern(sampl) = nansum( (simu - modern).^2 ./ (simu + modern) );
%     sampmodern(sampl) = nansum( (modern - samp).^2 ./ (modern + samp) );

%     Nsimu = hist(simu,bins); Nsimu = Nsimu / sum(Nsimu);
%     Nsamp = hist(samp,bins); Nsamp = Nsamp / sum(Nsamp);
%     Nmodern = hist(modern,bins); Nmodern = Nmodern / sum(Nmodern);

    simusamp(sampl) = sqrt( nansum( (simu - samp).^2 ) );
    simumodern(sampl) = sqrt( nansum( (simu - modern).^2 ) );
    sampmodern(sampl) = sqrt( nansum( (modern - samp).^2) );
    
%     [~, simusamp(sampl)] = testchi2hists(Nsimu,Nsamp,'noSamples','nonequal');
%     [~, simumodern(sampl)] = testchi2hists(Nsimu,Nmodern,'noSamples','nonequal');
%     [~, sampmodern(sampl)] = testchi2hists(Nmodern,Nsamp,'noSamples','nonequal');

    fprintf(1,' done\n');
end

idx_samp = 2;
dstmat = [0 simusamp(idx_samp) simumodern(idx_samp); ...
    simusamp(idx_samp) 0 sampmodern(idx_samp); ...
    simumodern(idx_samp) sampmodern(idx_samp) 0];
ph = seqneighjoin(dstmat,'equivar',{'simu','samp','mod'});
view(ph)


