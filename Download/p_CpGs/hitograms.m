% load data
load f_1116.mat         % ancient sample
load simu_1116.mat      % simulated sample
load RefBone5           % modern reference
clear simu_as
mms = mms.merge;
mms = mms.scale;

% recompute methylation
fas_mod = fas.reconstructmethylation('function','logistic');
simu_mod = simu_fas.reconstructmethylation('function','logistic');

% look at a specific chromosome
CHR_NAME = 'chr21';
anc = fas.getmethylation(CHR_NAME);
mod = mms.getmethylation(CHR_NAME)';
sim = simu_fas.getmethylation(CHR_NAME);
mod_smooth = mms.smooth(CHR_NAME,...
    simu_fas.methylation.win_size(simu_fas.indexofchr(CHR_NAME)));
anc_mod = fas_mod.getmethylation(CHR_NAME);
sim_mod = simu_mod.getmethylation(CHR_NAME);

% plot figures
bins = linspace(0,1,100);
figure, N = hist(anc,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('ancient')
figure, N = hist(mod,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('modern')
figure, N = hist(sim,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('simulated')
figure, N = hist(mod_smooth,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('modern smoothed')
figure, N = hist(anc_mod,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('ancient sigmoid')
figure, N = hist(sim_mod,bins); N = N / sum(N); bar(bins,N); ylim([0 0.45]); title('simulated sigmoid')

% compute distance matrix
order = {'modern','simulation','sample'};
dstmat = zeros(3);
DIST_TYPE = 'euclid';
if strcmp(DIST_TYPE,'euclid')
    dstmat(1,2) = sqrt(nansum( (mod_smooth - sim_mod).^2 ));
    dstmat(1,3) = sqrt(nansum( (mod_smooth - anc_mod).^2 ));
    dstmat(2,3) = sqrt(nansum( (sim_mod - anc_mod).^2 ));
    dstmat(2,1) = dstmat(1,2);
    dstmat(3,1) = dstmat(1,3);
    dstmat(3,2) = dstmat(2,3);
end
pt = seqneighjoin(dstmat,'equivar',order);
plot(pt)