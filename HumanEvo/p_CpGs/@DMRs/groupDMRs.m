function [obj, Qt_up, Qt_down] = groupDMRs(~,samples,groups,coord,varargin)

% groupDMRS detects DMRs between two groups of samples.
% -------------------------------------------------------------------
% [dm, Qt_up, Qt_down] = dm.groupDMRs(dm, samples, groups, coord,
%                                               i_name, i_value, ...)
% -------------------------------------------------------------------
% Description:  detects DMRs according to scanning statistics.
% Input:        {samples} cell array of samples. Each element is either an
%                   ancient sample (of type AMSAMPLE) or a modern sample
%                   (of type MMSAMPLE).
%               {groups} cell array of the group names, each name
%                   corresponding to a single sample.
%               {coord} coordinates of the CpGs in the genome
%                   (GCOORDINATES object).
%               <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'fname' name of file to print diagnostics into (def:
%                       'groupDMRs.txt').
%                   'win_size' window size for smoothing. If 'meth', it is
%                       taken as the value used to reconstruct the
%                       methylation in each sample. If 'auto', a
%                       recommended value is computed for every chromosome
%                       of each sample. Otherwise, it can be a scalar (used
%                       for all chromosomes in all samples), a vector over
%                       the samples (same window size is used for all
%                       chromosomes of each ancient individual, nan is
%                       substituted for each modern individual), or a
%                       cell array with values per individual and
%                       chromosome (def: 'meth').
%                   'win_algorithm' structure with parameters
%                       required to determine window size, see structure
%                       {p_alg} in AMSAMPLE/DETERMINEWINDOWSIZE. It is
%                       assumed that the same algorithm with the same
%                       paramters is applied to all samples.
%                   'lcf' low coverage factor. If 'meth', it is
%                       taken as the value used in reconstructing the
%                       methylation of each sample (def: 'meth').
%                   'delta' minimum methylation difference between the two
%                       groups (def: 0.5).
%                   'min_bases' the minimum length of each DMR in [bases].
%                       Shorter DMRs are filtered out (def: 100).
%                   'drate' deamination rates of each of the ancient
%                       samples. If provided, its length should match the
%                       number of ancient samples. NaN designate the
%                       default of a specific sample (def: taken from the
%                       amSamples objects).
%                   'min_meth' sets a lower bound to the methylation levels
%                       of the reference, such that in every position
%                           ref_meth = max(ref_meth,ref_min_meth)
%                       (def: 0).
%                   'max_meth' sets an upper bound to the methylation
%                       levels of the reference, such that in every
%                       position
%                           ref_meth = min(ref_meth,ref_max_meth)
%                       (def: 1).
%                   'min_Qt' DMRs with Qt < min_Qt are filtered out (def:
%                       0).
%                   'min_CpGs' DMRs whose number of CpGs is less than
%                       min_CpGs are filtered out (def: 10).
%                   'max_adj_dist' max distance between adjacent CpGs
%                       within the same DMR. If the distance between
%                       consecutive CpG positions is larger than
%                       max_adj_dist, the algorithm sets Qt to 0 (def:
%                       1000).
%                   'min_finite' a vector of length {no_groups} stating the
%                       minimum number of ancient samples for which we
%                       require data. If in a position there is not
%                       enough samples with data, a NaN is substituted in
%                       this position. It can also be a fraction between 0
%                       and 1, in which case it is understood as the
%                       minimum fraction of the total number of ancient
%                       samples in the group (def: 1 for each group).
%                   'k' minimum number of standard errors separating the
%                       groups (def: 2).
%                   'max_iterations' maximum number of iteration in the
%                       Newton-Raphson phase (def: 20).
%                   'tol' tolerance in the Newton-Raphson phase (def:
%                       1e-3).
%                   'to_repot' TRUE if reporting to the display is
%                       desired (default = TRUE).
%                   'chromosomes' (def: all chromosomes) names of
%                       chromosomes to be analyzed.
%                   'realORnot' false if the analysis run on simulations
%                   'no_perm' the number of permutation in case of
%                       simulation.
% Output:       {dm} modified DMRs object.
%               {Qt_up} cell array of chromosomes, containing the Qt_up
%                   values along the chromosome.
%               {Qt_down} cell array of chromosomes, containing the Qt_down
%                   values along the chromosome.

% (c) Liran Carmel & Yoav Mathov & Benny Yakir
% Classification: RoAM
% Last revision date: 01-Jan-2019

% error message
msg_id = 'DMRs:groupDMRs';

% parse input line
[p_alg, p_sim, to_report, fid, chromosomes, is_ancient, grp] = ...
    parseInput(samples,groups,coord,varargin{:});

% initializations
no_chr = length(chromosomes);
no_samples = length(samples);
S = length(samples);
iSa = find(is_ancient);
iSm = allbut(iSa,no_samples);
Qt_up = cell(1,no_chr);
Qt_down = cell(1,no_chr);
for chr = 1:no_chr
    Qt_up{chr} = nan(length(...
        coord.coordinates{coord.index(chromosomes{chr})}),1);
    Qt_down{chr} = Qt_up{chr};
end
samp_names = cell(1,no_samples);
species = cell(1,no_samples);
reference = samples{1}.reference;
for ss = 1:no_samples
    samp_names{ss} = samples{ss}.name;
    species{ss} = samples{ss}.species;
    if ~strcmp(samples{ss}.reference,reference)
        error(msg_id,'Sample %s does not use reference %s',...
            samples{ss}.name,reference);
    end
end
obj = DMRs('samples',samp_names,'groups',grp,'is_ancient',is_ancient,...
    'chromosomes',chromosomes,'species',species,'reference',reference);

% split samples between groups
no_groups = nogroups(grp);
if no_groups ~= 2
    error(msg_id,'Currently, the algorithm only works on two groups');
end
giSa = cell(1,no_groups);
giSm = cell(1,no_groups);
gSa = zeros(1,no_groups);
gSm = zeros(1,no_groups);
for gg = 1:no_groups
    giSa{gg} = intersect(iSa,grp.grp2samp(gg));
    giSm{gg} = intersect(iSm,grp.grp2samp(gg));
    gSa(gg) = length(giSa{gg});
    gSm(gg) = length(giSm{gg});
end

% loop on chromosomes
cdm = cDMRs(no_chr);
for chr = 1:no_chr
    % report
    if to_report
        fprintf(1,'Processing chromosome %s\n',chromosomes{chr});
        fprintf(fid,'Processing chromosome %s\n',chromosomes{chr});
    end
    % find matching chromosomes in the different samples
    % {chr_num} is the number of the current chromosome in each sample
    chr_num = zeros(1,S);
    for samp = 1:S
        chr_num(samp) = samples{samp}.indexofchr(chromosomes{chr});
    end
    % for ancient samples, compute tij and nij of windows
    no_positions = length(coord.coordinates{chr});
    % loop on groups
    m = zeros(no_groups,no_positions);  % methylation
    dm = zeros(no_groups,no_positions); % standard error in methylation
    for gg = 1:no_groups
        % compute tij and nij
        tij = zeros(gSa(gg),no_positions);
        nij = zeros(gSa(gg),no_positions);
        for samp  = 1:gSa(gg)
            sampid = giSa{gg}(samp);
            [tij(samp,:), nij(samp,:)] = ...
                samples{sampid}.smooth(chr_num(sampid),...
                p_alg.win_size(sampid,chr));
            % remove regions with particularly low coverage
            lct = findlowcoveragethreshold(nij(samp,:),p_alg.lcf(sampid));
            nij(samp,nij(samp,:)<lct) = nan;
        end
        % estimate methylation
        [m(gg,:), dm(gg,:)] = ancient_Newton_Raphson(...
            p_alg.max_iterations,p_alg.tol,...
            p_alg.drate(giSa{gg}),tij,nij);
        % substitute NaNs when there are too many NaNs in the original data
        finit = isfinite(tij) & isfinite(nij);
        finit = sum(finit,1);
        m(gg,finit < p_alg.min_finite(gg)) = nan;
    end
    % compute the two statistics
    m(m>1) = 1;
    diffi = m(1,:) - m(2,:);
    %%% old method
    %%% lt_up = diffi - p_alg.delta;
    %%% lt_down = -diffi - p_alg.delta;
    %%% end of old method
    idm = sqrt(dm(1,:).^2 + dm(2,:).^2);
    lt_up = diffi ./ idm - p_alg.k;
    lt_up(diffi < p_alg.delta) = 0;
    lt_down = -diffi ./ idm - p_alg.k;
    lt_down(-diffi < p_alg.delta) = 0;
    % remove NaNs
    not_nans = isfinite(lt_up);
    lt_up = lt_up(not_nans);
    lt_down = lt_down(not_nans);
    coordi = coord.coordinates{coord.index(chromosomes{chr})};
    coordi = coordi(not_nans);
    m_approxi = m(:,not_nans);
    no_finite_positions = length(coordi);
    % make index
    trunc2clean = 1:no_positions;
    trunc2clean = trunc2clean(not_nans);
    % compute the distance between adjacent CpG positions in {coordi}
    % mark all positions whose preceding CpG is at most
    % p_DMRs.max_adj_dist far by 1, and all the others by 0.
    coordi_plus1 = [coordi(1)-2; coordi(1:end-1)];
    coordi_diff = coordi - coordi_plus1 - 1;    % distance - 1
    coordi_diff = floor(coordi_diff/p_alg.max_adj_dist);
    coordi_diff = (coordi_diff == 0);
    % initialize {iQt_up}
    iQt_up = zeros(1,no_finite_positions+1);
    % set idx_binary=1 only for positions whose reference methylation
    % in group #2 is below 1-delta (otherwise, clearly there is no DMR)
    idx_binary = zeros(no_finite_positions,1);
    idx_binary(m_approxi(2,:)<1-p_alg.delta) = 1;
    % TODO: test the effect of idx_binary = ones(no_finite_positions,1);
    % compute {iQt_up} recursively
    for pos = 1:no_finite_positions
        iQt_up(pos+1) = max( ...
            coordi_diff(pos) * idx_binary(pos) * ...
            (iQt_up(pos) + lt_up(pos)), 0);
    end
    % initialize {iQt_down}
    iQt_down = zeros(1,no_finite_positions+1);
    % set idx_binary=1 only for positions whose methylation in group #2
    % is above delta
    %TODO: test idx_binary = ones(no_finite_positions,1);
    idx_binary = zeros(no_finite_positions,1);
    idx_binary(m_approxi(2,:)>p_alg.delta) = 1;
    % compute {iQt_down} recursively
    for pos = 1:no_finite_positions
        iQt_down(pos+1) = max( ...
            coordi_diff(pos) * idx_binary(pos) * ...
            (iQt_down(pos) + lt_down(pos)), 0);
    end
    % make a clean version of {Qt}, that takes into account the many CpG
    % positions we had removed earlier. The function reports this vector,
    % which is usful later for plotting
    Qt_up{chr}(not_nans) = iQt_up(2:end);
    Qt_down{chr}(not_nans) = iQt_down(2:end);
    % filter and characterize DMRs
    cdm(chr) = findDMRs(cdm(chr),iQt_up,coordi,p_alg,...
        trunc2clean,m);
    cdm(chr) = findDMRs(cdm(chr),iQt_down,coordi,p_alg,...
        trunc2clean,m);
    cdm(chr) = cdm(chr).set('chromosome',chromosomes{chr});
    % report
    if to_report
        fprintf(1,'\tdetected %d DMRs\n',cdm(chr).noDMRs);
        fprintf(fid,'\tdetected %d DMRs\n',cdm(chr).noDMRs);
    end
end

% substitute fields
p_alg.algorithm = 'groupDMRs';
obj = obj.set('cDMRs',cdm,'algorithm',p_alg);

% close file
if to_report
    fclose(fid);
end

% #########################################################################
function [p_alg,p_sim,to_report,fid,chromosomes,is_ancient,grp] = ...
    parseInput(samples,groups,coord,varargin)

% PARSEINPUT parses input line.
% -------------------------------------------------------------------
% [p_alg, p_sim, to_report, fid, chromosomes, is_ancient, grp] = ...
%                           parseInput(samples,groups,coord,varargin)
% -------------------------------------------------------------------
% Description:  parses the input line.
% Input:        {samples} cell of sample objects.
%               {groups} cell of group names.
%               {coord} GCOORDINATES object.
%               {varargin} original input line.
% Output:       {p_alg} structure with the parameters governing the DMR
%                   detection algorithm.
%               {p_sim} structure with parameters of the simulation.
%               {to_repot} either TRUE or FALSE.
%               {fid} FID of log file.
%               {chromosomes} names of chromosomes to be analyzed.
%               {is_ancient} designates which sample is ancient.
%               {grp} GROUPING object of the samples.

% process
no_samples = length(samples);
fname = 'groupDMRs.txt';
is_ancient = zeros(1,no_samples);
is_ancient(strcmp('amSample',cellfun(@class,samples,...
    'UniformOutput',false))) = 1;

% defaults
to_report = true;
fid = -1;
chromosomes = coord.chr_names;
is_auto_window = false;
is_meth_window = true;
win_size = nan;
lcf = nan * ones(1,no_samples);
is_meth_lcf = true;
p_sim = struct('is_sim',false,'no_simulation',0);
p_alg = struct('delta',0.5,'k',2,'min_bases',100,...
    'min_meth',0,'max_meth',1,'min_Qt',0,'min_CpGs',10,...
    'max_adj_dist',1000,'drate',nan(1,no_samples),'win_size',nan,...
    'min_finite',1,'lcf',nan,'max_iterations',20,'tol',1e-3);
winsize_algorithm = struct;

% read user-specific instructions
for ii = 1:2:(nargin-3)
    switch str2keyword(varargin{ii},7)
        case 'delta  '   % instruction: delta
            p_alg.delta = varargin{ii+1};
        case 'k      '   % instruction: k
            p_alg.k = varargin{ii+1};
        case 'min_bas'   % instruction: min_bases
            p_alg.min_bases = varargin{ii+1};
        case 'drate  '   % instruction: drate
            if length(varargin{ii+1}) == sum(is_ancient)
                p_alg.drate(is_ancient) = varargin{ii+1};
            else
                msg = sprintf('Length of drate');
                msg = sptrinf('%s (%d) does',msg,length(varargin{ii+1}));
                msg = sprintf('%s not match the number of ancient',msg);
                msg = sprintf('%s samples (%d)',msg,sum(is_ancient));
                error(msg_id,msg); %#ok<SPERR>
            end
        case 'to_repo'   % instruction: to_report
            to_report = varargin{ii+1};
        case 'min_met'   % instruction: min_meth
            p_alg.min_meth = varargin{ii+1};
        case 'max_met'   % instruction: max_meth
            p_alg.max_meth = varargin{ii+1};
        case 'min_qt '   % instruction: min_Qt
            p_alg.min_Qt = varargin{ii+1};
        case 'min_cpg'   % instruction: min_CpGs
            p_alg.min_CpGs = varargin{ii+1};
        case 'max_adj'   % instruction: max_adj_dist
            p_alg.max_adj_dist = varargin{ii+1};
        case 'chromos'  % instruction: chromosomes
            chromosomes = varargin{ii+1};
        case 'min_fin'   % instruction: min_finite
            p_alg.min_finite = varargin{ii+1};
        case 'realorn'  % instruction: realornot
            if varargin{ii+1}
                p_sim.is_sim = false;
            else
                p_sim.is_sim = true;
            end
        case 'no_perm'  % instruction: no_permutations
            p_sim.no_simulations = varargin {ii+1};
        case 'fname  '  % instruction: fname
            fname = varargin{ii+1};
        case 'win_siz'  % instruction: win_size
            if ischar(varargin{ii+1})
                switch str2keyword(varargin{ii+1},4)
                    case 'auto'
                        is_auto_window = true;
                        is_meth_window = false;
                    case 'meth'
                        is_auto_window = false;
                        is_meth_window = true;
                end
            else
                is_auto_window = false;
                is_meth_window = false;
                win_size = varargin{ii+1};
            end
        case 'lcf'        % instruction: lcf
            if ischar(varargin{ii+1})
                is_meth_lcf = true;
            else
                is_meth_lcf = false;
                lcf = varargin{ii+1};
            end
        case 'max_ite'        % instruction: max_iterations
            p_alg.max_iterations = varargin{ii+1};
        case 'tol    '        % instruction: tol
            p_alg.tol = varargin{ii+1};
        case 'winsize'    % instruction: winsize_algorithm
            fnames = fieldnames(varargin{ii+1});
            for ff = 1:length(fnames)
                winsize_algorithm.(fnames{ff}) = ...
                    varargin{ii+1}.(fnames{ff});
            end
    end
end

% guarantee methylation values are in the range [0,1]
if p_alg.min_meth > 1
    p_alg.min_meth = 0.01*p_alg.min_meth;
end
if p_alg.max_meth > 1
    p_alg.max_meth = 0.01*p_alg.max_meth;
end

% substitute default drates
for samp = 1:no_samples
    if is_ancient(samp) && isnan(p_alg.drate(samp))
        p_alg.drate(samp) = samples{samp}.drate.rate.global;
    end
end

% process low-coverage-filter
if is_meth_lcf
    % take data from samples
    for samp = 1:no_samples
        if is_ancient(samp)
            lcf(samp) = samples{samp}.methylation.lcf;
        end
    end
end
p_alg.lcf = lcf;

% process window size
no_chr = length(chromosomes);
if ~is_auto_window
    if is_meth_window
        % take data from samples
        win_size = nan * ones(no_samples,no_chr);
        for samp = 1:no_samples
            if is_ancient(samp)
                win_size(samp,:) = samples{samp}.methylation.win_size(...
                    samples{samp}.indexofchr(chromosomes));
            end
        end
    else
        % Option 1: same window size for all individuals/chromosomes
        if length(win_size) == 1
            if ~mod(win_size,2)
                win_size = win_size + 1;
            end
            win_size = win_size * ones(no_samples,no_chr);
            win_size(is_ancient==0,:) = nan;
        % same W for all chromosomes of an individual
        elseif nodims(win_size) == 1
            for samp = 1:no_samples
                if is_ancient(samp)
                    if ~mod(win_size(samp),2)
                        win_size(samp) = win_size(samp) + 1;
                    end
                end
            end
            win_size = win_size(:) * ones(1,no_chr);
            win_size(is_ancient==0,:) = nan;
        % different W for each chromosomes and individual
        else
            for samp = 1:no_samples
                for chr = 1:no_chr
                    if ~mod(win_size(samp,chr),2)
                        win_size(samp,chr) = win_size(samp,chr) + 1;
                    end
                end
            end
            win_size(is_ancient==0,:) = nan;
        end
    end
else        % if window size should be determined automatically
    win_size = zeros(no_samples,no_chr);
    baseline_left = 'win_size(';
    baseline_mid1 = ')= samples{';
    baseline_mid2 = '}.determinewinsize(';
    baseline_right = ');';
    fnames = fieldnames(winsize_algorithm);
    for ii = 1:length(fnames)
        baseline_right = sprintf('''%s'',%s,%s',fnames{ii},...
            winsize_algorithm.(fnames{ii}),baseline_right);
    end
    for samp = 1:no_samples
        if is_ancient(samp)
            for chr = 1:no_chr
                chrid = samples{samp}.indexofchr(chromosomes{chr});
                eval(sprintf('%s%d,%d%s%d%s%d%s',baseline_left,samp,chr,...
                    baseline_mid1,samp,baseline_mid2,chrid,...
                    baseline_right));
            end
        end
    end
    win_size(is_ancient==0,:) = nan;
end
p_alg.win_size = win_size;

% make grouping
[assg, naming] = group(groups);
grp = grouping(assg,naming);
no_groups = nogroups(grp);

% process min_finite
if isscalar(p_alg.min_finite)
    p_alg.min_finite = p_alg.min_finite * ones(1,no_groups);
end
gs = grp.groupsize;
for gg = 1:no_groups
    if p_alg.min_finite(gg) < 1 && p_alg.min_finite(gg) > 0
        p_alg.min_finite(gg) = floor(p_alg.min_finite(gg) * gs(gg));
    end
end

% write to log file - summary of input
if to_report
    fid = fopen(fname,'w');
    ttl = sprintf('Function GROUPDMRS was ran on %s',date);
    fprintf(fid,'%s:\n%s\n',ttl,titunderline(length(ttl)));
    fprintf(fid,'Comparing %d samples from %d groups:\n\n',...
        no_samples,no_groups);
    % make a summary of the input
    max_len = zeros(1,no_groups);
    for gg = 1:no_groups
        max_len(gg) = length(naming{gg});
        for ii = grp.grp2samp(gg)
            max_len(gg) = max(max_len(gg),length(samples{ii}.name));
        end
        max_len(gg) = max_len(gg) + 10;
    end
    % title
    ttl = '|';
    for gg = 1:no_groups
        ttl = sprintf('%s%s|',ttl,fitword(naming{gg},max_len(gg)));
    end
    ul = titunderline(length(ttl));
    fprintf(fid,'\t%s\n\t%s\n',ttl,ul);
    % line by line
    no_rows = max(grp.groupsize);
    for row = 1:no_rows
        ttl = '|';
        idx = nan(1,no_groups);
        for gg = 1:no_groups
            pos = find(assg==gg,1);
            if ~isempty(pos)
                idx(gg) = pos;
            end
            if isnan(idx(gg))
                word = '';
            else
                if is_ancient(idx(gg))
                    word = sprintf('%s (%.2f%%)',samples{idx(gg)}.name,...
                        100*p_alg.drate(idx(gg)));
                else
                    word = sprintf('%s',samples{idx(gg)}.name);
                end
                assg(idx(gg)) = nan;
            end
            ttl = sprintf('%s%s|',ttl,fitword(word,max_len(gg)));
        end
        fprintf(fid,'\t%s\n',ttl);
    end
    fprintf(fid,'\t%s\n\n',ul);
end

% write to log file parameters of the job
if to_report
    % delta
    fprintf(fid,'delta = %.2f\n',p_alg.delta);
    fprintf(fid,'\t[used to compute the lt statistics]\n');
    % k
    fprintf(fid,'k = %.1f\n',p_alg.k);
    fprintf(fid,'\t[used to compute the lt statistics]\n');
    % min_bases
    fprintf(fid,'min_bases = %d\n',p_alg.min_bases);
    fprintf(fid,'\t[minimum length of a DMR (bases). ');
    fprintf(fid,'Shorter DMRs are filtered out]\n');
    % 'max_adj_dist'
    fprintf(fid,'max_adj_dist = %d\n',p_alg.max_adj_dist);
    fprintf(fid,'\t[max distance between adjacent CpGs within the ');
    fprintf(fid,'same DMR (bases)]\n');
    % 'min_Qt'
    fprintf(fid,'min_Qt = %.2f\n',p_alg.min_Qt);
    fprintf(fid,'\t[a DMR must have Qt >= min_Qt]\n');
    % 'min_CpGs'
    fprintf(fid,'min_CpGs = %d\n',p_alg.min_CpGs);
    fprintf(fid,'\t[a DMR must contain at least min_CpGs CpGs]\n');
    % 'min_meth' and 'max_meth'
    fprintf(fid,'min_meth = %.2f\n',p_alg.min_meth);
    fprintf(fid,'\t[lower bound to the methylation levels of the ');
    fprintf(fid,'reference, ref_meth = max(ref_meth,min_meth)]\n');
    fprintf(fid,'max_meth = %.2f\n',p_alg.max_meth);
    fprintf(fid,'\t[upper bound to the methylation levels of the ');
    fprintf(fid,'reference, ref_meth = min(ref_meth,max_meth)]\n');
    % 'min_finite'
    fprintf(fid,'min_finite = [');
    fprintf(fid,'%d ',p_alg.min_finite(1:end-1));
    fprintf(fid,'%d]',p_alg.min_finite(end));
    fprintf(fid,'\t[minimum number of ancient samples per group');
    fprintf(fid,'for which we require data]\n');
    % 'max_iterations'
    fprintf(fid,'max_iterations = ');
    fprintf(fid,'%d',p_alg.max_iterations);
    fprintf(fid,'\t[maximum number of iterations in Newton-Raphson]\n');
    % 'tol'
    fprintf(fid,'tol = ');
    fprintf(fid,'%g',p_alg.tol);
    fprintf(fid,'\t[convergence toleratnce of Newton-Raphson]\n');
    fprintf(fid,'\n');
end

function idm = findDMRs(idm,iQt,icoord,p_alg,trunc2clean,imeth)

% -----------------------------------------------------------
% idm = findDMRs(idm, iQt, icoord, p_alg, trunc2clean, imeth)
% -----------------------------------------------------------
% Description:  Detects DMRs within the Q signals.
% Input:        {idm} DMRs object.
%               {iQt} the Q-vector.
%               {icoord} the coordinates vector.
%               {p_alg} structure of parameters.
%               {trunc2clean} index that maps truncated vectors (e.g.,
%                   icoord(not_nans)) to the full version (e.g., icoord).
%               {imeth} estimated methylation.
% Output:       {idm} updated DMRs object.

% binarize {Qt}
bQt = iQt;
bQt(iQt>0) = 1;

% find 0->1 transitions and 1->0 transitions
dbQt = diff(bQt);
idx0to1 = find(dbQt==1);
idx1to0 = find(dbQt==-1);
% add a last value to cover the case that bQt ends with a run of 1's
idx1to0 = [idx1to0 length(bQt)];

% collect all DMRs whose length is at least {minDMRlen}
for ii = 1:length(idx0to1)
    % beginning and end of a putative DMR (run of Q's)
    beg = idx0to1(ii)+1;
    fin = idx1to0(ii);
    
    % find the precise extent of the DMR
    [maxQt, CpGs_inDMR] = max(iQt(beg:fin));
    
    % compute the true beginning and end of the DMR in the vectors
    % (remember that iQt has an extra first value and is longer by one)
    tbeg = beg - 1;
    tfin = tbeg + CpGs_inDMR - 1;
    
    % check if putative DMR passes our filters. {dlen} is defined from
    % the first position of the first CpG to the second position of the
    % last CpG.
    dlen = icoord(tfin) - icoord(tbeg) + 2;
    if maxQt >= p_alg.min_Qt && dlen >= p_alg.min_bases && ...
            CpGs_inDMR >= p_alg.min_CpGs
        % compute methylation
        CpG_beg = trunc2clean(tbeg);
        CpG_fin = trunc2clean(tfin);
        meth = nanmean(imeth(:,CpG_beg:CpG_fin),2);
        % substitute all in {idm}
        idm = idm.add('gen_begin',icoord(tbeg),...
            'gen_end',icoord(tfin)+1,'no_bases',dlen,...
            'no_CpGs',CpGs_inDMR,'max_Qt',maxQt,...
            'CpG_begin',CpG_beg,'CpG_end',CpG_fin,'methylation',meth);
    end
end

function [m, dm, m0] = ancient_Newton_Raphson(max_iterations,min_tol,...
    pi,tij,nij)

% ANCIENT_NEWTON_RAPHSON solves for m when all samples are ancient
% -----------------------------------------------------------------------
% [m, dm, m0] = ancient_Newton_Raphson(max_iterations, min_tol, pi,
%                                                               tij, nij)
% -----------------------------------------------------------------------
% Description:  solves for m when all samples are ancient
% Input:        {max_iterations} maximum number of iterations.
%               {min_tol} convergence tolerance.
%               {pi} deamination rates of all samples.
%               {tij} Tij as in the algorithm.
%               {nij} Nij as in the algorithm.
% Output:       {m} methylation vector.
%               {dm} standard error of methylation.
%               {m0} initial guess (mainly for debugging purposes).

% compute useful magnitudes
Tj = nansum(tij,1);
no_samples = length(pi);

% transform Nans into 0's before matrix multiplication
idx = union(find(isnan(tij)),find(isnan(nij)));
tij_temp = tij;
tij_temp(idx) = 0;
nij_temp = nij;
nij_temp(idx) = 0;

% Compute more parameters using matrix multiplication
Tpij = pi * tij_temp;
Npij = pi * nij_temp;

% Compute initial guess
m = Tj ./ (Npij - Tpij);
m0 = m;

% make iterations
iter = 1;
tol = 1;
while iter < max_iterations & tol > min_tol %#ok<AND2>
    m_prev = m;
    dldm = Tj ./ m_prev;
    d2ldm2 = - Tj ./ m_prev.^2;
    for samp = 1:no_samples
        dldm = dldm - (nij_temp(samp,:) - tij_temp(samp,:)) * pi(samp) ...
            ./ (1 - m_prev * pi(samp));
        d2ldm2 = d2ldm2 - (nij_temp(samp,:) - tij_temp(samp,:)) * ...
            pi(samp).^2 ./ (1 - m_prev * pi(samp)).^2;
    end
    m = m_prev - dldm ./ d2ldm2;
    tol = abs(m - m_prev) ./ m_prev;
    iter = iter + 1;
end

% compute estimation of the variance
I = Tj ./ m.^2;
for samp = 1:no_samples
    I = I + (nij_temp(samp,:) - tij_temp(samp,:)) * ...
        pi(samp).^2 ./ (1 - m * pi(samp)).^2;
end
dm = sqrt(1./I);