function dm = detectDMRs(~,p_CpGs,varargin)

% DETECTDMRS detects DMRs according to the scanning statistics.
% -------------------------------------------------
% dm = detectDMRs(dm, p_CpGs, i_name, i_value, ...)
% -------------------------------------------------
% Description:  detects DMRs according to scanning statistics.
% Input:        {p_CpGs} vector over chromosomes of structres containing
%                   CpG-site information. Its fields are (not all have to
%                   be filled) 'coordinates','coordinates_per_position',
%                   'No_As','No_Cs','No_Gs','No_Ts','aOga', and 'tOct' (see
%                       COLLECTCPGINFORMATION).
%               <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'methodology' (default: 'wgbs') the type of reference
%                       to which the comparison will be done. 'ancient'
%                       designated another ancient reconstruction. 'wgbs'
%                       designates a measured WGBS map .
%                   'reference' cell array of chromosomes, containing the
%                       reference to which comparison is made. It can
%                       contain beta-values (if 'methodology'='WGBS'), in
%                       which case if they are not on a [0,1] scale, they
%                       are assumed to be on a [0,100] scale. It can also
%                       contain p_CpGs data (if 'methodology'='Ancient').
%                   'min_meth_diff' (default: 0.5) minimum difference
%                       between mean regional methylation that would count
%                       as DMR.
%                   'min_length' (default: 100) the minimum length of each
%                       DMR in [bases]. Shorter DMRs are filtered out.
%                   'drate' a 1-by-2 vector, designating the rate of
%                       deamination in the the sample (1) and in the
%                       reference (2). The second value is relevant only if
%                       'methodology'='ancient'.
%                   'ref_min_meth' (default: 0.39) sets a lower bound to
%                       the methylation levels of the reference, such that
%                       in every position
%                           ref_meth = max(ref_meth,ref_min_meth).
%                   'ref_max_meth' (default: 1) sets an upper boundary to
%                       the methylation levels of the reference, such that
%                       in every position
%                           ref_meth = min(ref_meth,ref_max_meth).
%                   'min_Qt' (default: 0) DMRs with Qt < min_Qt are
%                       filtered out .
%                   'min_CpGs' (default: 10) DMRs whose number of CpGs is
%                       less than min_CpGs are filtered out.
%                   'max_adj_dist' (default: 1000) max distance between
%                       adjacent CpGs within the same DMR. If the distance
%                       between consecutive CpG positions is larger than
%                       max_adj_dist, the algorithm sets Qt to 0.
%                   'to_repot' TRUE if reporting to the display is
%                       desired (default = TRUE).
%                   'upORdown' (default: 'both') either 'up', 'down' or
%                       'both', depending on whether we want to compute
%                       Qt_up or Qt_down.
%                   'l_samples' cell vector with code names of the first
%                       group of samples.
%                   'r_samples' cell vector with code names of the second
%                       group of samples.
% Output:       {DMRs} cell array of chromosomes, containing the DMR
%                   structures for each chromosome.
%               {Qt} cell array of chromosomes, containing the Qt
%                   values along the chromosome, in coordinates that match
%                   {coord}.

% (c) Liran Carmel & David Gokhman & Yoav Mathov
% Classification: RoAM
% Last revision date: 29-Jan-2018


% error message
msg_id = 'DMRs:detectDMRs';

% parse input line
fid = fopen('detectDMRs_log.txt','w');
fprintf(fid,'PARAMETERS OF THE JOB (%s):\n',date);
fprintf(fid,'----------------------------------\n');
[to_report, p_DMRs, l_samples, r_samples] = ...
    parseInput(fid,varargin{:});

% sanity checks
no_chr = length(p_CpGs);
if ~isempty(p_DMRs.ref)
    if length(p_DMRs.ref) ~= no_chr
        error(msg_id,'Expecting %d chromosomes in {ref}',no_chr);
    end
end

% initialize
dm = DMRs(no_chr);
if ~strcmp(p_DMRs.direction,'down')
    Qt_up = cell(1,no_chr);
    for chr = 1:no_chr
        Qt_up{chr} = nan(length(p_CpGs(chr).coordinates),1);
    end
end
if ~strcmp(p_DMRs.direction,'up')
    Qt_down = cell(1,no_chr);
    for chr = 1:no_chr
        Qt_down{chr} = nan(length(p_CpGs(chr).coordinates),1);
    end
end

% report
fprintf(fid,'Processing raw data vectors (%d chromosomes):\n',no_chr);
fprintf(fid,'---------------------------\n');
if to_report
    fprintf(1,'Number of chromosomes: %d\n',no_chr);
end

% loop on chromosomes
for chr = 1:no_chr
%     dm(chr).l_samples = l_samples;
%     dm(chr).r_samples = r_samples;
    % report
    fprintf(fid,'Detecting DMRs in chromosome #%d\n',chr);
    if to_report
        fprintf(1,'Detecting DMRs in chromosome #%d ...\n',chr);
    end
    
    % vectors of current chromosome
    no_ti = p_CpGs(chr).No_Ts;
    no_cti = no_ti + p_CpGs(chr).No_Cs;
    coordi = p_CpGs(chr).coordinates;
    if strcmp(p_DMRs.method,'WGBS')
        refi = p_DMRs.ref{chr};
        % ensure {refi} is a row-vector
        refi = refi(:)';
        % scale {refi} vector to be in the range [0,1]
        if max(refi) > 1
            refi = 0.01 * refi;
        end
        refi(refi<p_DMRs.min_ref) = p_DMRs.min_ref;
        refi(refi>p_DMRs.max_ref) = p_DMRs.max_ref;
    else
        ref_no_ti = p_DMRs.ref(chr).No_Ts;
        ref_no_cti = ref_no_ti + p_DMRs.ref(chr).No_Cs;
    end
    no_CpG_positions = length(no_ti);
    fprintf(fid,'\tTotal number of CpG positions: %s\n',...
        numcommas(no_CpG_positions));
       
    % final summary
    fprintf(fid,'\tIn total, this chromosome includes %s ',...
        numcommas(length(no_ti)));
    fprintf(fid,'CpG positions, of which %s (%.2f%%) are ',...
        numcommas(sum(isnan(no_ti))),sum(isnan(no_ti))/length(no_ti)*100);
    fprintf(fid,'non-informative.\n');
    fprintf(fid,'\tIts effective coverage is %.2f.\n',nanmean(no_cti));
    
    % compute the two statistics lt_up and lt_down
    if strcmp(p_DMRs.method,'WGBS')
        if ~strcmp(p_DMRs.direction,'down')
            % compute lt_up
            delta_up = min(p_DMRs.min_meth_diff, 1 - refi);
            lt_up = (no_ti) .* ...
                log( (refi + delta_up) .* (1 - refi*p_DMRs.drate(1)) ./ ...
                (refi.*(1 - (refi + delta_up)*p_DMRs.drate(1))) ) + ...
                (no_cti) .* ( log(1 - (refi + delta_up)*p_DMRs.drate(1)) - ...
                log(1 - refi*p_DMRs.drate(1)) );
        end
        if ~strcmp(p_DMRs.direction,'up')
            % compute lt_down
            delta_down = min(p_DMRs.min_meth_diff, refi-eps);        
            lt_down = no_ti .* ...
                log( (refi - delta_down) ...
                .* (1 - refi*p_DMRs.drate(1)) ./ ...
                (refi.*(1 - (refi - delta_down)*p_DMRs.drate(1))) ) + ...
                no_cti .* ( log(1 - (refi - delta_down)*p_DMRs.drate(1)) - ...
                log(1 - refi*p_DMRs.drate(1)) );
        end
    else
        % compute ratios to compare
        refi = ref_no_ti ./ ref_no_cti / p_DMRs.drate(2);
        refi(refi>1) = 1;
        sampi = no_ti ./ no_cti / p_DMRs.drate(1);
        sampi(sampi>1) = 1;
        diffi = sampi - refi;
        if ~strcmp(p_DMRs.direction,'down')
            % compute lt_up
            lt_up = diffi - p_DMRs.min_meth_diff;
        end
        if ~strcmp(p_DMRs.direction,'up')
            % compute lt_down
            lt_down = -diffi - p_DMRs.min_meth_diff;
        end
    end
    
    % remove NaNs
    if ~strcmp(p_DMRs.direction,'down')
        not_nans = isfinite(lt_up);
        lt_up = lt_up(not_nans);
    else
        not_nans = isfinite(lt_down);
    end
    if ~strcmp(p_DMRs.direction,'up')
        lt_down = lt_down(not_nans);
    end
    coordi = coordi(not_nans);
    refi = refi(not_nans);
    
    % compute the distance between adjacent CpG positions in {coordi}
    % mark all positions whose preceding CpG is at most
    % p_DMRs.max_adj_dist far by 1, and all the others by 0.
    coordi_plus1 = [coordi(1)-2; coordi(1:end-1)];
    coordi_diff = coordi - coordi_plus1 - 1;    % distance - 1
    coordi_diff = floor(coordi_diff/p_DMRs.max_adj_dist);
    coordi_diff = (coordi_diff == 0);
    
    % compute Qt
    no_finite_positions = length(coordi);
    if ~strcmp(p_DMRs.direction,'down')
        % initialize {iQt_up}
        iQt_up = zeros(1,no_finite_positions+1);
        % set idx_binary=1 only for positions whose reference methylation
        % is below 1-min_meth_diff
        %idx_binary = ones(no_finite_positions,1);
        idx_binary = zeros(no_finite_positions,1);
        idx_binary(refi<1-p_DMRs.min_meth_diff) = 1;
        % compute {iQt_up} recursively
        for pos = 1:no_finite_positions
            iQt_up(pos+1) = max( ...
                coordi_diff(pos) * idx_binary(pos) * ...
                (iQt_up(pos) + lt_up(pos)), 0);
        end
    end
    if ~strcmp(p_DMRs.direction,'up')
        % initialize {iQt_down}
        iQt_down = zeros(1,no_finite_positions+1);
        % set idx_binary=1 only for positions whose reference methylation
        % is above min_meth_diff
        %idx_binary = ones(no_finite_positions,1);
        idx_binary = zeros(no_finite_positions,1);
        idx_binary(refi>p_DMRs.min_meth_diff) = 1;
        % compute {iQt_down} recursively
        for pos = 1:no_finite_positions
            iQt_down(pos+1) = max( ...
                coordi_diff(pos) * idx_binary(pos) * ...
                (iQt_down(pos) + lt_down(pos)), 0);
        end
    end
    
    % make a clean version of {Qt}, that takes into account the many CpG
    % positions we had removed earlier. The function reports this vector,
    % which is usful later for plotting
    if ~strcmp(p_DMRs.direction,'down')
        Qt_up{chr}(not_nans) = iQt_up(2:end);
    end
    if ~strcmp(p_DMRs.direction,'up')
        Qt_down{chr}(not_nans) = iQt_down(2:end);
    end
    
    % filter and characterize DMRs
    switch p_DMRs.direction
        case 'up'
            dm(chr) = findDMRs(dm(chr),iQt_up,coordi,refi,p_DMRs,true);
        case 'down'
            dm(chr) = findDMRs(dm(chr),iQt_down,coordi,refi,p_DMRs,false);
        case 'both'
            dm(chr) = findDMRs(dm(chr),iQt_up,coordi,refi,p_DMRs,true);
            dm(chr) = findDMRs(dm(chr),iQt_down,coordi,refi,p_DMRs,false);
    end
    
    % report
    no_DMRs = dm(chr).no_DMRs;
    switch p_DMRs.direction
        case 'up'
            fprintf(fid,'\tFound %d up-DMRs\n\n',no_DMRs);
            if to_report
                fprintf(1,'\tFound %d up-DMRs\n\n',no_DMRs);
            end
        case 'down'
            fprintf(fid,'\tFound %d down-DMRs\n\n',no_DMRs);
            if to_report
                fprintf(1,'\tFound %d down-DMRs\n\n',no_DMRs);
            end
        case 'both'
            no_uDMRs = sum(dm(chr).is_up);
            no_dDMRs = no_DMRs - no_uDMRs;
            fprintf(fid,'\tFound %d up-DMRs (%.1f%%)\n',...
                no_uDMRs,100*no_uDMRs/no_DMRs);
            fprintf(fid,'\tFound %d down-DMRs (%.1f%%)\n',...
                no_dDMRs,100*no_dDMRs/no_DMRs);
            if to_report
                fprintf(1,'\tFound %d up-DMRs (%.1f%%)\n',...
                    no_uDMRs,100*no_uDMRs/no_DMRs);
                fprintf(1,'\tFound %d down-DMRs (%.1f%%)\n',...
                    no_dDMRs,100*no_dDMRs/no_DMRs);
            end
    end
    % update DMR vector
    dm(chr) = dm(chr).set('chromosome',chr);
    dm(chr) = dm(chr).populateCpGcoords(p_CpGs(chr).coordinates);
end

% close file
fclose(fid);

% #########################################################################
function [to_report, p_DMRs, l_samples, r_samples] = ...
    parseInput(fid, varargin)

% PARSEINPUT parses input line.
% -----------------------------------------------------------------------
% [to_report, p_DMRs, l_samples r_samples] = 
%                                               parseInput(fid, varargin)
% -----------------------------------------------------------------------
% Description:  parses the input line.
% Input:        {fid} FID of log file.
%               {varargin} original input line.
% Output:       {to_repot} either TRUE or FALSE.
%               {p_DMRs} structure with the parameters governing the DMR
%                   detection algorithm.
%               {l_samples} list of the first group of samples.
%               {r_samples} list of the second group of samples.

% defaults
to_report = true;
l_samples = '';
r_samples = '';
p_DMRs = struct('method','wgbs','ref',[],...
    'direction','both','min_meth_diff',0.5,'min_length',100,...
    'min_ref',0.39,'max_ref',1,'min_Qt',0,'min_CpGs',10,...
    'max_adj_dist',1000,'drate',[nan nan]);

% read user-specific instructions
for ii = 1:2:(nargin-1)
    switch str2keyword(varargin{ii},7)
        case 'methodo'   % instruction: methodology
            switch str2keyword(varargin{ii+1},3)
                case 'anc'
                    p_DMRs.method = 'Ancient';
                case 'wgb'
                    p_DMRs.method = 'WGBS';
                otherwise
                    error(sprints('Unknown methology ''%s''',...
                        varargin{ii+1}));
            end
        case 'referen'   % instruction: reference
            p_DMRs.ref = varargin{ii+1};
        case 'min_met'   % instruction: min_meth_diff
            p_DMRs.min_meth_diff = varargin{ii+1};
        case 'min_len'   % instruction: min_length
            p_DMRs.min_length = varargin{ii+1};
        case 'drate  '   % instruction: drate
            if length(varargin{ii+1}) == 1
                p_DMRs.drate(1) = varargin{ii+1};
            else
                p_DMRs.drate = varargin{ii+1};
            end
        case 'to_repo'   % instruction: to_report
            to_report = varargin{ii+1};
        case 'ref_min'   % instruction: ref_min_meth
            p_DMRs.min_ref = varargin{ii+1};
            if p_DMRs.min_ref > 1
                p_DMRs.min_ref = 0.01*p_DMRs.min_ref;
            end
        case 'ref_max'   % instruction: ref_max_meth
            p_DMRs.max_ref = varargin{ii+1};
            if p_DMRs.max_ref > 1
                p_DMRs.max_ref = 0.01*p_DMRs.max_ref;
            end
        case 'min_qt '   % instruction: min_Qt
            p_DMRs.min_Qt = varargin{ii+1};
        case 'mincpgs'   % instruction: minCpGs
            p_DMRs.min_CpGs = varargin{ii+1};
        case 'max_adj'   % instruction: max_adj_dist
            p_DMRs.max_adj_dist = varargin{ii+1};
        case 'upord  '   % instruction: upORdown
            switch str2keyword(varargin{ii+1},2)
                case 'up'
                    p_DMRs.direction = 'up';
                case 'do'
                    p_DMRs.direction = 'down';
                case 'bo'
                    p_DMRs.direction = 'both';
                otherwise
                    error(sprints(...
                        'Unknown value for ''upORdown'': ''%s''',...
                        varargin{ii+1}));
            end
        case 'l_sampl'    % instruction: l_samples
            l_samples = varargin{ii+1};
        case 'r_sampl'    % instruction: r_samples
            r_samples = varargin{ii+1};
    end
end

% write to log file
% methodology and drate
if ~isempty(l_samples)
    fprintf(fid,'Comparing (%s) to (%s)\n',l_samples{1},r_samples{1});
end
if strcmp(p_DMRs.method,'Ancient')
    fprintf(fid,'Methodology = ''Ancient''\n');
    fprintf(fid,'\t[compare sample to ancient reconstruction]\n');
    fprintf(fid,'Deamination rates are %.3f (%s), ',...
        p_DMRs.drate(1),l_samples{1});
    fprintf(fid,'and %.3f (%s)\n',...
        p_DMRs.drate(2),r_samples{1});
else
    fprintf(fid,'Methodology = ''WGBS''\n');
    fprintf(fid,'\t[compare sample to measured reference]\n');
    fprintf(fid,'Deamination rate of the sample (%s) is %.3f\n',...
        l_samples{1},p_DMRs.drate(1));
end
% min_meth_diff
fprintf(fid,'min_meth_diff = %.2f\n',p_DMRs.min_meth_diff);
fprintf(fid,'\t[minimum difference between mean regional methylation ');
fprintf(fid,'that would count as DMR]\n');
% min_length
fprintf(fid,'min_length = %d\n',p_DMRs.min_length);
fprintf(fid,'\t[minimum length of a DMR [bases]. ');
fprintf(fid,'Shorter DMRs are filtered out]\n');
% 'ref_min_meth' and 'ref_max_meth'
if strcmp(p_DMRs.method,'WGBS')
    fprintf(fid,'ref_min_meth = %.2f\n',p_DMRs.min_ref);
    fprintf(fid,'\t[lower bound to the methylation levels of the ');
    fprintf(fid,'reference, ref_meth = max(ref_meth,ref_min_meth).\n');
    fprintf(fid,'ref_max_meth = %.2f\n',p_DMRs.max_ref);
    fprintf(fid,'\t[upper bound to the methylation levels of the ');
    fprintf(fid,'reference, ref_meth = min(ref_meth,ref_max_meth).\n');
end
% 'min_Qt'
fprintf(fid,'min_Qt = %.2f\n',p_DMRs.min_Qt);
fprintf(fid,'\t[a DMR must have Qt >= min_Qt]\n');
% 'min_CpGs'
fprintf(fid,'min_CpGs = %d\n',p_DMRs.min_CpGs);
fprintf(fid,'\t[a DMR must contain at least min_CpGs CpGs]\n');
% 'max_adj_dist'
fprintf(fid,'max_adj_dist = %d\n',p_DMRs.max_adj_dist);
fprintf(fid,'\t[max distance between adjacent CpGs within the same DMR ');
fprintf(fid,'[bases]]\n');
% direction
if strcmp(p_DMRs.direction,'both')
    fprintf(fid,'Detecting both up- and down-DMRs\n');
elseif strcmp(p_DMRs.direction,'up')
    fprintf(fid,'Detecting up-DMRs\t(the sample is hypermethylated)\n');
else
    fprintf(fid,'Detecting down-DMRs\t(the sample is hypomethylated)\n');
end
fprintf(fid,'\n');

function dm_chr = findDMRs(dm_chr,iQt,icoord,iref,p_DMRs,is_up)

% --------------------------------------------
% dm_chr = findDMRs(dm_chr,Qt,coord,ref,p_DMRs,is_up)
% --------------------------------------------
% Description:  Detects DMRs within the Q signals.
% Input:        {dm_chr} DMRs object.
%               {Qt} the Q-vector.
%               {coord} the coordinates vector.
%               {ref} reference methylation levels.
%               {p_DMRs} structure of parameters.
%               {is_up} determines whether these are up- or down- DMRs.
% Output:       {dm_chr} updated DMRs object.

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
    [maxQt, imax_coor] = max(iQt(beg:fin));
    
    % check if putative DMR passes our filters. {dlen} is defined from
    % the first position of the first CpG to the second position of the
    % last CpG.
    dlen = icoord(beg+imax_coor-2) - icoord(beg-1) + 2;
    CpGs_inDMR = sum(~isnan(iref((beg-1):(beg+imax_coor-2))));  % TODO: not sure that nans are possible in the first place
    if maxQt >= p_DMRs.min_Qt && dlen >= p_DMRs.min_length && ...
            CpGs_inDMR >= p_DMRs.min_CpGs
        dm_chr = dm_chr.add('gen_begin',icoord(beg-1),...
            'gen_end',icoord(beg+imax_coor-2)+1,'no_bases',dlen,...
            'no_CpGs',CpGs_inDMR,'max_Qt',maxQt,'is_up',is_up,...
            'is_present',true,'is_computed',true,'is_missing',false);
    end
end