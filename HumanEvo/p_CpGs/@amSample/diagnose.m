function obj = diagnose(obj,varargin)

% DIAGNOSE computes basic statistics and recommends thresholds values.
% --------------------------------------------
% obj = diagnose(p_CpGs, i_name, i_value, ...)
% --------------------------------------------
% Description:  computes basic statistics on each input chromosome, and
%               recommends what thresholds to use when excluding PCR
%               duplicates and true mutations.
% Input:        <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'fname' (def: 'NAME_diagnostics.txt', where NAME is
%                       obj.name) name of file to print diagnostics into.
%                   'span' (def: 5) parameter for outlier removal in the
%                       vector {no_ct}.
%                   'strict' (def: true) TRUE/FALSE determining whether to
%                       apply a more strict threshold on {no_ct} on top of
%                       the looser one computed by the simple outlier
%                       removal procedure.
%                   'tol' (def: 1e-3) tolerance of BMM model convergence.
%                   'to_compare' (def: true) binary variable, indicating
%                       whether to compare DIANGNOSE's procedure of outlier
%                       removal to that of simple thresholds.
%                   'max_tOct' (def: 0.25) simple threshold, removing all
%                       positions where no_t/no_ct > max_tOct. Performance
%                       of DIAGNOSE is compared to this filter.
%                   'max_aOga' (def: 0.25) simple threshold applicable in
%                       double-stranded libraries, removing all positions 
%                       where no_a/no_ga > max_aOga. Performance of
%                       DIAGNOSE is compared to this filter. Usually, if
%                       the libray is single-stranded, set
%                           g2a_thresh = 2 * (1/(coverage/2)).
%                       otherwise, take the value of 'c2t_thresh'.
%                   'max_coverage' (def: 100) simple threshold determining
%                       the threhold for removing PCR duplicates.
%                       Performance of DIAGNOSE is compared to this filter.
% Output:       {obj} an amSample with an updated DIAGNOSTICS and P_FILTERS
%                   structures, which include the following fields:
%                   DIAGNOSTICS:
%                       'effective_coverage' vector over the chromosomes,
%                           with the effective coverage of each chromosome
%                           (after removal of PCR duplicates and
%                           mutations). 
%                   P_FILTERS:
%                       'method' of filtering, either 'both', 'tOct', or
%                           'aOga'. For now, this field is left empty, and
%                           is filled only when positions are filtered (see
%                           AMSAMPLE/FILTER).
%                       'max_coverage' vector over chromosomes, with the
%                           threshold coverage value for each chromosome
%                           (see AMSAMPLE/FILTER).
%                       'max_TsPerCoverage' cell array over chromosomes,
%                           with a vector of threshold max_TsPerCoverage
%                           values for each coverage level, in each
%                           chromosome, from 1 to {max_coverage} (see
%                           AMSAMPLE/FILTER).
%                       'max_aOga' maximum aOga per chromosome (see
%                           AMSAMPLE/FILTER).
%                       'max_No_As' maximum As per chromosome (see
%                           AMSAMPLE/FILTER).
%                       'is_filtered' true/false whether the data are
%                           actually filtered or not.

% (c) Liran Carmel
% Classification: RoAM
% Last revision date: 19-Sep-2018

% hard-coded parameters (should consider to include in the input line)
LOW_COVERAGE = 1;

% parse input line
[p_alg, p_compare] = parseInput(obj.name,varargin{:});

% initialize
no_chr = obj.no_chrs;
eff_coverage = zeros(1,no_chr);
coverage_thresh = zeros(1,no_chr);
C2T_thresh = cell(1,no_chr);
fid = fopen(p_alg.fname,'w');

% compute number of rows and columns in output figures
no_rows = floor(sqrt(no_chr));
no_cols = ceil(no_chr/no_rows);
fig_coverage = figure('name',obj.name);
fig_perc_removed_per_coverage = figure('name',obj.name);
fig_perc_removed = figure('name',obj.name);

% Loop on chromosomes
for chr = 1:no_chr   % 1:no_chr
    % report
    if isempty(obj.chr_names)
        fprintf(1,'Diagnozing chromosome #%d\n',chr);
        fprintf(fid,'Diagnozing chromosome #%d\n',chr);
    else
        fprintf(1,'Diagnozing %s\n',obj.chr_names{chr});
        fprintf(fid,'Diagnozing %s\n',obj.chr_names{chr});
    end
    
    % vectors of current chromosome
    no_ti = obj.getNo_Ts(chr);
    no_cti = no_ti + obj.getNo_Cs(chr);
    if strcmp(obj.library,'SS')
        no_ai = obj.getNo_As(chr);
    end
    if strcmp(obj.library,'SS')
        aOgai = obj.getaOga(chr);
        if isempty(aOgai)
            no_gi = obj.getNo_Gs(chr);
            if ~isempty(no_gi)
                aOgai = no_ai ./ (no_ai + no_gi);
            end
        end
    end
    fprintf(fid,'\tNumber of CpG positions: %s\n',...
        numcommas(length(no_cti)));
    fprintf(fid,'\tInitial coverage: %.1f\n',mean(no_cti));

    % compute the upper threshold of {nt} in a few steps
    % (i) crude outlier removal using very loose criterion
    prctls = prctile(no_cti,[25 50 75]); 
    fprintf(fid,'\t{C+T}: median = %d\n',prctls(2));
    fprintf(fid,'\t    [25th 75th] percentiles = [%d %d] ',prctls([1 3]));
    iqrange = prctls(3) - prctls(1);
    thresh = prctls(3) + p_alg.span*iqrange;
    fprintf(fid,'=> initial threshold was set to %d\n',thresh);
    to_remove = find(no_cti > thresh);
    fprintf(fid,'\t    %s positions (%.2f%%) removed as no_ct > %d\n',...
        numcommas(length(to_remove)),...
        100*length(to_remove)/length(no_cti),thresh);
    no_cti(to_remove) = nan;
    % (ii) evaluate the parameters of the normal distribution of {nt}
    nt_zeros = find(no_cti==0);
    nz_cti = no_cti;
    nz_cti(nt_zeros) = nan;
    N = hist(nz_cti,1:thresh);
    more_to_remove = [];
    if p_alg.strict
        [~, imx] = max(N);
        points = (imx-1):(imx+1);
        p2 = polyfit(points,N(points),2);
        mu = -0.5*p2(2)/p2(1);
        N0 = -0.25*p2(2)^2/p2(1) + p2(3);
        imu = ceil(mu);
        idx = (imu-1) + find(N(imu:end)<0.1*N0,1);
        delta = idx - mu;
        f = N(idx)/N0;
        sig = delta / sqrt(2*log(1/f));
        fprintf(fid,'\t    the distribution of {nt} is approximated as ');
        fprintf(fid,'N(mu,sig)=N(%.2f,%.2f)\n',mu,sig);
        % (iii) take the threshold as the first value where the expectation is
        % to get less than one count
        total_counts = sum(N);
        thresh = ceil(mu + sig*sqrt(2*log(total_counts/sig/sqrt(2*pi))));
        fprintf(fid,'\t    Final threshold was set to %d\n',thresh);
        more_to_remove = find(no_cti>thresh);
        fprintf(fid,'\t    %s positions (%.2f%%) removed ',...
            numcommas(length(more_to_remove)),...
            100*length(more_to_remove)/length(no_cti));
        fprintf(fid,'as nt > %d\n',thresh);
    end
    coverage_thresh(chr) = thresh;
    
    % plot the coverage 
    figure(fig_coverage);
    subplot(no_rows,no_cols,chr);
    bar(1:length(N),N);
    line([thresh thresh],get(gca,'ylim'),'color','k');
    if isempty(obj.chr_names)
        xlabel(sprintf('Chromosome #%d\n',chr));
    else
        xlabel(sprintf('%s\n',obj.chr_names{chr}));
    end
    
    % remove outliers
    to_remove = [to_remove; more_to_remove]; %#ok<AGROW>
    no_cti(to_remove) = nan;
    no_ti(to_remove) = nan;
    
    % effective coverage
    eff_coverage(chr) = nanmean(no_cti);
    fprintf(fid,'\t    Effective coverage: %.1f\n',eff_coverage(chr));
    
    % compare to previous PCR-duplicate threshold
    prev_to_remove = find(no_cti>p_compare.max_coverage);
    fprintf(fid,'\t    COMPARISON: %s positions (%.2f%%) would ',...
        numcommas(length(prev_to_remove)),...
        100*length(prev_to_remove)/length(no_cti));
    fprintf(fid,'have been removed if the PCR-duplicate threshold ');
    fprintf(fid,'was %d\n',p_compare.max_coverage);
    
    % report on the zero bin
    fprintf(fid,'\t%s positions (%.2f%%) have no_ct = 0\n',...
        numcommas(length(nt_zeros)),100*length(nt_zeros)/length(no_cti));
    
    % analyze each coverage level independently
    th_C2G = zeros(1,thresh);
    chr_cov = LOW_COVERAGE:thresh;
    no_positions_total = zeros(1,thresh);
    no_positions_removed = zeros(1,thresh);
    for cover = thresh:-1:LOW_COVERAGE
        fprintf(fid,'\tAnalyzing {xt} coming from coverage %d:\n',cover);
        idx = find(no_cti==cover);
        no_positions_cover = length(idx);
        no_positions_total(cover) = no_positions_cover;
        fprintf(fid,'\t\tTotal number of positions: %s\n',...
            numcommas(no_positions_cover));
        H = hist(no_ti(idx),0:cover);
        [p, w] = bmm(H,[0.01 0.5 0.99],[0.9 0.05 0.05],p_alg.tol,[0 1 0]);
        fprintf(fid,'\t\tEstimated homozygous mutations: ');
        fprintf(fid,'Pr(meth) = %.2f; Weight = %.3f\n',p(3),w(3));
        fprintf(fid,'\t\tEstimated heterozygous mutations: ');
        fprintf(fid,'Pr(meth) = %.2f; Weight = %.3f\n',p(2),w(2));
        fprintf(fid,'\t\tEstimated deaminated positions: ');
        fprintf(fid,'Pr(meth) = %.2f; Weight = %.3f\n',p(1),w(1));
        thresh = ceil( (log(w(2)/w(1)) + cover*log( 1/2/(1-p(1)) )) / ...
            log( p(1)/(1 - p(1)) )) - 1;
        th_C2G(cover) = thresh;
        more_to_remove = idx(no_ti(idx) > thresh);
        no_positions_removed(cover) = length(more_to_remove);
        fprintf(fid,'\t\t%s positions (%.2f%%) removed ',...
            numcommas(length(more_to_remove)),...
            100*length(more_to_remove)/length(no_cti));
        fprintf(fid,'as xt > %d ',thresh);
        fprintf(fid,'(corresponding to C2T-ratio of %.3f)\n',...
            (thresh+1)/cover);
        fprintf(fid,'\t\tThe threshold (%d) is expected to remove ',...
            thresh);
        fprintf(fid,'%.1f true positives\n',...
            w(1)*no_positions_cover*binocdf(thresh,cover,p(1),'upper'));
        fprintf(fid,'\t\t\tand to include %.1f false positives\n',...
            w(2)*no_positions_cover*binocdf(thresh-1,cover,0.5) + ...
            w(3)*no_positions_cover*binocdf(thresh-1,cover,p(3)));
        if thresh == cover
            % see what happens without a threshold at all
            fprintf(fid,'\t\tIf no threshold is used (no positions ');
            fprintf(fid,'removed), then %.1f false positives will be ',...
                (w(2)+w(3))*no_positions_cover);
            fprintf(fid,'retained\n');
        end
        
        % computing positions removed using G->A and C->T
        if p_compare.to_compare
            if isfinite(p_compare.max_tOct)
                c2t_crit = idx(no_ti(idx) > 1 & ...
                    no_ti(idx) ./ no_cti(idx) >= p_compare.max_tOct);
                fprintf(fid,'\t\tCOMPARISON: Using C->T ');
                fprintf(fid,'%s positions ',numcommas(length(c2t_crit)));
                fprintf(fid,'are removed.\n');
            end
            if strcmp(obj.library,'SS') && ~isempty(no_ai)
                g2a_crit = idx( (no_ai(idx) == 1 & ...
                    aOgai(idx) >= p_compare.max_aOga) | no_ai(idx) > 1);
                fprintf(fid,'\t\tCOMPARISON: Using G->A %s positions ',...
                    numcommas(length(g2a_crit)));
                fprintf(fid,'are removed.');
                if isfinite(p_compare.max_tOct)
                    g2a_additional = setdiff(g2a_crit,c2t_crit);
                    fprintf(fid,' Of them, %s ',...
                        numcommas(length(g2a_additional)));
                    fprintf(fid,'positions were not removed by the ');
                    fprintf(fid,'C->T ratio threshold.\n');
                    in_both = intersect([g2a_crit; c2t_crit],...
                        more_to_remove);
                    fprintf(fid,'\t\tCOMPARISON: ');
                    if length(in_both) == length(more_to_remove)
                        fprintf(fid,'All positions removed only by ');
                        fprintf(fid,'looking at {xt} (%s) are ',...
                            numcommas(length(more_to_remove)));
                        fprintf(fid,'also removed when looking ');
                        fprintf(fid,'at both C->T and G->A.\n');
                    else
                        fprintf(fid,'Of the %s positions ',...
                            numcommas(length(more_to_remove)));
                        fprintf(fid,'removed by only looking at {xt}, ');
                        fprintf(fid,'%s (%.1f%%) are also removed ',...
                            numcommas(length(in_both)),...
                            100*length(in_both)/length(more_to_remove));
                        fprintf(fid,'when looking at both C->T and ');
                        fprintf(fid,'G->A.\n');
                    end
                else
                    fprintf(fid,'\n');
                end
                in_crit = intersect(g2a_crit,more_to_remove);
                fprintf(fid,'\t\tCOMPARISON: ');
                if length(in_crit) == length(more_to_remove)
                    fprintf(fid,'All positions removed only by looking ');
                    fprintf(fid,'at {xt} (%s) are ',...
                        numcommas(length(more_to_remove)));
                    fprintf(fid,'removed by looking at G->A.\n');
                else
                    fprintf(fid,'Of the %s positions ',...
                        numcommas(length(more_to_remove)));
                    fprintf(fid,'removed by only looking at {xt}, ');
                    fprintf(fid,'%s (%.1f%%) are also removed when ',...
                        numcommas(length(in_crit)),...
                        100*length(in_crit)/length(more_to_remove));
                    fprintf(fid,'looking at G->A.\n');
                end
            end
        end
        % removing the positions
        no_cti(more_to_remove) = nan;
        no_ti(more_to_remove) = nan;
    end
    C2T_thresh{chr} = th_C2G;
    
    % plot the ratio of removed to total per coverage
    quotient = no_positions_removed ./ no_positions_total * 100;
    figure(fig_perc_removed_per_coverage);
    subplot(no_rows,no_cols,chr);
    plot(chr_cov,quotient,'bp');
    if isempty(obj.chr_names)
        xlabel(sprintf('Chromosome #%d\n',chr));
    else
        xlabel(sprintf('%s\n',obj.chr_names{chr}));
    end
    ylabel('%removed(coverage)/total(coverage)');
    grid on;
    
    % plot the ratio of removed to total per coverage
    tot = sum(no_positions_total);
    figure(fig_perc_removed);
    subplot(no_rows,no_cols,chr);
    plot(chr_cov,no_positions_removed/tot*100,'bp');
    if isempty(obj.chr_names)
        xlabel(sprintf('Chromosome #%d\n',chr));
    else
        xlabel(sprintf('%s\n',obj.chr_names{chr}));
    end
    ylabel('%removed(coverage)/total');
    grid on;
end
    
% close file
fclose(fid);

% substitute in returned variables
obj.diagnostics = struct('effective_coverage',eff_coverage);
obj.p_filters = struct('method','','max_coverage',coverage_thresh,...
    'max_TsPerCoverage',{C2T_thresh},'max_aOga',[],'max_No_As',[],...
    'is_filtered',false);

% #########################################################################
function [p_alg, p_compare] = parseInput(name,varargin)

% PARSEINPUT parses input line.
% -----------------------------------------------
% [p_alg, p_compare] = parseInput(name, varargin)
% -----------------------------------------------
% Description:  parses the input line.
% Input:        {name} obj.name.
%               {varargin} original input line.
% Output:       {p_alg} parameters of the algorithm. These include the
%                   fields: 'span', 'tol'.
%               {p_compare} parameters to compare the outlier removal
%                   recommendations of DIAGNOSE with simpler filters. These
%                   include the fields: 'to_compare', 'max_tOct',
%                   'max_aOga', 'max_coverage'.

% defaults
p_compare = struct('to_compare',true,...
    'max_tOct',0.25,'max_aOga',0.25,'max_coverage',100);
p_alg = struct('fname',sprintf('%s_diagnostics.txt',name),...
    'span',5,'strict',true,'tol',1e-3);

% read user-specific instructions
for ii = 1:2:(nargin-1)
    switch str2keyword(varargin{ii},5)
        case 'to_co'   % instruction: to_compare
            p_compare.to_compare = varargin{ii+1};
        case 'max_t'   % instruction: max_tOct
            thresh = varargin{ii+1};
            if thresh > 1
                thresh = 0.01 * thresh;
            end
            p_compare.max_tOct = thresh;
        case 'max_a'   % instruction: max_aOga
            thresh = varargin{ii+1};
            if thresh > 1
                thresh = 0.01 * thresh;
            end
            p_compare.max_aOga = thresh;
        case 'max_c'   % instruction: max_coverage
            p_compare.max_coverage = varargin{ii+1};
        case 'span '   % instruction: span
            p_alg.span = varargin{ii+1};
        case 'tol  '   % instruction: tol
            p_alg.tol = varargin{ii+1};
        case 'stric'   % instruction: strict
            p_alg.strict = varargin{ii+1};
        case 'fname'   % instruction: fname
            p_alg.fname = varargin{ii+1};
    end
end