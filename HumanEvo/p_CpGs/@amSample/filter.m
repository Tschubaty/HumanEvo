function f_ams = filter(ams,varargin)

% FILTER removes unreliable CpG positions.
% ----------------------------------------
% f_ams = ams.filter(i_name, i_value, ...)
% ----------------------------------------
% Description:  removes information from CpG sites that did not pass
%               various quality control tests.
% Input:        <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'max_coverage' (def: ams.p_filters.max_coverage. If
%                       empty, than 100) a coverage threshold used to 
%                       remove PCR duplicates. All positions for which
%                       No_Cs > {max_coverage} are removed. This
%                       parameter can be one of the following:
%                       (1) scalar, in which case it is used for all
%                           chromosomes.
%                       (2)vector over chromosomes, with a different
%                           threshold to each chromosome.
%                   'method' (def: 'both' for library='SS', 'tOct' for
%                       library='DS') method used to remove true 
%                       C->T mutations. It can be one of the following:
%                       (1) 'tOct' uses information from No_Ts and No_Cs
%                           only. Uses the parameter 'max_tOct' below.
%                       (2) 'aOga' uses information from No_As and No_Gs
%                           only. Uses the parameters 'max_aOga' and
%                           'max_No_As' below. Applicable only when the
%                           library field is 'SS'.
%                       (3) 'both' combines both options, using information
%                           from No_Ts, No_Cs, No_As, and No_Gs. Uses the
%                           parameters 'max_tOct', 'max_aOga' and
%                           'max_No_As' below. Applicable only when the
%                           library field is 'SS'.
%                   'max_TsPerCoverage' (def:
%                       ams.p_filters.max_TsPerCoverage. If empty, than
%                       computed from 'max_tOct') a threshold used to
%                       identify sites with a true C->T mutation. All
%                       positions with a certain coverage No_Ts+No_Cs=C,
%                       where No_Ts > max_TsPerCoverage(C) are removed.
%                       This parameter is a cell array over chromosomes,
%                       with a vector of max_TsPerCoverage for each
%                       coverage level in each chromosome.
%                   'max_tOct' (def: 0.25) a threshold used to identify
%                       sites with a true C->T mutation. All positions 
%                       where tOct >= max_tOct are removed. This
%                       parameter can be one of the following:
%                       (1) scalar, in which case it is used for all
%                           chromosomes, and all coverage levels.
%                       (2) vector over chromosomes, with a different
%                           threshold to each chromosome. This threshold is
%                           used for all coverage levels.
%                       (3) cell array over chromosomes, with a different
%                           threshold to each chromosome and each coverage
%                           level.
%                   'max_aOga' (def: 0.25) a threshold used to identify
%                       sites with a true C->T mutation. Applicable only
%                       when the library field is 'SS'. All positions 
%                       where aOga >= max_aOga are removed. A rule of
%                       thumb to determine it is to set
%                           max_aOga = 2 * (1/(coverage/2)).
%                       otherwise, take the value of 'max_tOct'. This
%                       parameter can be one of the following:
%                       (1) scalar, in which case it is used for all
%                           chromosomes, and all coverage levels.
%                       (2)vector over chromosomes, with a different
%                           threshold to each chromosome.
%                   'max_No_As' (def: 1) all positions where No_As >
%                       {max_No_As} are removed. Applicable only
%                       when the library field is 'SS'.
%                   'min_No_Ts' (def: 1) the filter 'max_tOct' is
%                       applied only to positions where No_Ts >
%                       {min_No_Ts}.
%                   'merge' (def: true) whether to merge the two
%                       consecutive coordinates of every CpG position.
%                   'fname' (def: 'NAME_filter.txt', where NAME is
%                       ams.name) log file name.
% Output:       {f_ams} a modified AMSAMPLE object, where sites that had
%                   been removed are marked with NaNs, and where fields
%                   that are not further required ('No_As', 'No_Gs',
%                   'aOga', and 'tOct') are emptied.

% (c) Liran Carmel & David Gokhman & Yoav Mathov
% Classification: RoAM
% Last revision date: 20-Jan-2019

% parse input line and initialize
no_chr = ams.no_chrs;
[p_alg, to_merge, fname] = parseInput(ams,varargin{:});
f_ams = ams;
f_ams.p_filters = p_alg;
fid = fopen(fname,'w');

% Loop on chromosomes
total_removed = 0;
for chr = 1:no_chr
    % report
    if isempty(ams.chr_names)
        fprintf(1,'Filtering chromosome #%d\n',chr);
        fprintf(fid,'Filtering chromosome #%d\n',chr);
    else
        fprintf(1,'Filtering %s\n',ams.chr_names{chr});
        fprintf(fid,'Filtering %s\n',ams.chr_names{chr});
    end
    
    % vectors of current chromosome
    no_ti = ams.getNo_Ts(chr);
    no_cti = no_ti + ams.getNo_Cs(chr);
    if ~strcmp(p_alg.method,'tOct')
        no_ai = ams.getNo_As(chr);
        aOgai = ams.getaOga(chr);
    end
    no_positions = length(no_ti);
    fprintf(1,'\tNumber of CpG positions ');
    fprintf(1,'(%d coordinates per ',ams.coord_per_position(chr));
    fprintf(1,'position): %s\n',numcommas(no_positions));
    fprintf(fid,'\tNumber of CpG positions ');
    fprintf(fid,'(%d coordinates per ',ams.coord_per_position(chr));
    fprintf(fid,'position): %s\n',numcommas(no_positions));

    % remove positions whose coverage is too high
    to_remove = find(no_cti > p_alg.max_coverage(chr));
    no_removed = length(to_remove);
    fprintf(fid,'\t%s positions (%.2f%%) removed as No_CTs > %d\n',...
        numcommas(no_removed),100*no_removed/no_positions,...
        p_alg.max_coverage(chr));
    
    % loop on each coverage level and remove positions with high No_Ts
    if ~strcmp(p_alg.method,'aOga')
        % loop over all coverage levels
        for cover = p_alg.max_coverage(chr):-1:1
            idx = find(no_cti==cover);
            fprintf(fid,'\tcoverage %d (%s positions): ',...
                cover,numcommas(length(idx)));
            more_to_remove = ...
                no_ti(idx) > p_alg.max_TsPerCoverage{chr}(cover);
            more_to_remove = extend_removed_positions(idx(more_to_remove)); 
            no_removed = length(to_remove);
            to_remove = unique([to_remove; more_to_remove]); 
            no_removed = length(to_remove) - no_removed;
            fprintf(fid,'\t%s (extended) positions ',...
                numcommas(no_removed));
            fprintf(fid,'were removed as No_Ts > %d\n',...
                p_alg.max_TsPerCoverage{chr}(cover));
        end
    end
    
    % remove more positions if data on A's and G's is available
    if ~strcmp(p_alg.method,'tOct')
        more_to_remove = find( (no_ai <= p_alg.max_No_As(chr) & ...
                aOgai >= p_alg.max_aOga(chr)) | ...
                no_ai > p_alg.max_No_As(chr) );
        more_to_remove = extend_removed_positions(more_to_remove);
        no_removed = length(to_remove);
        to_remove = unique([to_remove; more_to_remove]);
        no_removed = length(to_remove) - no_removed;
        fprintf(fid,'\t%s additional positions ',numcommas(no_removed));
        fprintf(fid,'removed as (1) No_As > %d or (2) No_As <= %d ',...
            p_alg.max_No_As(chr),p_alg.max_No_As(chr));
        fprintf(fid,' and aOga >= %.2f\n',p_alg.max_aOga(chr));
    end
    
    % remove positions and keep only relevant vectors in p_CpGs
    no_removed = length(to_remove);
    fprintf(1,'\tOverall %s positions (%.2f%%) were removed\n',...
        numcommas(no_removed),100*no_removed/no_positions);
    fprintf(fid,'\tOverall %s positions (%.2f%%) were removed\n',...
        numcommas(no_removed),100*no_removed/no_positions);
    no_ti(to_remove) = nan;
    no_cti(to_remove) = nan;
    no_ci = no_cti - no_ti;
    
    % merge positions
    if to_merge
        if f_ams.coord_per_position(chr) == 1
            if isempty(f_ams.chr_names)
                fprintf(1,'Chromosome #%d had already gone ',chr);
                fprintf(1,'through merger\n');
            else
                fprintf(1,'%s had already gone through merger\n',...
                    f_ams.chr_names{chr});
            end
        else
            f_ams.coord_per_position(chr) = 1;
            no_ti = nanmerge(no_ti,'sum');
            no_ci = nanmerge(no_ci,'sum');
        end
    end

    % substitute in {f_ams}
    f_ams.No_Ts{chr} = no_ti;
    f_ams.No_Cs{chr} = no_ci;
    f_ams.diagnostics.effective_coverage(chr) = ...
        (nansum(no_ti) + nansum(no_ci)) / length(no_ti);
    f_ams.No_As{chr} = [];
    f_ams.No_Gs{chr} = [];
    f_ams.aOga{chr} = [];
    f_ams.tOct{chr} = [];
    total_removed = total_removed + no_removed;
end
f_ams.is_filtered = true;
fprintf(1,'In total %s positions were removed\n',...
    numcommas(total_removed));

% close file
fclose(fid);
    
% #########################################################################
function [p_alg, to_merge, fname] = parseInput(ams,varargin)

% PARSEINPUT parses input line.
% ----------------------------------------------------
% [p_alg, to_merge, fname] = parseInput(ams, varargin)
% ----------------------------------------------------
% Description:  parses the input line.
% Input:        {ams} amSample object.
%               {varargin} original input line.
% Output:       {p_alg} structure with the fields 'max_coverage',
%                   'method','max_TsPerCoverage','max_aOga', and
%                   'max_No_As'.
%               {to_merge} true if it is desired to merge adjacent CpG
%                   positions.
%               {fname} log file name.

% defaults
is_tOct = false;
is_TsPerCov = false;
min_No_Ts = 1;
fname = sprintf('%s_filter.txt',ams.name);
to_merge = true;
if ~isempty(ams.p_filters.max_coverage)
    max_coverage = ams.p_filters.max_coverage;
else
    max_coverage = 100;
end
if strcmp(ams.library,'SS')
    method = 'both';
else
    method = 'tOct';
end
if ~isempty(ams.p_filters.max_TsPerCoverage)
    max_TsPerCoverage = ams.p_filters.max_TsPerCoverage;
else
    max_TsPerCoverage = 0.25;
end
p_alg = ams.p_filters;
p_alg.method = method;
p_alg.max_coverage = max_coverage;
p_alg.max_TsPerCoverage = max_TsPerCoverage;
p_alg.max_aOga = 0.25;
p_alg.max_No_As = 1;
p_alg.is_filtered = true;

% read user-specific instructions
for ii = 1:2:(nargin-2)
    switch str2keyword(varargin{ii},8)
        case 'method  '     % instruction: method
            switch str2keyword(varargin{ii+1},5)
                case 'max_t'
                    p_alg.method = 'max_tOct';
                case 'max_a'
                    p_alg.method = 'max_aOga';
                case 'both '
                    p_alg.method = 'both';
                otherwise
                    error(sprints(...
                        'Unknown value for ''method'': ''%s''',...
                        varargin{ii+1}));
            end
        case 'max_toct'     % instruction: max_tOct
            p_alg.max_TsPerCoverage = varargin{ii+1};
            is_tOct = true;
        case 'max_aoga'     % instruction: max_aOga
            p_alg.max_aOga = varargin{ii+1};
        case 'max_cove'     % instruction: max_coverage
            p_alg.max_coverage = varargin{ii+1};
        case 'max_no_a'     % instruction: max_No_As
            p_alg.max_No_As = varargin{ii+1};
        case 'min_no_t'    % instruction: min_No_Ts
            min_No_Ts = varargin{ii+1};
        case 'max_tspe'    % instruction: max_TsPerCoverage
            p_alg.max_TsPerCoverage = varargin{ii+1};
            is_TsPerCov = true;
        case 'merge   '    % instruction: merge
            to_merge = varargin{ii+1};
        case 'fname   '    % instruction: fname
            fname = varargin{ii+1};
    end
end

% 'max_tOct' and 'max_TsPerCoverage' cannot be used at the same time
if is_tOct && is_TsPerCov
    errstr = 'Both ''max_TsPerCoverage'' and ''max_toct''';
    error('%s were used in the input',errstr);
end

% bring input parameters into standard form - max_coverage
if isscalar(p_alg.max_coverage)
    p_alg.max_coverage = p_alg.max_coverage * ones(1,no_chr);
end

% bring input parameters into standard form - max_TsPerCoverage
if is_tOct            % max_tOct
    if ~iscell(p_alg.max_TsPerCoverage)          
        max_TsPerCoverage = cell(1,no_chr);
        if isscalar(p_alg.max_TsPerCoverage)    % a single ratio
            for chr = 1:no_chr
                max_TsPerCoverage{chr} = ceil( p_alg.max_TsPerCoverage ...
                    * (1:p_alg.max_coverage(chr))) - 1;
            end
        else                                    % ratio per chromosome
            for chr = 1:no_chr
                max_TsPerCoverage{chr} = ceil( ...
                    p_alg.max_TsPerCoverage(chr) ...
                    * (1:p_alg.max_coverage(chr))) - 1;
            end
        end
    else    % different ratio per coverage level in each chromosome
        for chr = 1:ams.no_chrs
            max_TsPerCoverage{chr} = ceil( p_alg.max_TsPerCoverage{chr} ...
                .* (1:p_alg.max_coverage(chr))) - 1;
        end
    end
    p_alg.max_TsPerCoverage = max_TsPerCoverage;
else            % max_TsPerCoverage or default
    if ~iscell(p_alg.max_TsPerCoverage)          
        max_TsPerCoverage = cell(1,no_chr);
        if isscalar(p_alg.max_TsPerCoverage)    % a single number
            for chr = 1:no_chr
                max_TsPerCoverage{chr} = p_alg.max_TsPerCoverage ...
                    * ones(1,p_alg.max_coverage(chr));
            end
        else                                    % number per chromosome
            for chr = 1:no_chr
                max_TsPerCoverage{chr} = p_alg.max_TsPerCoverage(chr) ...
                    * ones(1,p_alg.max_coverage(chr));
            end
        end
        p_alg.max_TsPerCoverage = max_TsPerCoverage;
    end
end

% bring input parameters into standard form - max_aOga
if isscalar(p_alg.max_aOga)
    p_alg.max_aOga = p_alg.max_aOga * ones(1,ams.no_chrs);
end

% bring input parameters into standard form - max_No_As
if isscalar(p_alg.max_No_As)
    p_alg.max_No_As = p_alg.max_No_As * ones(1,ams.no_chrs);
end

% bring input parameters into standard form - min_No_Ts
for chr = 1:ams.no_chrs
    vec = p_alg.max_TsPerCoverage{chr};
    vec(vec < min_No_Ts) = min_No_Ts;
    vec(1:min_No_Ts) = 1:min_No_Ts;
    p_alg.max_TsPerCoverage{chr} = vec;
end

% #########################################################################
function to_remove = extend_removed_positions(to_remove)

% EXTEND_REMOVED_POSITIONS removes both CpG positions
% -----------------------------------------------
% to_remove = extend_removed_positions(to_remove)
% -----------------------------------------------
% Description:  If a CpG position is marked to be removed, both consecutive
%               CpG positions are marked to be removed.
% Input:        {to_remove} original positions that are marked for removal.
% Output:       {to_remove} updated list.

odds = find(mod(to_remove,2));
evens = allbut(odds,length(to_remove));
to_remove = [to_remove; to_remove(odds)+1; to_remove(evens)-1];
to_remove = unique(to_remove);