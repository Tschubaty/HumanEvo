function obj = reconstructmethylation(obj,varargin)

% RECONSTRUCTMETHYLATION computes methylation from tOct data
% ------------------------------------------------------
% obj = obj.reconstructmethylation(i_name, i_value, ...)
% ------------------------------------------------------
% Description:  computes methylation from tOct data, based on the linear
%               transformation:
%                   meth = slope * no_ti / no_cti + intercept
% Input:        <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'win_size' (def: 'auto') window size for smoothing.
%                       If 'auto', a recommended value is computed for each
%                       chromosome. Otherwise, it can be a scalar (used for
%                       all chromosomes) or a vector over chromosomes.
%                   'winsize_algorithm' a structure with parameters
%                       required to determine window size, see structure
%                       {p_alg} in AMSAMPLE/DETERMINEWINDOWSIZE.
%                   'slope' (def: 1/obj.drate.rate.global) of the linear
%                       transformation. It can be a scalar (where it is
%                       used for all chromosomes) or a vector over
%                       chromosomes.
%                   'intercept' (def: 0) of the linear transformation,
%                       can be a scalar (where it is used for all
%                       chromosomes) or a vector over chromosomes.
%                   'lcf' (def: 0.05) low coverage factor.
%                   'to_report' (default: true) a binary variable,
%                       determining whether to report output to screen.
% Output:       {obj} AMSAMPLE object with a modified 'methylation' field.

% (c) Liran Carmel & David Gokhman & Yoav Mathov
% Classification: RoAM
% Last revision date: 15-Jan-2019

% parse input line
[p_alg, to_report] = parseInput(obj,varargin{:});

% initialize
no_chr = obj.no_chrs;
meth = cell(1,no_chr);

% loop on chromosomes
for chr = 1:no_chr
    % report
    if to_report
        if isempty(obj.chr_names)
            fprintf(1,'Computing methylation in chromosome #%d\n',chr);
        else
            fprintf(1,'Computing methylation in %s\n',obj.chr_names{chr});
        end
    end
    % get smoothed No_Ts and No_CTs
    [no_ti, no_cti] = smooth(obj,chr,p_alg.win_size(chr));
    % remove regions with particularly low coverage
    lct = findlowcoveragethreshold(no_cti,p_alg.lcf);
    no_cti(no_cti<lct) = nan;
    % compute methylation
    methi = p_alg.slope(chr) * no_ti ./ no_cti + p_alg.intercept(chr);
    % trim to the range [0,1], but keep NaNs
    idx_nan = isnan(methi);
    methi = min(max(methi,0),1);
    methi(idx_nan) = nan;
    if to_report
        fprintf(1,'\tAverage methylation: %.2f\n',nanmean(methi));
    end
    meth{chr} = methi;
end

% substitute in methylation field
obj.methylation = struct('methylation',{meth},'win_size',p_alg.win_size,...
    'slope',p_alg.slope,'intercept',p_alg.intercept,'lcf',p_alg.lcf);

% #########################################################################
function [p_alg, to_report] = parseInput(obj,varargin)

% PARSEINPUT parses input line.
% ------------------------------------------------
% [p_alg, to_report] = parseInput(no_chr,varargin)
% ------------------------------------------------
% Description:  parses the input line.
% Input:        {obj} AMSAMPLE object.
%               {varargin} original input line.
% Output:       {p_alg} structure with algorithm parameters.
%               {to_report} whether to display output on screen.

% defaults
no_chr = obj.no_chrs;
is_auto_window = true;
win_size = nan*ones(1,no_chr);
slope = 1/obj.drate.rate.global*ones(1,no_chr);
intercept = 0;
to_report = true;
lcf = 0.05;
winsize_algorithm = struct;

% read user-specific instructions
for ii = 1:2:(nargin-1)
    switch str2keyword(varargin{ii},5)
        case 'slope'    % instruction: slope
            slope = varargin{ii+1};
        case 'inter'    % instruction: intercept
            intercept = varargin{ii+1};
        case 'win_s'    % instruction: winsize
            if ischar(varargin{ii+1})
                is_auto_window = true;
            else
                is_auto_window = false;
                win_size = varargin{ii+1};
            end
        case 'lcf  '    % instruction: lcf
            lcf = varargin{ii+1};
        case 'winsi'    % instruction: winsize_algorithm
            fnames = fieldnames(varargin{ii+1});
            for ff = 1:length(fnames)
                winsize_algorithm.(fnames{ff}) = ...
                    varargin{ii+1}.(fnames{ff});
            end
        case 'to_re'    % instruction: to_report
            to_report = varargin{ii+1};
    end
end

% bring parameters into standard formt - win_size
if ~is_auto_window
    if length(win_size) == 1
        if ~mod(win_size,2)
            win_size = win_size + 1;
        end
        win_size = win_size * ones(1,no_chr);
    else
        for chr = 1:no_chr
            if ~mod(win_size(chr),2)
                win_size(chr) = win_size(chr) + 1;
            end
        end
    end
else        % if window size should be determined automatically
    win_size = zeros(1,no_chr);
    baseline_left = 'win_size(';
    baseline_middle = ')= obj.determinewinsize(';
    baseline_right = '';
    fnames = fieldnames(winsize_algorithm);
    for ii = 1:length(fnames)
        baseline_right = sprintf('%s,''%s'',%s',baseline_right,...
            fnames{ii},num2str(winsize_algorithm.(fnames{ii})));
    end
    baseline_right = sprintf('%s);',baseline_right);
    for chr = 1:no_chr
        eval(sprintf('%s%d%s%d%s',baseline_left,chr,baseline_middle,...
            chr,baseline_right));
    end
end

% bring parameters into standard formt - slope
if length(slope) == 1
    slope = slope * ones(1,no_chr);
end

% bring parameters into standard formt - intercept
if length(intercept) == 1
    intercept = intercept * ones(1,no_chr);
end

% aggregate all in one structure
p_alg = struct('win_size',win_size,'slope',slope,...
    'intercept',intercept,'lcf',lcf);