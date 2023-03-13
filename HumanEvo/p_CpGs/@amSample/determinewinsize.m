function [winsize, p_alg] = determinewinsize(obj,varargin)

% DETERMINEWINSIZE estimates the window size.
% ------------------------------------------------------------------
% [winsize, p_alg] = obj.determinewinsize(chr, i_name, i_value, ...)
% ------------------------------------------------------------------
% Description:  estimates the optimal window size in cases that collecting
%               data from a window is required.
% Input:        {chr} index of chromosome.
%               <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'method' (def: 'prob') either 'prob' (for probability)
%                       or 'relerror' (for relative error).
%                   'coverage' (def: EFFECTIVECOVERAGE(obj,chr)) effective
%                       coverage of the chromosome.
%                   'drate' (def: obj.drate.rate.global) demaination rate
%                       of the sample.
%                   'min_meth' (def: 0.2) minimum methylation level we want
%                       to detect.
%                   'max_width' (def: 31) maximum allowed width. Computed
%                       window is not allowed to be larger than
%                       'max_width'.
%                   'p0' (def: 1e-2) applicable if 'method'='prob'. It is
%                       the probability to get zero counts in the window if
%                       the methylation is {min_meth}.
%                   '1oK' (def: 1/2.5) applicable if 'method'='relerror'.
%                       It is one over the maximum relative error in
%                       estimating the methylation in a window whose true
%                       methylation is {min_meth}. It also means that the
%                       mean is far {k} standard deviations from zero.
% Output:       {winsize} recommended window size, forced to be an odd
%                   number.
%               {p_alg} parameters of the algorithm, a structure with
%                   fields:
%                   'method' either 'prob' or 'relerror'
%                   'coverage' effective coverage
%                   'drate' deamination rate
%                   'min_meth' minimum methylation we want to detect
%                   'max_width' maximum allowed window size
%                   'parameter' specific parameters of the algorithm. This
%                       can be {p0} if 'method'='prob', and {Kinv} if
%                       'method' = 'relerror'

% (c) Liran Carmel & David Gokhman
% Classification: RoAM
% Last revision date: 20-Jan-2019

% parse input
p_alg = parseInput(obj,varargin{:});

% compute the window size
p = p_alg.min_meth * p_alg.drate;
prm = p_alg.parameter;
switch p_alg.method
    case 'prob'
        winsize = ceil( log(prm) / log(1 - p) / p_alg.coverage );
    case 'relerror'
        winsize = ceil((1-p) / p_alg.coverage / p / prm^2);
end

% narrow window if too wide
winsize = min(winsize,p_alg.max_width);

% force winow to be of odd length
if ~mod(winsize,2)
    winsize = winsize + 1;
end

% #########################################################################
function p_alg = parseInput(obj,varargin)

% PARSEINPUT parses input line.
% ---------------------------------
% p_alg = parseInput(obj, varargin)
% ---------------------------------
% Description:  parses the input line.
% Input:        {obj} AMSAMPLE object.
%               {varargin} original input line.
% Output:       {p_alg} structure with all parameters.

% first argument is fixed
chr = varargin{1};

% defaults
method = 'prob';
coverage = obj.effectivecoverage(chr);
drate = obj.drate.rate.global;
min_meth = 0.2;
p0 = 1e-2;
Kinv = 1/2.5;
max_width = 31;

% read user-specific instructions
for ii = 2:2:(nargin-1)
    switch str2keyword(varargin{ii},5)
        case 'metho'    % instruction: method
            switch str2keyword(varargin{ii+1},3)
                case 'pro'
                    method = 'prob';
                case 'rel'
                    method = 'relerror';
                otherwise
                    error('%s: unfamiliar method',varargin{ii+1});
            end
        case 'cover'    % instruction: coverage
            coverage = varargin{ii+1};
        case 'drate'    % instruction: drate
            drate = varargin{ii+1};
        case 'min_m'    % instruction: min_meth
            min_meth = varargin{ii+1};
            if min_meth > 1
                min_meth = 0.01 * min_meth;
            end
        case 'p0   '    % instruction: p0
            p0 = varargin{ii+1};
        case '1ok  '    % instruction: 1oK
            Kinv = varargin{ii+1};
        case 'max_w'    % instruction: max_width
            max_width = varargin{ii+1};
        otherwise
            error('%s: unfamiliar instruction',varargin{ii});
    end
end

% aggregate parameters in a structure
p_alg = struct('method','prob','coverage',coverage,'drate',drate,...
        'min_meth',min_meth,'max_width',max_width,'parameter',[]);
if strcmp(method,'prob')
    p_alg.parameter = p0;
else
    p_alg.parameter = Kinv;
end