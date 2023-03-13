function obj = estimatedrate(obj,varargin)

% ESTIMATEDRATE estimates deamination rate.
% -----------------------------------------------
% drate = obj.estimatedrate(i_name, i_value, ...)
% -----------------------------------------------
% Description:  estimates deamination rate.
% Input:        <{i_name, i_value}> pairs that modify the algorithm.
%                   Currently available:
%                   'method' (default: 'reference') the method used to
%                       perform the estimation. Can be either:
%                       'reference' the preferred and most accurate method,
%                           that should be used when we have a vector of 
%                           beta-values as a reference.
%                       'global' should be used in cases we have no
%                           measured reference, and the estimation is based
%                           on the pre-assumption of the total level of
%                           methylation in the genome.
%                       'estref' should be used in cases we have no
%                           measured reference, and the estimation is based
%                           on estimating positions with the beta-values
%                           close to one.
%                   'global_methylation' the estimated value of global
%                       genomic methylation, either as a fraction or as
%                       percentage. Applicable for 'method'='global'.
%                   'reference' mmSample, containing the beta-values of the
%                       reference. Applicable for 'method'='reference'.
%                   'min_beta' (def: 1) minimum beta value to take for
%                       the estimation. Applicable for
%                       'method'='reference'.
%                   'min_coverage' (def: 1) minimum coverage of sites
%                       that are used for the estimation.
% Output:       {obj} an amSample object with udpated 'drate' field. This
%                   field is a structure with several fields.
%                   'global' global deamination rate (computed from all
%                       chromosomes).
%                   'dglobal' STD of the global deamination rate (assuming
%                       binomial distribution).
%                   'local' local deamination rate, computed for each
%                       chromosome separately.
%                   'dlocal' STD of the local deamination rate (assuming
%                       binomial distribution).
%                   'no_position' in each chromosome, upon which
%                       deamination rate was computed.

% (c) Liran Carmel & David Gokhman
% Classification: RoAM
% Last revision date: 01-Sep-2018

% parse input line
[method, p_meth] = parseInput(varargin{:});

% sanity check
if ~obj.p_filters.is_filtered
    error('%s is not filtered',inputname(1));
end
if ~all(obj.coord_per_position == 1)
    error('%s is not merged',inputname(1));
end

% verify reference is merged and scaled
if strcmp(method,'reference')
    p_meth.ref = p_meth.ref.merge();
    p_meth.ref = p_meth.ref.scale();
end

% initialize
drate = struct('global',nan,'dglobal',nan,...
    'local',nan(1,obj.no_chrs),'dlocal',nan(1,obj.no_chrs),...
    'no_positions',nan(1,obj.no_chrs));
fact = 1;
if strcmp(method,'global')
    fact = 1 / p_meth.global_methylation;
end
totT = 0;
totCT = 0;

% Loop on chromosomes
for chr = 1:obj.no_chrs
    % report
    fprintf(1,'Estimating deamination rate in %s\n',obj.chr_names{chr});
    
    % vectors of current chromosome
    no_ti = obj.getNo_Ts(chr);
    no_cti = no_ti + obj.getNo_Cs(chr);
    
    % remove positions for which the reference has low beta-values
    if strcmp(method,'reference')
        % Take only positions where {ref}>=min_beta
        to_include = find(p_meth.ref.getmethylation(obj.chr_names{chr}) ...
            >= p_meth.min_beta);
        no_ti = no_ti(to_include);
        no_cti = no_cti(to_include);
    end
    
    % remove positions that are not covered high enough
    to_include = find(no_cti >= p_meth.min_coverage);
    no_ti = no_ti(to_include);
    no_cti = no_cti(to_include);
    
    % compute estimates based on the current chromosome
    drate.local(chr) = fact * nansum(no_ti) / nansum(no_cti);
    drate.no_positions(chr) = sum(isfinite(no_ti));
    drate.dlocal(chr) = fact * sqrt( drate.local(chr) * ...
        (1 - drate.local(chr)) / drate.no_positions(chr) );
    
    % accumulate sums
    totT = totT + nansum(no_ti);
    totCT = totCT + nansum(no_cti);
end

% evaluate global degradation rate
drate.global = fact * totT / totCT;
drate.dglobal = fact * sqrt( drate.global * (1 - drate.global) ...
    / sum(drate.no_positions) );

% substitute into field
switch method
    case 'reference'
        obj.drate = struct('method','reference','rate',drate,...
            'reference',p_meth.ref.name,'min_beta',p_meth.min_beta,...
            'min_coverage',p_meth.min_coverage);
    case 'global'
        obj.drate = struct('method','global','rate',drate,...
            'global_methylation',p_meth.global_methylation,...
            'min_coverage',p_meth.min_coverage);
    case 'estref'
        obj.drate = struct('method','estref','rate',drate);
end

% #########################################################################
function [method, meth_params] = parseInput(varargin)

% PARSEINPUT parses input line.
% --------------------------------------------
% [method, meth_params] = parseInput(varargin)
% --------------------------------------------
% Description:  parses the input line.
% Input:        {varargin} original input line.
% Output:       {method} is either 'reference', 'global', or 'estref'.
%               {meth_params} the parameters specific to the methodology.

% defaults
method = 'reference';
p_global = struct('global_methylation',nan,'min_coverage',1);
p_ref = struct('ref',[],'min_beta',1,'min_coverage',1);

% read user-specific instructions
for ii = 1:2:nargin
    switch str2keyword(varargin{ii},5)
        case 'metho'   % instruction: method
            switch str2keyword(varargin{ii+1},3)
                case 'ref'
                    method = 'reference';
                case 'glo'
                    method = 'global';
                case 'est'
                    method = 'estref';
                otherwise
                    error(sprints('Unknown method ''%s''',...
                        varargin{ii+1}));
            end
        case 'globa'   % instruction: global_methylation
            p_global.global_methylation = varargin{ii+1};
            if p_global.global_methylation > 1   % in case %% is provided
                p_global.global_methylation = 0.01 * ...
                    p_global.global_methylation;
            end
        case 'refer'   % instruction: reference
            p_ref.ref = varargin{ii+1};
        case 'min_b'   % instruction: min_beta
            p_ref.min_beta = varargin{ii+1};
            if p_ref.min_beta > 1
                p_ref.min_beta = 0.01*p_ref.min_beta;
            end
        case 'min_c'   % instruction: min_coverage
            p_ref.min_coverage = varargin{ii+1};
            p_global.min_coverage = varargin{ii+1};
    end
end

% test consistency
if strcmp(method,'global')
    meth_params = p_global;
    if isnan(meth_params.global_methylation)
        err_str = ['If using the ''global'' method, the parameter '...
            '''global_methylation'' should be provided.'];
        error(err_str);
    end
elseif strcmp(method,'reference')
    meth_params = p_ref;
    if isempty(meth_params.ref)
        err_str = ['If using the ''reference'' method, a reference '...
            'should be provided.'];
        error(err_str);
    end
end

% display parameters
fprintf(1,'Estimating deamination rate using the ''%s'' method:\n',method);
if strcmp(method,'global')
    fprintf(1,'\tgloabl_methylation: %.2f\n',...
        meth_params.global_methylation);
elseif strcmp(method,'reference')
    fprintf(1,'\tmin_beta: %.2f\n',meth_params.min_beta);
end
fprintf(1,'\tmin_coverage: %d\n',meth_params.min_coverage);