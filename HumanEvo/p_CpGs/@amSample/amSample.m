classdef amSample

% AMSAMPLE constructor method
% -----------------------------------------------
% (1) obj = amSample()
% (2) obj = amSample(size)
% (3) obj = amSample(obj0)
% (4) obj = amSample(field_name, field_value, ...)
% (5) amSample('help')
% -----------------------------------------------
% Description:  constructs aSample instance(s).
% Input:        (1) An empty default AMSAMPLE is formed, identified by
%                   having sample data.
%               (2) {size} if scalar, {obj} would be an empty AMSAMPLE
%                  vector of length {size}. If a vector, {obj} would be
%                  an empty AMSAMPLE n-dimensional matrix of
%                  dimensions dictated by {size}.
%               (3) {obj0} AMSAMPLE instance. In this case {obj} would be
%                  an identical copy.
%               (4) pairs of field names accompanied by their
%                  corresponding values.
%               (5) opens a web page with help.
% Output:       {obj} AMSAMPLE instance(s).

% (c) Liran Carmel
% Classification: Classdef
% Last revision date: 28-Nov-2018

    properties (SetAccess = public)
        name = 'unknown';         % full name of the sample
        abbrev = 'unk';           % abbreviated name
        species = 'unknown';      % species of sample
        reference = '';           % genome reference
        library = '';             % either 'single' or 'double'
    end
    
    properties (SetAccess = protected)
        chr_names = [];     % names of chromosomes. Cell of strings (1xC)
        coord_per_position = []; % vector (1xC) with 1 or 2
        No_As = [];         % Cell of vectors (1xC)
        No_Cs = [];         % Cell of vectors (1xC)
        No_Gs = [];         % Cell of vectors (1xC)
        No_Ts = [];         % Cell of vectors (1xC)
        aOga = [];          % Cell of vectors (1xC)
        tOct = [];          % Cell of vectors (1xC)
        diagnostics = [];   % diagnostics structure
        p_filters = [];     % filtration parameters structure.
        is_filtered = false;    % TRUE if filtered
        is_simulated = false;   % TRUE if simulated
        methylation = [];   % methylation level
        drate = [];         % deamination rate
        metadata = [];      % structure with sample metadata (like sex)
    end
    
    properties (SetAccess = protected, Dependent = true)
        no_chrs             % number of chromosomes
    end
    
    methods
        
        % constructor
        function obj = amSample(varargin)
            % decide on which kind of constructor should be used
            if nargin > 0
                switch class(varargin{1})
                    case 'amSample'         % case (3)
                        obj = varargin{1};
                    case 'char'
                        if nargin == 1  % case (5)
                            web(['http://carmelab.huji.ac.il/' ...
                                'software/RoAM/amSample_obj.html']);
                        else            % case (4)
                            % initialize a default instance
                            obj = amSample;
                            % loop on all variables
                            for ii = 1:2:nargin 
                                obj.(varargin{ii}) = varargin{ii+1};
                            end
                        end
                    case 'double'           % case (2)
                        obj = [];
                        for ii = 1:prod(varargin{1})
                            obj = [obj amSample];   %#ok<AGROW>
                        end
                        if length(varargin{1}) > 1
                            obj = reshape(obj,varargin{1});
                        end
                end
            end
        end
        
        % modifications to get
        function no_chrs = get.no_chrs(obj)
            no_chrs = length(obj.coord_per_position);
        end
        
    end
end