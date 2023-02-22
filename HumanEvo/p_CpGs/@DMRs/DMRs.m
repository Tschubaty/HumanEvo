classdef DMRs

% DMRs constructor method
% --------------------------------------------
% (1) obj = DMRs()
% (2) obj = DMRs(size)
% (3) obj = DMRs(obj0)
% (4) obj = DMRs(field_name, field_value, ...)
% (5) DMRs('help')
% --------------------------------------------
% Description:  constructs DMRs instance(s).
% Input:        (1) An empty default DMRs is formed, identified by having
%                   no DMRs.
%               (2) {size} if scalar, {obj} would be an empty DMRs
%                  vector of length {size}. If a vector, {obj} would be
%                  an empty DMRs n-dimensional matrix of
%                  dimensions dictated by {size}.
%               (3) {obj0} DMRs instance. In this case {obj} would be
%                  an identical copy.
%               (4) pairs of field names accompanied by their
%                  corresponding values.
%               (5) opens a web page with help.
% Output:       {obj} DMRs instance(s).

% ? Liran Carmel & Yoav Mathov & David Gokhman
% Classification: Classdef
% Last revision date: 20-Mar-2018

    properties (SetAccess = public)
        samples = '';       % name of samples
        groups = '';        % links samples to groups (GROUPING)
        species = '';       % species of samples
        reference = '';     % genomic reference of the samples
    end
    
    properties (SetAccess = protected)
        chromosomes = {};   % list of chromosome / genomic section names
        cDMRs = [];         % list of chromosome-specific DMRs (cDMRs)
        is_ancient = [];    % TRUE for samples that are ancient
        algorithm = [];     % DMR-detection algorithm
    end
    
    properties (SetAccess = protected, Dependent = true)
        no_samples          % number of samples (s)
        no_chromosomes      % number of chromosomes (c)
    end
    
    methods
        
        % constructor
        function obj = DMRs(varargin)
            % decide on which kind of constructor should be used
            if nargin > 0
                switch class(varargin{1})
                    case 'DMRs'         % case (3)
                        obj = varargin{1};
                    case 'char'
                        if nargin == 1  % case (5)
                            web(['http://carmelab.huji.ac.il/' ...
                                'software/MVA/DMRs_obj.html']);
                        else            % case (4)
                            % initialize a default instance
                            obj = DMRs;
                            % loop on all variables
                            for ii = 1:2:nargin 
                                obj.(varargin{ii}) = varargin{ii+1};
                            end
                        end
                    case 'double'           % case (2)
                        obj = [];
                        for ii = 1:prod(varargin{1})
                            obj = [obj DMRs];   %#ok<AGROW>
                        end
                        if length(varargin{1}) > 1
                            obj = reshape(obj,varargin{1});
                        end
                end
            end
        end
        
        % modifications to get
        function no_DMRs = get.no_chromosomes(obj)
            no_DMRs = length(obj.chromosomes);
        end
        % modifications to get
        function no_samples = get.no_samples(obj)
            no_samples = length(obj.samples);
        end
        
    end
end