function obj = set(obj,varargin)

% SET set method
% -------------------------------------------------
% obj = set(obj, property_name, property_value,...)
% -------------------------------------------------
% Description: sets field values. If {obj} is not a scalar, the same
%              property values are substituted for all instances.
% Input:       {obj} AMSAMPLE instance.
%              {property_name},{property_value} legal pairs.
% Output:      {obj} updated AMSAMPLE instance.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 26-Aug-2018

% use the set methods of the individual parameters
for jj = 1:numel(obj)
    for ii = 1:2:length(varargin)
        obj(jj).(varargin{ii}) = varargin{ii+1};
    end
end