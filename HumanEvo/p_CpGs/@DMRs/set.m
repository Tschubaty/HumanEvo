function obj = set(obj,varargin)

% SET set method
% -------------------------------------------------
% obj = set(obj, property_name, property_value,...)
% -------------------------------------------------
% Description: sets field values. If {obj} is not a scalar, the same
%              property values are substituted for all instances.
% Input:       {obj} DMRs instance.
%              {property_name},{property_value} legal pairs.
% Output:      {obj} updated DMRs instance.

% © Liran Carmel
% Classification: SET/GET functions
% Last revision date: 07-Jun-2017

% use the set methods of the individual parameters
for jj = 1:numel(obj)
    for ii = 1:2:length(varargin)
        obj(jj).(varargin{ii}) = varargin{ii+1};
    end
end