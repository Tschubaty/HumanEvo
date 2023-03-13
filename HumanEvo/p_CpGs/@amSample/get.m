function val = get(obj,p_name)

% GET get method
% ----------------------------------------
% property_value = get(obj, property_name)
% ----------------------------------------
% Description: gets field values.
% Input:       {obj} AMSAMPLE instance(s).
%              {property_name} legal property.
% Output:      {property_value} value of {property_name}. If {obj} is not a
%                   scalar, {property_value} is a cell of the dimensions of
%                   {obj}.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 26-Aug-2018

% initialize {val}
val = cell(size(obj));

% get field value for each element in {obj}
for ii = 1:numel(obj)
    val{ii} = obj(ii).(p_name);
end

% if {obj} is a scalar, don't return cell
if isscalar(obj)
    val = val{1};
end