function bool = isfiltered(obj)

% ISFILTERED reports whether AMSAMPLE is filtered.
% -----------------------
% bool = obj.isfiltered()
% -----------------------
% Description: TRUE/FALSE indication of whether AMSAMPLE is filtered.
% Input:       {obj} AMSAMPLE object.
% Output:      {bool} boolean indicator.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 08-Sep-2018

bool = obj.p_filters.is_filtered;