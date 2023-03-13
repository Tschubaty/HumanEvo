function aOga = getaOga(obj,chr)

% GETAOGA gets the aOga vector for a specific chromosome
% -------------------------
% aOga = obj.getaOga(chr)
% -------------------------
% Description: gets the aOga vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {aOga} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 29-Aug-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

aOga = obj.aOga{chr};