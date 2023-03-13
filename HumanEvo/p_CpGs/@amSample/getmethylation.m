function meth = getmethylation(obj,chr)

% GETMETHYLATION gets the methylation vector for a specific chromosome
% -------------------------------
% No_As = obj.getmethylation(chr)
% -------------------------------
% Description: gets the methylation vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {meth} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 08-Sep-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

meth = obj.methylation.methylation{chr};