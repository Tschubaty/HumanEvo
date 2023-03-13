function methylation = getmethylation(obj,chr)

% GETMETHYLATION gets the methylation vector for a specific chromosome
% -------------------------------------
% methylation = obj.getmethylation(chr)
% -------------------------------------
% Description:  gets the methylation vector for a specific chromosome.
% Input:        {chr} index or name of chromosome.
% Output:       {methylation} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 01-Sep-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

methylation = obj.methylation{chr};