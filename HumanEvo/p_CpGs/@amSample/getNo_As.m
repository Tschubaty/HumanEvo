function No_As = getNo_As(obj,chr)

% GETNO_AS gets the No_As vector for a specific chromosome
% -------------------------
% No_As = obj.getNo_As(chr)
% -------------------------
% Description: gets the No_As vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {No_As} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 29-Aug-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

No_As = obj.No_As{chr};