function No_Gs = getNo_Gs(obj,chr)

% GETNO_GS gets the No_Gs vector for a specific chromosome
% -------------------------
% No_Gs = obj.getNo_Gs(chr)
% -------------------------
% Description: gets the No_Gs vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {No_Gs} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 29-Aug-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

No_Gs = obj.No_Gs{chr};