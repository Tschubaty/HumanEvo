function No_Cs = getNo_Cs(obj,chr)

% GETNo_CS gets the No_Cs vector for a specific chromosome
% -------------------------
% No_Cs = obj.getNo_Cs(chr)
% -------------------------
% Description: gets the No_Cs vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {No_Cs} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 29-Aug-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

No_Cs = obj.No_Cs{chr};