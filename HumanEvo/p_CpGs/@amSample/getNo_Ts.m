function No_Ts = getNo_Ts(obj,chr)

% GETNo_TS gets the No_Ts vector for a specific chromosome
% -------------------------
% No_Ts = obj.getNo_Ts(chr)
% -------------------------
% Description: gets the No_Ts vector for a specific chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {No_Ts} the corresponding vector.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 26-Aug-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

No_Ts = obj.No_Ts{chr};