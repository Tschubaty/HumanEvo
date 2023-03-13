function cvrg = effectivecoverage(obj,chr)

% EFFECTIVECOVERAGE gets the effective coverage of a chromosome
% ---------------------------------
% cvrg = obj.effectivecoverage(chr)
% ---------------------------------
% Description: gets the effective coverage of a chromosome.
% Input:       {chr} index or name of chromosome.
% Output:      {cvrg} effective coverage.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 08-Sep-2018

if ischar(chr)
    chr = obj.indexofchr(chr);
end

cvrg = obj.diagnostics.effective_coverage(chr);