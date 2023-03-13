function meth = gregionmethylation(obj,reg,gc)

% GREGIONMETHYLATION computes methylation in a region
% ---------------------------------
% obj = obj.gregionmethylation(reg)
% ---------------------------------
% Description:  Uses the same algorithm to reconstruct methylation as is
%               used in RECONSTRUCTMETHYLATION
%                   meth = slope * no_ti / no_cti + intercept
%               to compute the average methylation in a specific region
%               given by genomic coordinates.
% Input:        {reg} any region format.
%               {gc} GCOORIDNATES object of CpG positions.
% Output:       {meth} methylation value in the region.

% (c) Liran Carmel
% Classification: RoAM
% Last revision date: 10-Dec-2018

% bring region into a standard format
reg = standardregion(reg);

% find beginning and end in CpG positions
chr = gc.index(reg.chromosome);
CpG_beg = find(gc.coordinates{chr}>=reg.beg,1);
CpG_end = find(gc.coordinates{chr}<=reg.end,'last');

% find chromosome in {obj}
chr = obj.indexofchr(reg.chromosome);
 
% get smoothed No_Ts and No_CTs
no_ti = nansum(obj.No_Ts{chr}(CpG_beg:CpG_end));
no_cti = no_ti + nansum(obj.No_Cs{chr}(CpG_beg:CpG_end));
% compute methylation
meth = obj.methylation.slope(chr) * no_ti ./ no_cti + ...
    obj.methylation.intercept(chr);
% trim to the range [0,1], but keep NaNs
meth = min(max(meth,0),1);