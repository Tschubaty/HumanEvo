function [no_DMRs, tot_DMRs] = noDMRs(dm)

% NODMRS reports the number of DMRs in each chromosome.
% --------------------------------
% [no_DMRs, tot_DMRs] = noDMRs(dm)
% --------------------------------
% Description:  NODMRS reports the number of DMRs in the object.
% Input:        {dm} cDMRs instance.
% Output:       {no_DMRs} number of DMRs per chromosome.
%               {tot_DMRs} total number of DMRs.

% © Liran Carmel
% Classification: SET/GET functions
% Last revision date: 30-Nov-2018

% loop on cDMRs
no_DMRs = nan(1,dm.no_chromosomes);
for ii = 1:dm.no_chromosomes
    no_DMRs(ii) = dm.cDMRs(ii).no_DMRs;
end
tot_DMRs = sum(no_DMRs);