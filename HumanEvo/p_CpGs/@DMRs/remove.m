function dm = remove(dm,chr,idx)

% REMOVE Removes DMRs from the current list
% -----------------------------------
% Description:  Removes specific DMRs from the list in a specific
%               chromosome.
% Input:        {dm} DMRs instance.
%               {chr} chromosome name.
%               {idx} indexes of the DMRs to remove.

% © Yoav Mathov & Liran Carmel
% Classification: SET/GET functions
% Last revision date: 27-Nov-2018

chr = dm.index(chr);
dm.cDMRs(chr) = cDMRs(chr).remove(idx);