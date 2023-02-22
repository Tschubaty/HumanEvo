function dm = populateCpGcoords(dm,coordi)

% POPULATECPGCOORDS populates the CpG coordinates of the DMRs.
% ---------------------------------
% dm = populateCpGcoords(dm, coord)
% ---------------------------------
% Description:  populates the CpG coordinates of the DMRs.
% Input:        {dm} DMR object.
%               {coord} vector of CpG coordinates.
% Output:       {dm} updated DMR object.

% © Liran Carmel
% Classification: Properties
% Last revision date: 18-Jan-2018

% initialize
CpG_begin = zeros(1,dm.no_DMRs);
CpG_end = zeros(1,dm.no_DMRs);
% populate
for ii = 1:dm.no_DMRs
    CpG_begin(ii) = find(coordi==dm.gen_begin(ii));
    CpG_end(ii) = find(coordi==dm.gen_end(ii)-1);
end
% set in object
dm = set(dm,'CpG_begin',CpG_begin,'CpG_end',CpG_end);