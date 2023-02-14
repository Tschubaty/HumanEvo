function idx = index(obj,chr_name)

% INDEX gets the index of chromosomes
% -------------------------
% idx = obj.index(chr_name)
% -------------------------
% Description:  gets the index of chromosomes.
% Input:        {chr_name} cell of name(s) of chromosome. If a single
%                   chromosome is required, it can be a string.
% Output:       {idx} index(ices) of the chromosome. If a chromosome name
%                   was not found, NaN is returned.

% (c) Liran Carmel
% Classification: Information extraction
% Last revision date: 27-Nov-2018

% force input to be cell
chr_name = cellstr(chr_name);

% find requested indices
no_chromosomes = length(chr_name);
idx = nan(1,no_chromosomes);
for ii = 1:no_chromosomes
    chr_idx = find(strcmp(chr_name{ii},obj.chromosomes));
    if ~isempty(chr_idx)
        idx(ii) = chr_idx;
    end
end