function compareto(dm1,dm2)

% COMPARETO compares two DMR lists.
% ------------------
% compareto(dm1,dm2)
% ------------------
% Description:  compares two lists of DMRs.
% Input:        {dm1} DMRs object.
%               {dm2} another DMR object.

% © Liran Carmel
% Classification: Properties
% Last revision date: 18-Jan-2018

% get names
name1 = inputname(1);
name2 = inputname(2);

% sanity check
no_chr = length(dm1);
if length(dm2) ~= no_chr
    error('DMR objects have a different number of chromosomes');
end

% loop on chromosomes
for chr = 1:no_chr
    % define a vector of zeros whose length spans all DMRs
    max_coord = max([dm1(chr).CpG_end dm2(chr).CpG_end]);
    vec = zeros(1,max_coord);
    % mark all DMRs that belong to the first list
    for dm = 1:dm1(chr).no_DMRs
        vec(dm1(chr).CpG_begin(dm):dm1(chr).CpG_end(dm)) = 1;
    end
    % mark all DMRs that belong to the second list
    for dm = 1:dm2(chr).no_DMRs
        vec(dm2(chr).CpG_begin(dm):dm2(chr).CpG_end(dm)) = ...
            vec(dm2(chr).CpG_begin(dm):dm2(chr).CpG_end(dm)) + 2;
    end
    % compute total lengths and overlaps
    tot_length1 = sum(vec==1 | vec==3);
    tot_length2 = sum(vec==2 | vec==3);
    overlap = sum(vec==3);
    % print
    fprintf(1,'Chromosome %d:\n',chr);
    fprintf(1,'\tTotal length of DMRs in %s: %s\n',...
        name1,numcommas(tot_length1));
    fprintf(1,'\tTotal length of DMRs in %s: %s ',...
        name2,numcommas(tot_length2));
    fprintf(1,'(%.2f%% of %s)\n',...
        (tot_length2-tot_length1)/tot_length1*100,name1);
    fprintf(1,'\tLength of overlap: %s\n',numcommas(overlap));
    fprintf(1,'\tOverlap is %.2f%% from %s, and %.2f%% from %s\n',...
        overlap/tot_length1*100,name1,overlap/tot_length2*100,name2);
end