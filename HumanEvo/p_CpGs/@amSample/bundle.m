function obj = bundle(obj,p_CpGs)

% BUNDLE converts p_CpG structure into amSample object
% ------------------------
% obj = obj.bundle(p_CpGs)
% ------------------------
% Description: converts p_CpG structure into amSample object.
% Input:       {p_CpGs} p_CpGs structure.
% Output:      {samp} aSample object.

% (c) Liran Carmel
% Classification: SET/GET functions
% Last revision date: 29-Aug-2018

% clean previous information
obj.name = 'unknown';
obj.abbrev = 'unk';
obj.species = 'unknown';
obj.metadata = [];
obj.chr_names = [];

% get number of chromosomes
no_chr = length(p_CpGs);

% initialize
coords_pp = zeros(1,no_chr);
No_As = cell(1,no_chr);
No_Cs = cell(1,no_chr);
No_Gs = cell(1,no_chr);
No_Ts = cell(1,no_chr);
aOga = cell(1,no_chr);
tOct = cell(1,no_chr);

% populate fields
for chr = 1:no_chr
    coords_pp(chr) = p_CpGs(chr).coordinates_per_position;
    No_As{chr} = p_CpGs(chr).No_As;
    No_Cs{chr} = p_CpGs(chr).No_Cs;
    No_Gs{chr} = p_CpGs(chr).No_Gs;
    No_Ts{chr} = p_CpGs(chr).No_Ts;
    aOga{chr} = p_CpGs(chr).aOga;
    tOct{chr} = p_CpGs(chr).tOct;
end

% populate class
obj.coord_per_position = coords_pp;
obj.No_As = No_As;
obj.No_Cs = No_Cs;
obj.No_Gs = No_Gs;
obj.No_Ts = No_Ts;
obj.aOga = aOga;
obj.tOct = tOct;