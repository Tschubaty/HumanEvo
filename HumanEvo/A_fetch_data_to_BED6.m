
% v 01 18.04.2023
% read methylation files
%
%% include github repo
git_folder = 'C:\Users\Daniel Batyrev\Documents\GitHub\RoAM';
addpath(genpath(git_folder));

%% consrtans and paramter

INPUT_FOLDER = "p_CpGs";
OUTPUT_FOLDER = "methylation_data";

%% load genomic coordinates
load(fullfile(INPUT_FOLDER,"gc_CpGs.mat"));

%%
reg_expression =  fullfile(INPUT_FOLDER,"*.mat");
files = dir(reg_expression);
pat = "f_" + digitsPattern(4) + ".mat";


samples = {};


for i =1:numel(files)
    if matches(files(i).name,pat)
        fprintf("start loading %s \n",files(i).name)
        samples{end +1} = files(i).name;
        load(fullfile(INPUT_FOLDER,files(i).name));
        sample_name = strcat("I",  erase(erase(files(i).name,"f_"),".mat"));
        mkdir(fullfile(OUTPUT_FOLDER, sample_name));

        for chromosome_number = 1:22

            bed_filename = fullfile(OUTPUT_FOLDER, sample_name,...
                strcat(sample_name,".",gc_CpGs.chr_names{chromosome_number},".bed"));
            if ~isfile(bed_filename)
                % chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random)
                % chromStart - The starting position of the feature in the chromosome or scaffold.
                %               The first base in a chromosome is numbered 0.
                % chromEnd -The ending position of the feature in the chromosome or scaffold.
                %               The chromEnd base is not included in the display of the feature,
                %               however, the number in position format will be represented.
                %               For example, the first 100 bases of chromosome 1 are defined as
                %               chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100),
                %               but will represent the position notation chr1:1-100.
                % name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window
                %           when the track is open to full display mode or directly to the left of the item in pack mode.
                % score - A score between 0 and 1.
                % strand - Defines the strand. Either "." (=no strand) or "+" or "-".
                varTypes = ["string","int32","int32","string","double","string"];
                varNames = ["chrom","start","end","name","score","strand"];
                bed_length = numel(gc_CpGs.coordinates{chromosome_number});
                sz = [bed_length numel(varTypes)];


                df = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

                df.chrom(:) = gc_CpGs.chr_names{chromosome_number};
                df.start(:) = gc_CpGs.coordinates{chromosome_number};
                df.end(:) = gc_CpGs.coordinates{chromosome_number};
                df.name(:) = sample_name;
                df.score(:) = fas.getmethylation(chromosome_number);
                df.strand(:) = "+";

                writetable(df, bed_filename, ...
                    'Delimiter','\t','FileType',"text");
            end
            fprintf("chromosome %d saved \n",chromosome_number);
        end
    end
end


disp(numel(samples))